"""
Simulation orchestration for Chion.

`Simulation` owns reference and current model state (`ref` and `now`), while
`SimulationIntegrator` owns initialized runtime state for explicit stepping.
"""

"""
    Simulation(model; forcing, state=nothing, options=RunOptions(), ...)

Couple a model configuration with forcing, reference state, current state, and
execution/output options.
"""
mutable struct Simulation{M, R, S}
    model::M
    forcing::SnowpackForcing
    ref::R
    now::S
    options::RunOptions
end

function Simulation(
    model;
    forcing::SnowpackForcing,
    state=nothing,
    options::RunOptions=RunOptions(years=1, write_netcdf=false, netcdf_variables=Symbol[]),
    years::Union{Nothing, Integer}=nothing,
    backend=nothing,
    history_year_stride::Union{Nothing, Integer}=nothing,
    netcdf_variables=nothing,
    output_dir::Union{Nothing, AbstractString}=nothing,
    netcdf_path::Union{Nothing, AbstractString}=nothing,
    write_netcdf::Union{Nothing, Bool}=nothing,
    name::Union{Nothing, AbstractString}=nothing,
    input_label::Union{Nothing, AbstractString}=nothing,
)
    resolved_netcdf_variables = isnothing(netcdf_variables) ? options.netcdf_variables : normalize_netcdf_variables(netcdf_variables)
    resolved_write_netcdf = isnothing(write_netcdf) ? options.write_netcdf : write_netcdf
    run_options = RunOptions(
        name=isnothing(name) ? options.name : name,
        input_label=isnothing(input_label) ? options.input_label : input_label,
        output_dir=isnothing(output_dir) ? options.output_dir : output_dir,
        netcdf_path=isnothing(netcdf_path) ? options.netcdf_path : netcdf_path,
        write_netcdf=resolved_write_netcdf,
        netcdf_variables=resolved_netcdf_variables,
        years=isnothing(years) ? options.years : years,
        backend=isnothing(backend) ? options.backend : backend,
        history_year_stride=isnothing(history_year_stride) ? options.history_year_stride : history_year_stride,
    )
    now_state = isnothing(state) ? initial_state(model) : state
    return Simulation(model, forcing, reference_state(now_state), now_state, run_options)
end

function _normalize_model_name(model)
    name = lowercase(strip(String(model)))
    name in ("bessi", "bessimodel") && return :bessi
    name in ("pdd", "pddmodel") && return :pdd
    error("Unsupported model '$model'. Use `:bessi` or `:pdd`.")
end

function build_model(model, grid::SnowpackGrid; kwargs...)
    name = _normalize_model_name(model)
    name == :bessi && return BESSIModel(grid; kwargs...)
    name == :pdd && return PDDModel(grid; kwargs...)
    error("Unsupported model '$model'.")
end

function Simulation(
    model,
    grid::SnowpackGrid;
    model_kwargs=NamedTuple(),
    kwargs...,
)
    built_model = build_model(model, grid; model_kwargs...)
    return Simulation(built_model; kwargs...)
end

function Base.show(io::IO, ::MIME"text/plain", sim::Simulation)
    println(io, "Simulation")
    println(io, "  model: ", typeof(sim.model))
    println(io, "  ref: ", typeof(sim.ref))
    println(io, "  now: ", typeof(sim.now))
    println(io, "  columns: ", ncols(sim.model.grid))
    println(io, "  forcing steps: ", length(sim.forcing.time_values))
    println(io, "  backend: ", sim.options.backend)
    println(io, "  years: ", sim.options.years)
    println(io, "  netcdf vars: ", isempty(sim.options.netcdf_variables) ? "(none)" : join(string.(sim.options.netcdf_variables), ", "))
end

"""
    init_integrator(sim; io=stdout)

Initialize a runtime that can be advanced with `step!`, `run!`, and
`finalize!`.
"""
function init_integrator(
    sim::Simulation;
    io::IO=stdout,
)
    timings = StepTimingStats()
    validate_integrator_setup!(sim, sim.options)
    model_runtime = init_model_runtime!(sim, sim.options, timings)
    return _new_integrator(sim, io, timings, model_runtime)
end

finished(integrator::SimulationIntegrator) = _finished(integrator)

step!(integrator::SimulationIntegrator) = _step_scheduled!(integrator)

step!(integrator::SimulationIntegrator, n::Integer) = _step_n!(integrator, n)

step!(integrator::SimulationIntegrator, Δt_days::Real, force_dt::Bool=true) =
    _step_external!(integrator, Δt_days, force_dt)

_model_smb_ice_vector(::BESSIModel, ::BESSIState, runtime) =
    _host_vector(runtime.state.smb_ice; copy_array=true)

_model_smb_ice_vector(::PDDModel, ::PDDState, runtime) =
    _host_vector(runtime.smb_ice; copy_array=true)

function _yearly_grid(values::Vector{Float64}, grid)
    any(isnothing, (grid.x, grid.y, grid.js, grid.is, grid.mask)) &&
        return Matrix{Float64}(undef, 0, 0)
    return scatter_to_grid(values, grid.js, grid.is, size(grid.mask))
end

function _mask_inactive_yearly_outputs!(
    ice_sheet_net_forcing_yearly::Vector{Float64},
    mean_T_srf_K::Vector{Float64},
    active::Vector{Bool},
)
    length(active) == length(ice_sheet_net_forcing_yearly) || error("Active mask length must match yearly output length.")
    @inbounds for idx in eachindex(active)
        if !active[idx]
            ice_sheet_net_forcing_yearly[idx] = 0.0
            mean_T_srf_K[idx] = NaN
        end
    end
    return nothing
end

@inline _uses_monthly_output(options::RunOptions) =
    :monthly in options.netcdf_variables ||
    any(key -> key in MONTHLY_OUTPUT_VARS && !(key in DEFAULT_STATE_OUTPUT_VARS), options.netcdf_variables)

function _bessi_output_from_options(sim::Simulation, options::RunOptions)
    options.write_netcdf || return nothing
    monthly_mode = _uses_monthly_output(options)
    monthly_mode && length(options.netcdf_variables) > 1 &&
        error("`monthly` NetCDF output cannot currently be combined with other selectors.")
    vars = monthly_mode && :monthly in options.netcdf_variables ? MONTHLY_OUTPUT_VARS :
        monthly_mode ? intersect(options.netcdf_variables, MONTHLY_OUTPUT_VARS) :
        state_output_vars(options.netcdf_variables)
    isempty(vars) && return nothing
    return init_state_netcdf(
        resolve_netcdf_path(options),
        options,
        sim.forcing.time_values,
        sim.model.grid,
        sim.now,
        vars,
        ntime=monthly_mode ? options.years * 12 : options.years * length(sim.forcing.time_values),
        nlayer=sim.now.Ntot,
    )
end

function _bessi_step_kwargs(model::BESSIModel)
    return (
        diurnal_shortwave_substeps=model.diurnal_shortwave_substeps,
        diurnal_shortwave_threshold=model.diurnal_shortwave_threshold,
        diurnal_shortwave_max_substeps=model.diurnal_shortwave_max_substeps,
        diurnal_shortwave_min_air_temperature=model.diurnal_shortwave_min_air_temperature,
        diurnal_temperature_cycle=model.diurnal_temperature_cycle,
        diurnal_temperature_amplitude=model.diurnal_temperature_amplitude,
    )
end

function _sync_bessi_state!(sim::Simulation{<:BESSIModel}, backend_state, is_gpu::Bool, timings::StepTimingStats)
    if is_gpu
        time_block!(timings, :gpu_transfer) do
            _copy_bessi_state!(sim.now, cpu_state(backend_state))
        end
    end
    update_diagnostics!(sim.now)
    return nothing
end

@inline function _is_month_boundary(time_values::Vector{DateTime}, k::Int)
    k == length(time_values) && return true
    return month(time_values[k + 1]) != month(time_values[k])
end

function run!(
    sim::Simulation{<:BESSIModel};
    io::IO=stdout,
)
    integrator = init_integrator(sim; io=io)
    run!(integrator)
    return finalize!(integrator)
end

@inline _can_block_bessi_scheduled_steps(monthly_mode::Bool, nc) =
    !monthly_mode && nc === nothing

@inline function _bessi_scheduled_step_stop(integrator::SimulationIntegrator, nsteps::Int)
    runtime = integrator.model_runtime.data
    block_steps = _step_time_block_steps(_ka_backend(runtime.state.mass))
    return min(nsteps, integrator.time_index + block_steps - 1)
end

function _run_bessi_integrator!(integrator::SimulationIntegrator)
    run_options = integrator.sim.options

    sim = integrator.sim
    timings = integrator.timings
    runtime = integrator.model_runtime.data
    backend_state = runtime.state
    is_gpu = runtime.is_gpu
    nc = _bessi_output_from_options(sim, run_options)
    monthly_mode = nc !== nothing && _uses_monthly_output(run_options)
    monthly_state = monthly_mode ? MonthlyState(backend_state) : nothing
    monthly_output = monthly_mode ? MonthlyOutputBuffer(backend_state; nmonth=run_options.years * 12) : nothing
    record_index = 0
    nsteps = length(sim.forcing.time_values)
    year_summary = _summary_buffers(backend_state, (:thickness, :wet_mass, :bulk_density, :base_mass))
    time_block!(timings, :summarize_columns_initial) do
        summarize_year_state!(year_summary.thickness, year_summary.wet_mass, year_summary.bulk_density, year_summary.base_mass, backend_state)
    end
    prev_year_summary = map(copy, year_summary)
    delta_thickness = similar(year_summary.thickness)
    delta_wet_mass = similar(year_summary.wet_mass)
    delta_base_mass = similar(year_summary.base_mass)
    history = NamedTuple[]
    progress = Progress(run_options.years; desc="Running years: ", output=integrator.io, showspeed=true)
    while !_finished(integrator)
        year = integrator.completed_years + 1
        allow_step_blocking = _can_block_bessi_scheduled_steps(monthly_mode, nc)
        while integrator.completed_years < year
            if allow_step_blocking
                first_step = integrator.time_index
                last_step = _bessi_scheduled_step_stop(integrator, nsteps)
                if first_step == last_step
                    _step_scheduled!(integrator)
                else
                    _step_scheduled_range!(integrator, first_step:last_step)
                end
            else
                k = integrator.time_index
                _step_scheduled!(integrator)
                if monthly_mode
                    accumulate_monthly!(monthly_state, backend_state, sim.forcing.dt_days[Int(k)])
                    if _is_month_boundary(sim.forcing.time_values, Int(k))
                        time_block!(timings, :write_netcdf) do
                            finalize_monthly!(monthly_state, backend_state)
                            store_monthly!(monthly_output, monthly_state)
                        end
                        reset_monthly!(monthly_state)
                    end
                elseif nc !== nothing
                    _sync_bessi_state!(sim, backend_state, is_gpu, timings)
                    record_index += 1
                    time_block!(timings, :write_netcdf) do
                        write_nc!(nc, sim.now, record_index, sim.model.grid)
                    end
                end
            end
        end
        time_block!(timings, :summarize_columns_year) do
            summarize_year_state!(year_summary.thickness, year_summary.wet_mass, year_summary.bulk_density, year_summary.base_mass, backend_state)
        end
        record = time_block!(timings, :year_metrics) do
            make_year_record_and_deltas!(
                year,
                delta_thickness,
                delta_wet_mass,
                delta_base_mass,
                year_summary.thickness,
                year_summary.wet_mass,
                year_summary.bulk_density,
                year_summary.base_mass,
                prev_year_summary.thickness,
                prev_year_summary.wet_mass,
                prev_year_summary.base_mass,
            )
        end
        if should_record_year_metrics(year, run_options.years, run_options.history_year_stride)
            push!(history, record)
            time_block!(timings, :year_logging) do
                println(integrator.io, year_log_line(record))
                flush(integrator.io)
            end
        end
        copyto!(prev_year_summary.thickness, year_summary.thickness)
        copyto!(prev_year_summary.wet_mass, year_summary.wet_mass)
        copyto!(prev_year_summary.base_mass, year_summary.base_mass)
        next!(progress)
    end
    if monthly_mode
        time_block!(timings, :write_netcdf) do
            write_monthly_output_nc!(nc, monthly_output, 1, sim.model.grid)
        end
        record_index = monthly_output.count
    end
    _sync_bessi_state!(sim, backend_state, is_gpu, timings)
    status = :complete
    nc_path = nc === nothing ? "" : resolve_netcdf_path(run_options)
    nc !== nothing && close_output!(nc, status, record_index)

    integrator.history = history
    integrator.netcdf_path = nc_path

    run_wall_sec = (time_ns() - integrator.wall_t0) * 1.0e-9
    print_run_report(
        integrator.io,
        run_options,
        sim.forcing.time_values,
        history,
        status,
        run_wall_sec,
        timings;
        nc_path=nc_path,
    )
    return nothing
end

"""
    yearly_step!(integrator)

Advance a scheduled BESSI simulation by one full forcing year and return annual
coupling fields. `ice_sheet_net_forcing_yearly` is the yearly `smb_ice` delta in
native Chion mass units, and `mean_T_srf_K` is the mean surface temperature.
"""
function yearly_step!(integrator::SimulationIntegrator)
    integrator.sim.model isa BESSIModel ||
        error("yearly_step! currently supports BESSIModel simulations.")
    integrator.time_index == 1 ||
        error("yearly_step! must be called at the start of a forcing year.")

    nsteps = length(integrator.sim.forcing.time_values)
    nsteps > 0 || error("Cannot advance a yearly step without forcing time values.")
    weights_days = integrator.sim.forcing.dt_days
    total_days = sum(weights_days)
    total_days > 0.0 || error("Cannot compute yearly means with non-positive total forcing duration.")

    model = integrator.sim.model
    state = integrator.sim.now
    runtime = integrator.model_runtime.data
    grid = model.grid

    smb_before = _model_smb_ice_vector(model, state, runtime)
    Tsrf_sum = similar(runtime.state.Tsrf)
    fill!(Tsrf_sum, zero(eltype(Tsrf_sum)))

    for time_index in 1:nsteps
        step!(integrator)
        @. Tsrf_sum += runtime.state.Tsrf * weights_days[time_index]
    end

    ice_sheet_net_forcing_yearly = _model_smb_ice_vector(model, state, runtime) .- smb_before
    mean_T_srf_K = _host_vector(Tsrf_sum; copy_array=true) ./ total_days
    _mask_inactive_yearly_outputs!(ice_sheet_net_forcing_yearly, mean_T_srf_K, integrator.model_runtime.active)

    return (
        year=integrator.completed_years,
        mean_T_srf_K=mean_T_srf_K,
        ice_sheet_net_forcing_yearly=ice_sheet_net_forcing_yearly,
        mean_T_srf_K_grid=_yearly_grid(mean_T_srf_K, grid),
        ice_sheet_net_forcing_yearly_grid=_yearly_grid(ice_sheet_net_forcing_yearly, grid),
    )
end

run!(integrator::SimulationIntegrator) = _run_integrator!(integrator)

finalize!(integrator::SimulationIntegrator) = _finalize_integrator!(integrator)

set_active_mask!(integrator::SimulationIntegrator, mask; kwargs...) =
    _set_active_mask!(integrator, mask; kwargs...)

function run!(
    sim::Simulation;
    io::IO=stdout,
)
    integrator = init_integrator(sim; io=io)
    run!(integrator)
    return finalize!(integrator)
end
