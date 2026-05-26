"""
Simulation orchestration for Chion.

`Simulation` owns reference and current model state (`ref` and `now`), while
`SimulationIntegrator` owns initialized runtime state for explicit stepping.
"""

"""
    Simulation(model; forcing, state=nothing, options=SimulationOptions(), output=OutputOptions(), ...)

Couple a model configuration with forcing, reference state, current state, and
execution/output options.
"""
mutable struct Simulation{M <: AbstractSnowModel, D, R <: AbstractState, S <: AbstractState}
    model::M
    domain::D
    forcing::SnowpackForcing
    ref::R
    now::S
    options::SimulationOptions
    output::OutputOptions
end

function Simulation(
    model::AbstractSnowModel;
    forcing::SnowpackForcing,
    state::Union{Nothing, AbstractState}=nothing,
    options::SimulationOptions=SimulationOptions(),
    output::OutputOptions=OutputOptions(),
    years::Union{Nothing, Integer}=nothing,
    backend=nothing,
    history_year_stride::Union{Nothing, Integer}=nothing,
    save=nothing,
    output_dir::Union{Nothing, AbstractString}=nothing,
    netcdf_path::Union{Nothing, AbstractString}=nothing,
    write_outputs::Union{Nothing, Bool}=nothing,
    name::Union{Nothing, AbstractString}=nothing,
    input_label::Union{Nothing, AbstractString}=nothing,
)
    resolved_options = options
    if !isnothing(years) || !isnothing(backend) || !isnothing(history_year_stride) || !isnothing(name) || !isnothing(input_label)
        resolved_options = SimulationOptions(
            name=isnothing(name) ? options.name : name,
            input_label=isnothing(input_label) ? options.input_label : input_label,
            years=isnothing(years) ? options.years : years,
            backend=isnothing(backend) ? options.backend : backend,
            history_year_stride=isnothing(history_year_stride) ? options.history_year_stride : history_year_stride,
        )
    end

    resolved_output = output
    if !isnothing(save) || !isnothing(output_dir) || !isnothing(netcdf_path) || !isnothing(write_outputs)
        resolved_output = OutputOptions(
            save=isnothing(save) ? output.variables : save,
            output_dir=isnothing(output_dir) ? output.output_dir : output_dir,
            netcdf_path=isnothing(netcdf_path) ? output.netcdf_path : netcdf_path,
            write_outputs=isnothing(write_outputs) ? output.write_outputs : write_outputs,
        )
    end

    now_state = isnothing(state) ? initial_state(model) : state
    return Simulation(model, model_domain(model), forcing, reference_state(now_state), now_state, resolved_options, resolved_output)
end

state(sim::Simulation) = sim.now
get_state(sim::Simulation, idx::Int=1) = get_state(sim.now, idx)
print_state(sim::Simulation, idx::Int=1) = print_state(sim.now, idx)

function _normalize_model_name(model)
    name = lowercase(strip(String(model)))
    name in ("bessi", "bessimodel") && return :bessi
    name in ("pdd", "pddmodel") && return :pdd
    name in ("itm", "itmmodel") && return :itm
    error("Unsupported model '$model'. Use `:bessi`, `:pdd`, or `:itm`.")
end

function build_model(model, grid::AbstractSnowpackGrid; kwargs...)
    name = _normalize_model_name(model)
    name == :bessi && return BESSIModel(grid; kwargs...)
    name == :pdd && return PDDModel(grid; kwargs...)
    name == :itm && return ITMModel(grid; kwargs...)
    error("Unsupported model '$model'.")
end

function Simulation(
    model,
    grid::AbstractSnowpackGrid;
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
    println(io, "  domain: ", typeof(sim.domain))
    println(io, "  columns: ", ncols(sim.model.grid))
    println(io, "  forcing steps: ", length(sim.forcing.time_values))
    println(io, "  backend: ", sim.options.backend)
    println(io, "  years: ", sim.options.years)
    println(io, "  netcdf vars: ", isempty(sim.output.variables) ? "(none)" : join(string.(sim.output.variables), ", "))
end

"""
    init_problem!(sim, options)

Validate and prepare the model/forcing context before an initialized run.
"""
function init_problem!(sim::Simulation, options::RunOptions)
    validate_integrator_setup!(sim, options)
    return nothing
end

function init_problem!(
    sim::Simulation;
    options::SimulationOptions=sim.options,
    output::OutputOptions=sim.output,
)
    return init_problem!(sim, _run_options(options, output))
end

"""
    init_integrator(sim; options=sim.options, output=sim.output, io=stdout)

Initialize a runtime that can be advanced with `step!`, `run!`, and
`finalize!`.
"""
function init_integrator(
    sim::Simulation;
    options::SimulationOptions=sim.options,
    output::OutputOptions=sim.output,
    io::IO=stdout,
)
    run_options = _run_options(options, output)
    timings = StepTimingStats()
    init_problem!(sim, run_options)
    model_runtime = init_model_runtime!(sim, run_options, timings)
    diagnostics = init_diagnostics!(sim, model_runtime, run_options, timings)
    output_runtime = init_io!(sim, model_runtime, diagnostics, run_options, timings)
    return _new_integrator(sim, run_options, io, timings, model_runtime, diagnostics, output_runtime)
end

finished(integrator::SimulationIntegrator) = _finished(integrator)

step!(integrator::SimulationIntegrator) = _step_scheduled!(integrator)

step!(integrator::SimulationIntegrator, n::Integer) = _step_n!(integrator, n)

step!(integrator::SimulationIntegrator, Δt_days::Real, force_dt::Bool=true) =
    _step_external!(integrator, Δt_days, force_dt)

_surface_temperature_vector(::BESSIModel, ::CurrentState, runtime) =
    _host_vector(runtime.state.Tsrf; copy_array=true)

function _yearly_grid(values::Vector{Float64}, grid::AbstractSnowpackGrid)
    has_spatial_coords(grid) || return Matrix{Float64}(undef, 0, 0)
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

function _bessi_output_from_options(sim::Simulation, options::RunOptions)
    options.write_netcdf || return nothing
    monthly_mode = :monthly in options.netcdf_variables
    monthly_mode && length(options.netcdf_variables) > 1 &&
        error("`monthly` NetCDF output cannot currently be combined with other selectors.")
    vars = monthly_mode ? MONTHLY_OUTPUT_VARS : state_output_vars(options.netcdf_variables)
    isempty(vars) && return nothing
    return init_state_netcdf(
        resolve_netcdf_path(options),
        options,
        sim.forcing.time_values,
        sim.domain,
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
            _copy_current_state!(sim.now, cpu_state(backend_state))
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
    options::SimulationOptions=sim.options,
    output::OutputOptions=sim.output,
    io::IO=stdout,
    checkpoint_path::AbstractString="",
    checkpoint_year_stride::Integer=1,
)
    run_options = _run_options(options, output)
    checkpoint_path == "" || error("Checkpointing was removed with the simplified BESSI state/output runtime. Re-run without `--checkpoint-path`.")
    checkpoint_year_stride >= 1 || error("`checkpoint_year_stride` must be positive.")
    init_problem!(sim, run_options)
    timings = StepTimingStats()
    clocks = IntegratorClocks(time_ns(), time_ns())
    backend_state = sim.now
    backend_forcing = sim.forcing
    is_gpu = run_options.backend == :gpu
    if is_gpu
        backend_state = time_block!(timings, :gpu_transfer) do
            gpu_state(sim.now)
        end
        backend_forcing = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), sim.forcing)
        end
    end
    workspace = time_block!(timings, :create_workspaces) do
        ColumnarStepWorkspace(backend_state)
    end
    active_indices = is_gpu ? adapt(gpu_storage_type(), collect(1:sim.domain.ncol)) : 1:sim.domain.ncol
    nc = _bessi_output_from_options(sim, run_options)
    monthly_mode = nc !== nothing && (:monthly in run_options.netcdf_variables)
    monthly_state = monthly_mode ? MonthlyState(backend_state) : nothing
    monthly_year_state = monthly_mode ? MonthlyYearState(backend_state; nmonth=run_options.years * 12) : nothing
    record_index = 0
    kwargs = _bessi_step_kwargs(sim.model)
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
    progress = Progress(run_options.years; desc="Running years: ", output=io, showspeed=true)
    for year in 1:run_options.years
        for k in eachindex(sim.forcing.time_values)
            time_counted_block!(timings, :model_step_wall, sim.domain.ncol) do
                step!(backend_state, backend_forcing, k, workspace, active_indices; kwargs...)
            end
            if monthly_mode
                accumulate_monthly!(monthly_state, backend_state)
                if _is_month_boundary(sim.forcing.time_values, Int(k))
                    time_block!(timings, :write_netcdf) do
                        finalize_monthly!(monthly_state, backend_state)
                        store_monthly!(monthly_year_state, monthly_state)
                    end
                    reset_monthly!(monthly_state)
                end
            elseif nc !== nothing
                _sync_bessi_state!(sim, backend_state, is_gpu, timings)
                record_index += 1
                time_block!(timings, :write_netcdf) do
                    write_nc!(nc, sim.now, record_index, sim.domain)
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
                println(io, year_log_line(record))
                flush(io)
            end
        end
        copyto!(prev_year_summary.thickness, year_summary.thickness)
        copyto!(prev_year_summary.wet_mass, year_summary.wet_mass)
        copyto!(prev_year_summary.base_mass, year_summary.base_mass)
        next!(progress)
    end
    if monthly_mode
        time_block!(timings, :write_netcdf) do
            write_monthly_year_nc!(nc, monthly_year_state, 1, sim.domain)
        end
        record_index = monthly_year_state.count
    end
    _sync_bessi_state!(sim, backend_state, is_gpu, timings)
    status = :complete
    nc_path = nc === nothing ? "" : resolve_netcdf_path(run_options)
    nc !== nothing && close_output!(nc, status, record_index)
    run_wall_sec = (time_ns() - clocks.run_wall_t0) * 1.0e-9
    simulation_wall_sec = (time_ns() - clocks.simulation_wall_t0) * 1.0e-9
    result = SimulationResult(
        history,
        status,
        run_options.years,
        timings,
        simulation_wall_sec,
        run_wall_sec,
        nc_path,
        "",
        "",
    )
    print_run_report(
        io,
        run_options,
        sim.forcing.time_values,
        history,
        status,
        simulation_wall_sec,
        run_wall_sec,
        timings;
        nc_path=nc_path,
        summary_path="",
        history_csv_path="",
    )
    return result
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
    runtime = integrator.model_runtime.backend
    grid = integrator.model_runtime.grid

    smb_before = _model_smb_ice_vector(model, state, runtime)
    Tsrf_sum = zeros(Float64, integrator.model_runtime.ncol)

    for time_index in 1:nsteps
        step!(integrator)
        Tsrf = _surface_temperature_vector(model, state, runtime)
        @. Tsrf_sum += Tsrf * weights_days[time_index]
    end

    ice_sheet_net_forcing_yearly = _model_smb_ice_vector(model, state, runtime) .- smb_before
    mean_T_srf_K = Tsrf_sum ./ total_days
    _mask_inactive_yearly_outputs!(ice_sheet_net_forcing_yearly, mean_T_srf_K, integrator.model_runtime.active)

    return (
        year=integrator.completed_years,
        mean_T_srf_K=mean_T_srf_K,
        ice_sheet_net_forcing_yearly=ice_sheet_net_forcing_yearly,
        mean_T_srf_K_grid=_yearly_grid(mean_T_srf_K, grid),
        ice_sheet_net_forcing_yearly_grid=_yearly_grid(ice_sheet_net_forcing_yearly, grid),
    )
end

run!(
    integrator::SimulationIntegrator;
    checkpoint_path::AbstractString="",
    checkpoint_year_stride::Integer=1,
) = _run_integrator!(
    integrator;
    checkpoint_path=checkpoint_path,
    checkpoint_year_stride=checkpoint_year_stride,
)

finalize!(integrator::SimulationIntegrator) = _finalize_integrator!(integrator)

set_forcing!(integrator::SimulationIntegrator; kwargs...) =
    _set_forcing!(integrator; kwargs...)

set_active_mask!(integrator::SimulationIntegrator, mask; kwargs...) =
    _set_active_mask!(integrator, mask; kwargs...)

function run!(
    sim::Simulation;
    options::SimulationOptions=sim.options,
    output::OutputOptions=sim.output,
    io::IO=stdout,
    checkpoint_path::AbstractString="",
    checkpoint_year_stride::Integer=1,
)
    integrator = init_integrator(sim; options=options, output=output, io=io)
    run!(
        integrator;
        checkpoint_path=checkpoint_path,
        checkpoint_year_stride=checkpoint_year_stride,
    )
    return finalize!(integrator)
end
