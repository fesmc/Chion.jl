"""
FastIsostasy-style simulation orchestration for Chion.

`Simulation` owns reference and current model state (`ref` and `now`), while
`SimulationIntegrator` owns initialized runtime state for explicit stepping.
"""

"""
    Simulation(model; forcing, state=nothing, options=SimulationOptions(), output=OutputOptions(), ...)

Couple a model configuration with forcing, reference state, current state, and
execution/output options.
"""
mutable struct Simulation{M <: AbstractSnowModel, S <: AbstractSnowModelState}
    model::M
    forcing::SnowpackForcing
    ref::S
    now::S
    options::SimulationOptions
    output::OutputOptions
end

function Simulation(
    model::AbstractSnowModel;
    forcing::SnowpackForcing,
    state::Union{Nothing, AbstractSnowModelState}=nothing,
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
    return Simulation(model, forcing, deepcopy(now_state), now_state, resolved_options, resolved_output)
end

state(sim::Simulation) = sim.now

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

_surface_temperature_vector(::BESSIModel, ::BESSIState, runtime) =
    _host_vector(runtime.domain.Tsrf; copy_array=true)

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
