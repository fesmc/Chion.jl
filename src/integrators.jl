"""Initialized integrator type and lifecycle helpers for `Simulation`."""

struct IntegratorClocks
    run_wall_t0::Int
    simulation_wall_t0::Int
end

struct IntegratorModelRuntime
    runtime
    ncol::Int
    grid
end

mutable struct DiagnosticsRuntime
    selected::Set{Symbol}
    need_step_outputs::Bool
    need_monthly_outputs::Bool
    need_layer_outputs::Bool
    need_last_year_smb_delta::Bool
    need_step_diagnostics::Bool
    prev
    final
    backend_year_summary
    step_summary
    backend_step_summary
    deltas
    previous
    previous_year_smb_ice::Vector{Float64}
    history::Vector{NamedTuple}
end

mutable struct OutputRuntime
    schedule
    writer
    nc_path::String
    monthly_sums
    monthly_count::Vector{Int32}
    step_vectors
    steps_written::Int
end

mutable struct StepperRuntime
    time_index::Int
    completed_years::Int
    current_forcing::SnowpackForcing
    progress
end

mutable struct SimulationIntegrator
    sim
    options::RunOptions
    io::IO
    timings::StepTimingStats
    clocks::IntegratorClocks
    model_runtime::IntegratorModelRuntime
    diagnostics::DiagnosticsRuntime
    output::OutputRuntime
    stepper::StepperRuntime
    finalized::Bool
    result::Union{Nothing, SimulationResult}
end

function _single_step_forcing_template(forcing::SnowpackForcing)
    return SnowpackForcing(
        dt_days=[first(forcing.dt_days)],
        ncol=size(forcing.air_temperature, 1),
        air_temperature=forcing.air_temperature[:, 1:1],
        snowfall_rate=forcing.snowfall_rate[:, 1:1],
        rainfall_rate=forcing.rainfall_rate[:, 1:1],
        shortwave_down=forcing.shortwave_down[:, 1:1],
        latitude_deg=forcing.latitude_deg[:, 1:1],
        wind_speed=forcing.wind_speed[:, 1:1],
        q_lw_down=forcing.q_lw_down[:, 1:1],
        has_q_lw_down=forcing.has_q_lw_down[:, 1:1],
        q_sh=forcing.q_sh[:, 1:1],
        has_q_sh=forcing.has_q_sh[:, 1:1],
        q_lh=forcing.q_lh[:, 1:1],
        has_q_lh=forcing.has_q_lh[:, 1:1],
        time_values=[first(forcing.time_values)],
    )
end

function init_stepper_state!(context, options::RunOptions, io::IO)
    progress = Progress(options.years; desc="Running years: ", output=io, showspeed=true)
    return StepperRuntime(1, 0, _single_step_forcing_template(context.forcing), progress)
end

function _new_integrator(sim, options::RunOptions, io::IO, timings::StepTimingStats, model_runtime, diagnostics, output, stepper)
    return SimulationIntegrator(
        sim,
        options,
        io,
        timings,
        IntegratorClocks(time_ns(), time_ns()),
        model_runtime,
        diagnostics,
        output,
        stepper,
        false,
        nothing,
    )
end

_finished(integrator::SimulationIntegrator) =
    integrator.stepper.completed_years >= integrator.options.years

function _external_forcing_for_runtime(integrator::SimulationIntegrator)
    forcing = integrator.stepper.current_forcing
    integrator.options.backend == :gpu || return forcing
    return time_block!(integrator.timings, :gpu_transfer) do
        adapt(gpu_storage_type(), forcing)
    end
end

_scheduled_forcing_for_runtime(integrator::SimulationIntegrator) =
    integrator.model_runtime.runtime.step_fields

function _advance_with_forcing!(integrator::SimulationIntegrator, forcing::SnowpackForcing, time_index::Int)
    integrator.finalized && error("Cannot step a finalized integrator.")
    _finished(integrator) && error("Cannot step an integrator that has already completed all years.")

    model = integrator.sim.model
    state = integrator.sim.now
    model_runtime = integrator.model_runtime
    runtime = model_runtime.runtime

    time_counted_block!(integrator.timings, :model_step_wall, model_runtime.ncol) do
        step_model!(model, state, runtime, forcing, time_index)
    end

    accumulate_step_diagnostics!(integrator)
    maybe_write_step_outputs!(integrator)

    if integrator.stepper.time_index == length(integrator.sim.forcing.time_values)
        _complete_year!(integrator)
        integrator.stepper.time_index = 1
    else
        integrator.stepper.time_index += 1
    end
    return nothing
end

function _step_scheduled!(integrator::SimulationIntegrator)
    forcing = _scheduled_forcing_for_runtime(integrator)
    return _advance_with_forcing!(integrator, forcing, integrator.stepper.time_index)
end

function _step_n!(integrator::SimulationIntegrator, n::Integer)
    n >= 0 || error("Step count must be non-negative.")
    for _ in 1:Int(n)
        _step_scheduled!(integrator)
    end
    return nothing
end

function _step_external!(integrator::SimulationIntegrator, Δt_days::Real, force_dt::Bool=true)
    force_dt || error("Chion's initialized stepper requires `force_dt=true`, matching the FastIsostasy coupling pattern.")
    Δt_days > 0 || error("`Δt_days` must be positive.")
    integrator.stepper.current_forcing.dt_days[1] = Float64(Δt_days)
    forcing = _external_forcing_for_runtime(integrator)
    return _advance_with_forcing!(integrator, forcing, 1)
end

function _run_integrator!(integrator::SimulationIntegrator)
    while !_finished(integrator)
        _step_scheduled!(integrator)
    end
    return nothing
end

function _assign_numeric_step_field!(dest::AbstractMatrix, value, ncol::Int, name::AbstractString; transform=identity)
    isnothing(value) && return false
    dest[:, :] .= transform.(_forcing_numeric_matrix(value, ncol, 1, name))
    return true
end

function _assign_bool_step_field!(dest::AbstractMatrix{Bool}, value, ncol::Int, name::AbstractString)
    isnothing(value) && return false
    dest[:, :] .= _forcing_bool_matrix(value, ncol, 1, name)
    return true
end

function _set_forcing!(
    integrator::SimulationIntegrator;
    air_temperature=nothing,
    snowfall_rate=nothing,
    rainfall_rate=nothing,
    air_temperature_c=nothing,
    snowfall_mm_day=nothing,
    rainfall_mm_day=nothing,
    shortwave_down=nothing,
    wind_speed=nothing,
    q_lw_down=nothing,
    has_q_lw_down=nothing,
    q_sh=nothing,
    has_q_sh=nothing,
    q_lh=nothing,
    has_q_lh=nothing,
    latitude_deg=nothing,
    time_value=nothing,
)
    f = integrator.stepper.current_forcing
    ncol = integrator.model_runtime.ncol
    has_native = !isnothing(air_temperature) || !isnothing(snowfall_rate) || !isnothing(rainfall_rate)
    has_user = !isnothing(air_temperature_c) || !isnothing(snowfall_mm_day) || !isnothing(rainfall_mm_day)
    has_native && has_user && error("Pass either model-native forcing fields or user-facing fields, not both.")
    if has_user
        _assign_numeric_step_field!(f.air_temperature, air_temperature_c, ncol, "air_temperature_c"; transform=x -> x + 273.15)
        _assign_numeric_step_field!(f.snowfall_rate, snowfall_mm_day, ncol, "snowfall_mm_day"; transform=x -> x / 86_400.0)
        _assign_numeric_step_field!(f.rainfall_rate, rainfall_mm_day, ncol, "rainfall_mm_day"; transform=x -> x / 86_400.0)
    else
        _assign_numeric_step_field!(f.air_temperature, air_temperature, ncol, "air_temperature")
        _assign_numeric_step_field!(f.snowfall_rate, snowfall_rate, ncol, "snowfall_rate")
        _assign_numeric_step_field!(f.rainfall_rate, rainfall_rate, ncol, "rainfall_rate")
    end
    _assign_numeric_step_field!(f.shortwave_down, shortwave_down, ncol, "shortwave_down")
    !isnothing(latitude_deg) && (f.latitude_deg[:, :] .= _forcing_column_metadata_matrix(latitude_deg, ncol, 1, "latitude_deg"))
    _assign_numeric_step_field!(f.wind_speed, wind_speed, ncol, "wind_speed")
    if _assign_numeric_step_field!(f.q_lw_down, q_lw_down, ncol, "q_lw_down") && isnothing(has_q_lw_down)
        fill!(f.has_q_lw_down, true)
    end
    _assign_bool_step_field!(f.has_q_lw_down, has_q_lw_down, ncol, "has_q_lw_down")
    if _assign_numeric_step_field!(f.q_sh, q_sh, ncol, "q_sh") && isnothing(has_q_sh)
        fill!(f.has_q_sh, true)
    end
    _assign_bool_step_field!(f.has_q_sh, has_q_sh, ncol, "has_q_sh")
    if _assign_numeric_step_field!(f.q_lh, q_lh, ncol, "q_lh") && isnothing(has_q_lh)
        fill!(f.has_q_lh, true)
    end
    _assign_bool_step_field!(f.has_q_lh, has_q_lh, ncol, "has_q_lh")
    if !isnothing(time_value)
        f.time_values[1] = DateTime(time_value)
        f.day_of_year[1] = _calendar_day_of_year(f.time_values[1])
        f.solar_longitude_deg[1] = _solar_longitude_deg_from_calendar_day(f.day_of_year[1])
    end
    return integrator
end

function _finalize_integrator!(integrator::SimulationIntegrator)
    integrator.finalized && return integrator.result

    history = integrator.diagnostics.history
    simulation_wall_sec = (time_ns() - integrator.clocks.simulation_wall_t0) * 1.0e-9
    status = _finished(integrator) ? :complete : :incomplete
    years_completed = completed_year_count(history, status, integrator.options.years)
    summary_path = ""
    history_csv_path = ""

    if integrator.options.write_outputs
        mkpath(integrator.options.output_dir)
        summary_path = joinpath(integrator.options.output_dir, "$(integrator.options.name)_summary.txt")
        history_csv_path = joinpath(integrator.options.output_dir, "$(integrator.options.name)_history.csv")
        time_block!(integrator.timings, :write_summary_text) do
            write_run_summary(
                summary_path,
                integrator.options,
                integrator.sim.forcing.time_values,
                integrator.model_runtime.ncol,
                history,
                status,
                integrator.timings,
            )
        end
        time_block!(integrator.timings, :write_history_csv) do
            write_run_history_csv(history_csv_path, history)
        end
    end

    finalize_output_runtime!(integrator, status, years_completed)
    finalize_state!(integrator.sim.model, integrator.sim.now, integrator.model_runtime.runtime, integrator.options, integrator.timings)

    run_wall_sec = (time_ns() - integrator.clocks.run_wall_t0) * 1.0e-9
    print_run_report(
        integrator.io,
        integrator.options,
        integrator.sim.forcing.time_values,
        history,
        status,
        simulation_wall_sec,
        run_wall_sec,
        integrator.timings;
        nc_path=integrator.output.nc_path,
        summary_path=summary_path,
        history_csv_path=history_csv_path,
    )
    result = SimulationResult(
        history,
        status,
        years_completed,
        integrator.timings,
        simulation_wall_sec,
        run_wall_sec,
        integrator.output.nc_path,
        summary_path,
        history_csv_path,
    )
    integrator.finalized = true
    integrator.result = result
    return result
end
