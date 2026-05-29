"""Initialized integrator type and lifecycle helpers for `Simulation`."""

function _single_step_forcing_view(forcing::SnowpackForcing, dt_days::Real)
    return SnowpackForcing(
        dt_days=[Float64(dt_days)],
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
        relative_humidity=forcing.relative_humidity[:, 1:1],
        has_relative_humidity=forcing.has_relative_humidity[:, 1:1],
        surface_height=forcing.surface_height[:, 1:1],
        air_pressure=forcing.air_pressure[:, 1:1],
        prescribed_albedo=forcing.prescribed_albedo[:, 1:1],
        has_prescribed_albedo=forcing.has_prescribed_albedo[:, 1:1],
        time_values=[first(forcing.time_values)],
    )
end

function _new_integrator(
    sim,
    options::RunOptions,
    io::IO,
    timings::StepTimingStats,
    model_runtime,
    diagnostics,
    output;
    time_index::Integer=1,
    completed_years::Integer=0,
)
    return SimulationIntegrator(
        sim,
        options,
        io,
        timings,
        IntegratorClocks(time_ns(), time_ns()),
        model_runtime,
        diagnostics,
        output,
        Int(time_index),
        Int(completed_years),
        nothing,
        false,
        nothing,
    )
end

_finished(integrator::SimulationIntegrator) =
    integrator.completed_years >= integrator.options.years

_scheduled_forcing_for_runtime(integrator::SimulationIntegrator) =
    integrator.model_runtime.backend.step_fields

function _copy_forcing!(dest::SnowpackForcing, src::SnowpackForcing)
    dest.dt_days .= src.dt_days
    dest.day_of_year .= src.day_of_year
    dest.solar_longitude_deg .= src.solar_longitude_deg
    dest.air_temperature .= src.air_temperature
    dest.snowfall_rate .= src.snowfall_rate
    dest.rainfall_rate .= src.rainfall_rate
    dest.shortwave_down .= src.shortwave_down
    dest.latitude_deg .= src.latitude_deg
    dest.wind_speed .= src.wind_speed
    dest.q_lw_down .= src.q_lw_down
    dest.has_q_lw_down .= src.has_q_lw_down
    dest.q_sh .= src.q_sh
    dest.has_q_sh .= src.has_q_sh
    dest.q_lh .= src.q_lh
    dest.has_q_lh .= src.has_q_lh
    dest.relative_humidity .= src.relative_humidity
    dest.has_relative_humidity .= src.has_relative_humidity
    dest.surface_height .= src.surface_height
    dest.air_pressure .= src.air_pressure
    dest.prescribed_albedo .= src.prescribed_albedo
    dest.has_prescribed_albedo .= src.has_prescribed_albedo
    return dest
end

function sync_forcing!(integrator::SimulationIntegrator)
    backend = integrator.model_runtime.backend
    getproperty(backend, :is_gpu) || return integrator
    time_block!(integrator.timings, :gpu_transfer) do
        _copy_forcing!(backend.step_fields, integrator.forcing)
    end
    return integrator
end

function _advance_with_forcing!(integrator::SimulationIntegrator, forcing::SnowpackForcing, time_index::Int)
    integrator.finalized && error("Cannot step a finalized integrator.")
    _finished(integrator) && error("Cannot step an integrator that has already completed all years.")

    model = integrator.sim.model
    state = integrator.sim.now
    model_runtime = integrator.model_runtime
    runtime = model_runtime.backend

    time_counted_block!(integrator.timings, :model_step_wall, model_runtime.ncol) do
        step_model!(model, state, model_runtime, forcing, time_index)
    end

    if integrator.time_index == length(integrator.sim.forcing.time_values)
        integrator.completed_years += 1
        integrator.time_index = 1
    else
        integrator.time_index += 1
    end
    return nothing
end

function _step_scheduled!(integrator::SimulationIntegrator)
    forcing = _scheduled_forcing_for_runtime(integrator)
    return _advance_with_forcing!(integrator, forcing, integrator.time_index)
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
    forcing = _single_step_forcing_view(integrator.forcing, Δt_days)
    if integrator.options.backend == :gpu
        forcing = time_block!(integrator.timings, :gpu_transfer) do
            adapt(gpu_storage_type(), forcing)
        end
    end
    return _advance_with_forcing!(integrator, forcing, 1)
end

function _should_checkpoint_after_year(integrator::SimulationIntegrator, checkpoint_path::AbstractString, checkpoint_year_stride::Integer)
    isempty(checkpoint_path) && return false
    stride = Int(checkpoint_year_stride)
    stride >= 0 || error("`checkpoint_year_stride` must be >= 0.")
    stride == 0 && return false
    _finished(integrator) && return false
    return mod(integrator.completed_years, stride) == 0
end

function _run_integrator!(
    integrator::SimulationIntegrator;
    checkpoint_path::AbstractString="",
    checkpoint_year_stride::Integer=1,
)
    if integrator.sim.model isa BESSIModel
        return _run_bessi_integrator!(
            integrator;
            checkpoint_path=checkpoint_path,
            checkpoint_year_stride=checkpoint_year_stride,
        )
    end
    while !_finished(integrator)
        _step_scheduled!(integrator)
        checkpoint_path == "" || error("Checkpointing was removed with the simplified state/output runtime.")
    end
    return nothing
end

@inline _active_mask_value(value::Bool) = value
@inline _active_mask_value(value) = isfinite(Float64(value)) && Float64(value) != 0.0

function _active_mask_vector(mask, grid, ncol::Int)
    data = collect(mask)
    if ndims(data) == 1
        length(data) == ncol || error("Active mask vector length must match the model column count.")
        return [_active_mask_value(data[idx]) for idx in eachindex(data)]
    elseif ndims(data) == 2
        any(isnothing, (grid.x, grid.y, grid.js, grid.is, grid.mask)) &&
            error("2-D active masks require a grid with spatial coordinates.")
        size(data) == size(grid.mask) || error("Active mask matrix must have size $(size(grid.mask)), got $(size(data)).")
        return [_active_mask_value(data[grid.js[col], grid.is[col]]) for col in 1:ncol]
    end
    error("Active mask must be a vector of model columns or a y-x matrix matching `grid.mask`.")
end

function _set_active_mask!(
    integrator::SimulationIntegrator,
    mask;
    reset_newly_inactive::Bool=true,
)
    integrator.finalized && error("Cannot update the active mask of a finalized integrator.")
    active = _active_mask_vector(mask, integrator.model_runtime.grid, integrator.model_runtime.ncol)
    old_active = integrator.model_runtime.active
    newly_inactive = reset_newly_inactive ? findall(old_active .& .!active) : Int[]
    if !isempty(newly_inactive)
        _reset_model_columns!(
            integrator.sim.model,
            integrator.sim.now,
            integrator.model_runtime.backend,
            newly_inactive,
        )
    end
    _set_model_runtime_active_indices!(integrator.model_runtime, active)
    return integrator
end

function _finalize_integrator!(integrator::SimulationIntegrator)
    integrator.finalized && return integrator.result

    simulation_wall_sec = (time_ns() - integrator.clocks.simulation_wall_t0) * 1.0e-9
    status = _finished(integrator) ? :complete : :incomplete
    years_completed = integrator.completed_years
    finalize_state!(integrator.sim.model, integrator.sim.now, integrator.model_runtime.backend, integrator.options, integrator.timings)

    run_wall_sec = (time_ns() - integrator.clocks.run_wall_t0) * 1.0e-9
    result = SimulationResult(
        isnothing(integrator.diagnostics) ? NamedTuple[] : integrator.diagnostics,
        status,
        years_completed,
        integrator.timings,
        simulation_wall_sec,
        run_wall_sec,
        isnothing(integrator.output) ? "" : integrator.output.nc_path,
        isnothing(integrator.output) ? "" : integrator.output.summary_path,
        isnothing(integrator.output) ? "" : integrator.output.history_csv_path,
    )
    integrator.finalized = true
    integrator.result = result
    return result
end
