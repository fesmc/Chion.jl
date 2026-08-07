"""Initialized integrator type and lifecycle helpers for `Simulation`."""

function _single_step_forcing_view(forcing::SnowpackForcing, dt_days::Real)
    return SnowpackForcing(
        ;
        dt_days=[Float64(dt_days)],
        ncol=size(forcing.air_temperature, 1),
        time_values=[first(forcing.time_values)],
        (name => getproperty(forcing, name)[:, 1:1] for name in _FORCING_MATRIX_FIELD_NAMES)...,
    )
end

@inline _copy_forcing_field!(dest, src) = copyto!(dest, src)
@inline _copy_forcing_field!(::ConstantForcingMatrix, ::ConstantForcingMatrix) = nothing
@inline _copy_forcing_field!(dest::ColumnForcingMatrix, src::ColumnForcingMatrix) =
    copyto!(dest.values, src.values)
@inline _copy_forcing_field!(dest::TimeForcingMatrix, src::TimeForcingMatrix) =
    copyto!(dest.values, src.values)

function _new_integrator(
    sim,
    io::IO,
    timings::StepTimingStats,
    model_runtime;
    time_index::Integer=1,
    completed_years::Integer=0,
)
    return SimulationIntegrator(
        sim,
        io,
        timings,
        Int(time_ns()),
        model_runtime,
        YearMetrics[],
        "",
        Int(time_index),
        Int(completed_years),
        nothing,
    )
end

_finished(integrator::SimulationIntegrator) =
    integrator.completed_years >= integrator.sim.options.years

_scheduled_forcing_for_runtime(integrator::SimulationIntegrator) =
    integrator.model_runtime.data.step_fields

function _copy_forcing!(dest::SnowpackForcing, src::SnowpackForcing)
    for name in _FORCING_COPY_FIELD_NAMES
        _copy_forcing_field!(getproperty(dest, name), getproperty(src, name))
    end
    return dest
end
function _copy_forcing!(dest::PDDForcing, src::SnowpackForcing)
    copyto!(dest.dt_days, src.dt_days)
    _copy_forcing_field!(dest.air_temperature, src.air_temperature)
    _copy_forcing_field!(dest.snowfall_rate, src.snowfall_rate)
    _copy_forcing_field!(dest.rainfall_rate, src.rainfall_rate)
    return dest
end

function _copy_forcing!(dest::ITMForcing, src::SnowpackForcing)
    copyto!(dest.dt_days, src.dt_days)
    for name in (:air_temperature, :snowfall_rate, :rainfall_rate, :shortwave_down,
                 :q_sw_net, :has_q_sw_net, :latitude_deg, :surface_height,
                 :ice_thickness, :annual_pdd)
        _copy_forcing_field!(getproperty(dest, name), getproperty(src, name))
    end
    return dest
end

function sync_forcing!(integrator::SimulationIntegrator)
    data = integrator.model_runtime.data
    getproperty(data, :is_gpu) || return integrator
    time_block!(integrator.timings, :gpu_transfer) do
        _copy_forcing!(data.step_fields, integrator.sim.forcing)
    end
    return integrator
end

function _advance_with_forcing!(integrator::SimulationIntegrator, forcing, time_index::Int)
    return _advance_with_forcing_range!(
        integrator, forcing, time_index, time_index, integrator.time_index,
    )
end

function _advance_with_forcing!(integrator::SimulationIntegrator, forcing, time_range)
    first_time = Int(first(time_range))
    last_time = Int(last(time_range))
    first_time <= last_time || return nothing
    first_time == integrator.time_index || error("Scheduled forcing range must start at the integrator time index.")
    return _advance_with_forcing_range!(integrator, forcing, first_time, last_time, last_time)
end

function _advance_with_forcing_range!(
    integrator::SimulationIntegrator,
    forcing,
    first_time::Int,
    last_time::Int,
    schedule_stop::Int,
)
    integrator.result !== nothing && error("Cannot step a finalized integrator.")
    _finished(integrator) && error("Cannot step an integrator that has already completed all years.")
    nsteps = length(integrator.sim.forcing.time_values)
    last_time <= nsteps || error("Scheduled forcing range cannot cross a forcing year boundary.")

    model = integrator.sim.model
    state = integrator.sim.now
    model_runtime = integrator.model_runtime
    step_count = last_time - first_time + 1

    time_counted_block!(integrator.timings, :model_step_wall, step_count) do
        step_input = first_time == last_time ? first_time : first_time:last_time
        step_model!(model, state, model_runtime, forcing, step_input)
    end

    if schedule_stop == nsteps
        integrator.completed_years += 1
        integrator.time_index = 1
    else
        integrator.time_index = schedule_stop + 1
    end
    return nothing
end

function _step_scheduled!(integrator::SimulationIntegrator)
    forcing = _scheduled_forcing_for_runtime(integrator)
    return _advance_with_forcing!(integrator, forcing, integrator.time_index)
end

function _step_scheduled_range!(integrator::SimulationIntegrator, time_range)
    forcing = _scheduled_forcing_for_runtime(integrator)
    return _advance_with_forcing!(integrator, forcing, time_range)
end

function _step_n!(integrator::SimulationIntegrator, n::Integer)
    n >= 0 || error("Step count must be non-negative.")
    for _ in 1:Int(n)
        _step_scheduled!(integrator)
    end
    return nothing
end

@inline _external_gpu_forcing(::BESSIModel, forcing::SnowpackForcing) =
    adapt(gpu_storage_type(), get_fields(forcing))
@inline _external_gpu_forcing(::PDDModel, forcing::SnowpackForcing) =
    _gpu_forcing_view(PDDForcing(forcing))
@inline _external_gpu_forcing(::ITMModel, forcing::SnowpackForcing) =
    _gpu_forcing_view(ITMForcing(forcing))

function _step_external!(integrator::SimulationIntegrator, Δt_days::Real, force_dt::Bool=true)
    force_dt || error("Chion's initialized stepper requires `force_dt=true`, matching the FastIsostasy coupling pattern.")
    Δt_days > 0 || error("`Δt_days` must be positive.")
    forcing = _single_step_forcing_view(integrator.sim.forcing, Δt_days)
    if integrator.sim.options.backend == :gpu
        forcing = time_block!(integrator.timings, :gpu_transfer) do
            _external_gpu_forcing(integrator.sim.model, forcing)
        end
    end
    return _advance_with_forcing!(integrator, forcing, 1)
end

function _run_integrator!(integrator::SimulationIntegrator)
    while !_finished(integrator)
        _step_scheduled!(integrator)
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
    integrator.result !== nothing && error("Cannot update the active mask of a finalized integrator.")
    grid = integrator.sim.model.grid
    active = _active_mask_vector(mask, grid, ncols(grid))
    old_active = integrator.model_runtime.active
    newly_inactive = reset_newly_inactive ? findall(old_active .& .!active) : Int[]
    if !isempty(newly_inactive)
        _reset_model_columns!(
            integrator.sim.model,
            integrator.sim.now,
            integrator.model_runtime.data,
            newly_inactive,
        )
    end
    _set_model_runtime_active_indices!(integrator.model_runtime, active)
    return integrator
end

function _finalize_integrator!(integrator::SimulationIntegrator)
    integrator.result !== nothing && return integrator.result

    status = _finished(integrator) ? :complete : :incomplete
    years_completed = integrator.completed_years
    finalize_state!(integrator.sim.model, integrator.sim.now, integrator.model_runtime.data, integrator.sim.options, integrator.timings)

    run_wall_sec = (time_ns() - integrator.wall_t0) * 1.0e-9
    result = SimulationResult(
        integrator.history,
        status,
        years_completed,
        integrator.timings,
        run_wall_sec,
        integrator.netcdf_path,
    )
    integrator.result = result
    return result
end
