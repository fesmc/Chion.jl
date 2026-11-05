"""Model runtime preparation, backend transfer, and state finalization."""

function _validate_model_outputs!(model, options::RunOptions)
    options.write_netcdf || return nothing
    allowed = supports_monthly_output(model) ?
        (output_variables(model)..., monthly_output_variables(model)..., :all, :none, :monthly) :
        (output_variables(model)..., :all, :none)
    unsupported = setdiff(options.netcdf_variables, allowed)
    isempty(unsupported) || error(
        "Unsupported NetCDF variables for $(nameof(typeof(model))): $(join(string.(unsupported), ", ")).",
    )
    return nothing
end

function _validate_model_forcing!(model::BESSIModel, forcing::SnowpackForcing)
    if model.diurnal_shortwave_substeps
        all(isfinite, forcing.latitude_deg) || error("`latitude_deg` is required in `SnowpackForcing` for the default BESSI diurnal substeps; " *
            "provide it or set `diurnal_shortwave_substeps=false` (e.g. for sub-daily forcing).")
    end
    _validate_albedo_forcing!(albedo_scheme(model.c), forcing)
    return nothing
end

_validate_albedo_forcing!(::AbstractAlbedoScheme, ::SnowpackForcing) = nothing
function _validate_albedo_forcing!(::PrescribedAlbedo, forcing::SnowpackForcing)
    all(forcing.has_prescribed_albedo) || error("`prescribed_albedo` is required for every column and timestep when BESSI uses `albedo=:prescribed`.")
    return nothing
end

_validate_model_forcing!(::PDDModel, ::SnowpackForcing) = nothing

function _validate_model_forcing!(::ITMModel, forcing::SnowpackForcing)
    all(isfinite, forcing.surface_height) || error("`surface_height` is required for ITMModel.")
    all(isfinite, forcing.ice_thickness) || error("`ice_thickness` is required for ITMModel.")
    all(isfinite, forcing.annual_pdd) || error("`annual_pdd` is required for ITMModel.")
    all(isfinite, forcing.latitude_deg) || error("`latitude_deg` is required for ITMModel.")
    all(>=(0.0), forcing.ice_thickness) || error("`ice_thickness` must be non-negative for ITMModel.")
    all(>=(0.0), forcing.annual_pdd) || error("`annual_pdd` must be non-negative for ITMModel.")
    return nothing
end
function validate_integrator_setup!(sim, options::RunOptions)
    model = sim.model
    forcing = sim.forcing
    grid = model.grid
    ncol = ncols(grid)
    size(forcing.air_temperature, 1) == ncol || error("Forcing column count must match the model column count.")
    _validate_model_forcing!(model, forcing)
    spatial_grid = has_spatial_coords(grid)
    options.write_netcdf && !spatial_grid && error("NetCDF output requires a grid with spatial coordinates.")
    spatial_grid && length(grid.js) != ncol && error("Grid point count must match the domain column count.")
    _validate_model_outputs!(model, options)
    return nothing
end

_model_forcing_view(::BESSIModel, forcing::SnowpackForcing) = get_fields(forcing)
_model_forcing_view(::PDDModel, forcing::SnowpackForcing) = PDDForcing(forcing)
_model_forcing_view(::ITMModel, forcing::SnowpackForcing) = ITMForcing(forcing)

"""Forcing fields stepped by `model`, in the model's number type (no copy for Float64)."""
model_forcing(model::AbstractSnowModel, forcing::SnowpackForcing) =
    forcing_fields_as(number_type(model.c), _model_forcing_view(model, forcing))

"""Copy `state` to `backend`; BESSI layer fields are transposed for coalesced GPU access."""
to_backend(::ThreadsBackend, state) = state
to_backend(backend::GPUBackend, state::BESSIState) = gpu_state(state, storage_type(backend))
to_backend(backend::GPUBackend, state) = adapt(storage_type(backend), state)

# Forcing fields on `backend`. PDD and ITM keep `dt_days` on the host because their
# kernels receive the step length as a scalar.
_backend_forcing(::ThreadsBackend, step_fields) = step_fields
_backend_forcing(backend::GPUBackend, step_fields::NamedTuple) = adapt(storage_type(backend), step_fields)
function _backend_forcing(backend::GPUBackend, forcing::PDDForcing)
    to_device(field) = adapt(storage_type(backend), field)
    return PDDForcing(
        forcing.dt_days,
        to_device(forcing.air_temperature),
        to_device(forcing.snowfall_rate),
        to_device(forcing.rainfall_rate),
    )
end
function _backend_forcing(backend::GPUBackend, forcing::ITMForcing)
    to_device(field) = adapt(storage_type(backend), field)
    return ITMForcing(
        forcing.dt_days,
        to_device(forcing.air_temperature),
        to_device(forcing.snowfall_rate),
        to_device(forcing.rainfall_rate),
        to_device(forcing.shortwave_down),
        to_device(forcing.q_sw_net),
        to_device(forcing.has_q_sw_net),
        to_device(forcing.latitude_deg),
        to_device(forcing.surface_height),
        to_device(forcing.ice_thickness),
        to_device(forcing.annual_pdd),
    )
end

_runtime_workspace(state::BESSIState) = ColumnarStepWorkspace(state)
_runtime_workspace(state) = nothing

function prepare_runtime!(model::AbstractSnowModel, state, forcing::SnowpackForcing, options::RunOptions, timings::StepTimingStats)
    return _prepare_runtime(options.backend, timings, state, model_forcing(model, forcing))
end

function _prepare_runtime(::ThreadsBackend, timings::StepTimingStats, state, step_fields)
    workspace = time_block!(timings, :create_workspaces) do
        _runtime_workspace(state)
    end
    return ModelRuntimeData(state, step_fields, workspace, false)
end

function _prepare_runtime(backend::GPUBackend, timings::StepTimingStats, state, step_fields)
    _require_functional(backend)
    state, step_fields, workspace = time_block!(timings, :gpu_transfer) do
        backend_state = to_backend(backend, state)
        backend_state, _backend_forcing(backend, step_fields), _runtime_workspace(backend_state)
    end
    return ModelRuntimeData(state, step_fields, workspace, is_gpu(backend))
end

_backend_active_indices(indices::Vector{Int}, data) =
    data.is_gpu ? adapt(gpu_storage_type(), indices) : indices

@kernel function _reset_bessi_columns_kernel!(
    fields,
    Ntot::Int,
    inactive_indices,
    density_init,
    temperature_init,
    surface_temperature_init,
    albedo_init,
)
    active_idx = @index(Global)
    @inbounds begin
        idx = inactive_indices[active_idx]
        fields.N[idx] = 0
        for layer_index in 1:Ntot
            fields.mass[layer_index, idx] = zero(density_init)
            fields.mass_w[layer_index, idx] = zero(density_init)
            fields.density[layer_index, idx] = density_init
            fields.temperature[layer_index, idx] = temperature_init
        end
        fields.mass_base[idx] = zero(density_init)
        fields.smb_ice[idx] = zero(density_init)
        fields.runoff[idx] = zero(density_init)
        fields.melt[idx] = zero(density_init)
        fields.refreezing[idx] = zero(density_init)
        fields.vapor_mass[idx] = zero(density_init)
        fields.sublimation[idx] = zero(density_init)
        fields.latent_heat_flux_sum[idx] = zero(density_init)
        fields.Tsrf[idx] = surface_temperature_init
        fields.albedo[idx] = albedo_init
        fields.snow_age_days[idx] = zero(density_init)
    end
end

function _reset_model_columns!(model::BESSIModel, ::BESSIState, runtime, inactive_indices::Vector{Int})
    isempty(inactive_indices) && return nothing
    backend_indices = _backend_active_indices(inactive_indices, runtime)
    kernel! = _reset_bessi_columns_kernel!(_ka_backend(runtime.state.mass))
    kernel!(
        get_fields(runtime.state),
        runtime.state.Ntot,
        backend_indices,
        convert(eltype(runtime.state.mass), model.density_init),
        convert(eltype(runtime.state.mass), model.temperature_init),
        runtime.state.c.T0,
        _initial_snow_albedo(runtime.state.c);
        ndrange=length(inactive_indices),
    )
    return nothing
end

@kernel function _reset_pdd_columns_kernel!(
    fields,
    inactive_indices,
)
    active_idx = @index(Global)
    @inbounds begin
        idx = inactive_indices[active_idx]
        fields.snowpack_swe[idx] = zero(eltype(fields.snowpack_swe))
        fields.smb_ice[idx] = zero(eltype(fields.smb_ice))
        fields.runoff[idx] = zero(eltype(fields.runoff))
        fields.pdd_sum[idx] = zero(eltype(fields.pdd_sum))
    end
end

function _reset_model_columns!(::PDDModel, ::PDDState, runtime, inactive_indices::Vector{Int})
    isempty(inactive_indices) && return nothing
    backend_indices = _backend_active_indices(inactive_indices, runtime)
    kernel! = _reset_pdd_columns_kernel!(_ka_backend(runtime.state.snowpack_swe))
    kernel!(
        get_fields(runtime.state),
        backend_indices;
        ndrange=length(inactive_indices),
    )
    return nothing
end

@kernel function _reset_itm_columns_kernel!(fields, inactive_indices, H_snow_init, albedo_init, T0)
    active_idx = @index(Global)
    @inbounds begin
        idx = inactive_indices[active_idx]
        zero_value = zero(H_snow_init)
        fields.H_snow[idx] = H_snow_init
        fields.alb_s[idx] = albedo_init
        fields.smb[idx] = zero_value
        fields.smbi[idx] = zero_value
        fields.melt[idx] = zero_value
        fields.runoff[idx] = zero_value
        fields.refreezing[idx] = zero_value
        fields.Tsrf[idx] = T0
        fields.melt_net[idx] = zero_value
        fields.smb_cum[idx] = zero_value
        fields.smb_ice[idx] = zero_value
        fields.melt_cum[idx] = zero_value
        fields.runoff_cum[idx] = zero_value
        fields.refreezing_cum[idx] = zero_value
    end
end

function _reset_model_columns!(model::ITMModel, ::ITMState, runtime, inactive_indices::Vector{Int})
    isempty(inactive_indices) && return nothing
    indices = _backend_active_indices(inactive_indices, runtime)
    _reset_itm_columns_kernel!(_ka_backend(runtime.state.H_snow))(get_fields(runtime.state), indices,
        model.H_snow_max, model.alb_snow_dry, model.c.T0; ndrange=length(inactive_indices))
    return nothing
end

function _set_model_runtime_active_indices!(model_runtime::ModelRuntime, active::AbstractVector{Bool})
    length(active) == length(model_runtime.active) || error("Active mask length must match the model column count.")
    active_v = Vector{Bool}(active)
    active_indices = findall(active_v)
    isempty(active_indices) && error("Active mask kept no Chion columns.")
    model_runtime.active = active_v
    model_runtime.active_indices = _backend_active_indices(active_indices, model_runtime.data)
    return model_runtime
end

function init_model_runtime!(sim, options::RunOptions, timings::StepTimingStats)
    data = prepare_runtime!(sim.model, sim.now, sim.forcing, options, timings)
    ncol = ncols(sim.model.grid)
    active = fill(true, ncol)
    active_indices = _backend_active_indices(collect(1:ncol), data)
    return ModelRuntime(data, active, active_indices)
end

function step_model!(model::BESSIModel, state::BESSIState, model_runtime::ModelRuntime, forcing, time_index::Int)
    return step_model!(model, state, model_runtime, forcing, time_index:time_index)
end

function step_model!(model::BESSIModel, ::BESSIState, model_runtime::ModelRuntime, forcing, time_range)
    runtime = model_runtime.data
    return _step_range!(
        runtime.state,
        forcing,
        time_range,
        runtime.workspace,
        model_runtime.active_indices,
        _bessi_step_kwargs(model),
    )
end

function step_model!(model::PDDModel, ::PDDState, model_runtime::ModelRuntime, forcing, time_index::Int)
    runtime = model_runtime.data
    return _pdd_step_arrays!(
        runtime.state,
        forcing,
        time_index,
        model,
        model_runtime.active_indices,
    )
end

function step_model!(model::ITMModel, ::ITMState, model_runtime::ModelRuntime, forcing, time_index::Int)
    _itm_step_arrays!(model_runtime.data.state, forcing, time_index, model, model_runtime.active_indices)
    return nothing
end

# PDD and ITM step one time index at a time.
function step_model!(model::Union{PDDModel, ITMModel}, state, model_runtime::ModelRuntime, forcing, time_range::AbstractUnitRange)
    for k in time_range
        step_model!(model, state, model_runtime, forcing, Int(k))
    end
    return nothing
end

finalize_state!(::BESSIModel, state::BESSIState, runtime, ::RunOptions, timings::StepTimingStats) =
    _sync_bessi_state!(state, runtime, timings)

"""Copy the backend state into the host `state` and refresh its diagnostics."""
function _sync_bessi_state!(state::BESSIState, runtime, timings::StepTimingStats)
    if runtime.is_gpu
        time_block!(timings, :gpu_transfer) do
            _copy_bessi_state!(state, runtime.state)
        end
    end
    update_diagnostics!(state)
    return nothing
end

function _copy_vector_state_from_runtime!(state, runtime)
    for name in fieldnames(typeof(state))
        copyto!(getfield(state, name), Array(getfield(runtime.state, name)))
    end
    return state
end

function finalize_state!(::Union{PDDModel,ITMModel}, state, runtime, ::RunOptions, timings::StepTimingStats)
    runtime.is_gpu || return nothing
    time_block!(timings, :gpu_transfer) do
        _copy_vector_state_from_runtime!(state, runtime)
    end
    return nothing
end
