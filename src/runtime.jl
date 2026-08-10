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
        all(isfinite, forcing.latitude_deg) || error("`latitude_deg` is required in `SnowpackForcing` when BESSI diurnal shortwave options are enabled.")
    end
    if _uses_prescribed_albedo(model.c)
        all(forcing.has_prescribed_albedo) || error("`prescribed_albedo` is required for every column and timestep when BESSI uses `albedo=:prescribed`.")
    end
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

function _prepare_backend!(timings::StepTimingStats, state::BESSIState, forcing::SnowpackForcing; is_gpu::Bool)
    step_fields = get_fields(forcing)
    if is_gpu
        cuda_available() || error("`backend=gpu` requested, but CUDA is not functional in the current environment.")
        state = time_block!(timings, :gpu_transfer) do
            gpu_state(state)
        end
        step_fields = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), step_fields)
        end
        workspace = time_block!(timings, :gpu_transfer) do
            ColumnarStepWorkspace(state)
        end
        return ModelRuntimeData(state, step_fields, workspace, true)
    end
    workspace = time_block!(timings, :create_workspaces) do
        ColumnarStepWorkspace(state)
    end
    return ModelRuntimeData(state, step_fields, workspace, false)
end

function prepare_runtime!(model::BESSIModel, state::BESSIState, forcing::SnowpackForcing, options::RunOptions, timings::StepTimingStats)
    return _prepare_backend!(timings, state, forcing; is_gpu=options.backend == :gpu)
end

function _gpu_forcing_view(forcing::PDDForcing)
    to_gpu(field) = adapt(gpu_storage_type(), field)
    return PDDForcing(
        forcing.dt_days,
        to_gpu(forcing.air_temperature),
        to_gpu(forcing.snowfall_rate),
        to_gpu(forcing.rainfall_rate),
    )
end

function _gpu_forcing_view(forcing::ITMForcing)
    to_gpu(field) = adapt(gpu_storage_type(), field)
    return ITMForcing(
        forcing.dt_days,
        to_gpu(forcing.air_temperature),
        to_gpu(forcing.snowfall_rate),
        to_gpu(forcing.rainfall_rate),
        to_gpu(forcing.shortwave_down),
        to_gpu(forcing.q_sw_net),
        to_gpu(forcing.has_q_sw_net),
        to_gpu(forcing.latitude_deg),
        to_gpu(forcing.surface_height),
        to_gpu(forcing.ice_thickness),
        to_gpu(forcing.annual_pdd),
    )
end

function prepare_runtime!(model::PDDModel, state::PDDState, forcing::SnowpackForcing, options::RunOptions, timings::StepTimingStats)
    is_gpu = options.backend == :gpu
    step_fields = PDDForcing(forcing)
    if is_gpu
        cuda_available() || error("`backend=gpu` requested, but CUDA is not functional in the current environment.")
        backend_state = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), state)
        end
        step_fields = time_block!(timings, :gpu_transfer) do
            _gpu_forcing_view(step_fields)
        end
        return ModelRuntimeData(backend_state, step_fields, nothing, true)
    end
    return ModelRuntimeData(state, step_fields, nothing, false)
end

function prepare_runtime!(model::ITMModel, state::ITMState, forcing::SnowpackForcing, options::RunOptions, timings::StepTimingStats)
    is_gpu = options.backend == :gpu
    is_gpu && !cuda_available() && error("`backend=gpu` requested, but CUDA is not functional in the current environment.")
    backend_state = is_gpu ? time_block!(timings, :gpu_transfer) do
        adapt(gpu_storage_type(), state)
    end : state
    step_fields = ITMForcing(forcing)
    step_fields = is_gpu ? time_block!(timings, :gpu_transfer) do
        _gpu_forcing_view(step_fields)
    end : step_fields
    return ModelRuntimeData(backend_state, step_fields, nothing, is_gpu)
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
    event = kernel!(
        get_fields(runtime.state),
        runtime.state.Ntot,
        backend_indices,
        convert(eltype(runtime.state.mass), model.density_init),
        convert(eltype(runtime.state.mass), model.temperature_init),
        runtime.state.c.T0,
        _initial_snow_albedo(runtime.state.c);
        ndrange=length(inactive_indices),
    )
    _wait_kernel(event)
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
    event = kernel!(
        get_fields(runtime.state),
        backend_indices;
        ndrange=length(inactive_indices),
    )
    _wait_kernel(event)
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
    event = _reset_itm_columns_kernel!(_ka_backend(runtime.state.H_snow))(get_fields(runtime.state), indices,
        model.H_snow_max, model.alb_snow_dry, model.c.T0; ndrange=length(inactive_indices))
    _wait_kernel(event)
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

function finalize_state!(::BESSIModel, state::BESSIState, runtime, ::RunOptions, timings::StepTimingStats)
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

function finalize_state!(::PDDModel, state::PDDState, runtime, ::RunOptions, timings::StepTimingStats)
    runtime.is_gpu || return nothing
    time_block!(timings, :gpu_transfer) do
        _copy_vector_state_from_runtime!(state, runtime)
    end
    return nothing
end

function finalize_state!(::ITMModel, state::ITMState, runtime, ::RunOptions, timings::StepTimingStats)
    runtime.is_gpu || return nothing
    time_block!(timings, :gpu_transfer) do
        _copy_vector_state_from_runtime!(state, runtime)
    end
    return nothing
end
