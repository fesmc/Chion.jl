"""Model runtime preparation, backend transfer, and state finalization."""

_model_grid(model::AbstractSnowModel) = model.grid
_model_column_count(model::AbstractSnowModel) = ncols(_model_grid(model))
model_layer_count(::AbstractSnowModel, ::AbstractState, runtime) = 0
model_layer_count(::BESSIModel, ::CurrentState, runtime) = runtime.state.Ntot

model_output_groups(::BESSIModel) = (:final, :layers, :history, :monthly, :step, :daily)
model_output_groups(::PDDModel) = (:final, :history, :step)
model_output_groups(::ITMModel) = ()

function _validate_model_outputs!(model::AbstractSnowModel, options::RunOptions)
    options.write_netcdf || return nothing
    allowed = copy(NETCDF_VARIABLES)
    append!(allowed, keys(STATE_OUTPUT_ALIASES))
    append!(allowed, (:all, :none, :final, :layers, :monthly, :step, :daily, :history))
    unsupported = setdiff(options.netcdf_variables, allowed)
    isempty(unsupported) || error("Unsupported NetCDF variables: $(join(string.(unsupported), ", ")).")
    return nothing
end

function validate_integrator_setup!(sim, options::RunOptions)
    model = sim.model
    forcing = sim.forcing
    grid = sim.domain
    ncol = _model_column_count(model)
    size(forcing.air_temperature, 1) == ncol || error("Forcing column count must match the model column count.")
    if model isa BESSIModel && model.diurnal_shortwave_substeps
        all(isfinite, forcing.latitude_deg) || error("`latitude_deg` is required in `SnowpackForcing` when BESSI diurnal shortwave options are enabled.")
    end
    if model isa BESSIModel && _uses_prescribed_albedo(model.c)
        all(forcing.has_prescribed_albedo) || error("`prescribed_albedo` is required for every column and timestep when BESSI uses `PrescribedAlbedo`.")
    end
    spatial_grid = has_spatial_coords(grid)
    options.write_netcdf && !spatial_grid && error("NetCDF output requires a grid with spatial coordinates.")
    spatial_grid && length(grid.js) != ncol && error("Grid point count must match the domain column count.")
    _validate_model_outputs!(model, options)
    return nothing
end

function _prepare_backend!(timings::StepTimingStats, state::CurrentState, forcing::SnowpackForcing; is_gpu::Bool)
    step_fields = forcing
    if is_gpu
        cuda_available() || error("`backend=gpu` requested, but CUDA is not functional in the current environment.")
        state = time_block!(timings, :gpu_transfer) do
            gpu_state(state)
        end
        step_fields = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), forcing)
        end
        workspace = time_block!(timings, :gpu_transfer) do
            ColumnarStepWorkspace(state)
        end
        return (state=state, step_fields=step_fields, workspace=workspace, is_gpu=true)
    end
    workspace = time_block!(timings, :create_workspaces) do
        ColumnarStepWorkspace(state)
    end
    return (state=state, step_fields=step_fields, workspace=workspace, is_gpu=false)
end

function prepare_runtime!(model::BESSIModel, state::CurrentState, forcing::SnowpackForcing, options::RunOptions, timings::StepTimingStats)
    return _prepare_backend!(timings, state, forcing; is_gpu=options.backend == :gpu)
end

function prepare_runtime!(model::PDDModel, state::PDDState, forcing::SnowpackForcing, options::RunOptions, timings::StepTimingStats)
    is_gpu = options.backend == :gpu
    if is_gpu
        cuda_available() || error("`backend=gpu` requested, but CUDA is not functional in the current environment.")
        snowpack_swe = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), state.snowpack_swe)
        end
        smb_ice = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), state.smb_ice)
        end
        runoff = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), state.runoff)
        end
        pdd_sum = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), state.pdd_sum)
        end
        step_fields = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), forcing)
        end
        return (snowpack_swe=snowpack_swe, smb_ice=smb_ice, runoff=runoff, pdd_sum=pdd_sum, step_fields=step_fields, is_gpu=true)
    end
    ncol = ncols(model.grid)
    scratch = (
        a=Vector{Float64}(undef, ncol),
        b=Vector{Float64}(undef, ncol),
        c=Vector{Float64}(undef, ncol),
        d=Vector{Float64}(undef, ncol),
        e=Vector{Float64}(undef, ncol),
        f=Vector{Float64}(undef, ncol),
    )
    return (snowpack_swe=state.snowpack_swe, smb_ice=state.smb_ice, runoff=state.runoff, pdd_sum=state.pdd_sum, step_fields=forcing, scratch=scratch, is_gpu=false)
end

function prepare_runtime!(::ITMModel, ::AbstractState, ::SnowpackForcing, ::RunOptions, ::StepTimingStats)
    error("ITMModel is not yet implemented. Physics coming soon.")
end

_backend_active_indices(indices::Vector{Int}, backend) =
    getproperty(backend, :is_gpu) ? adapt(gpu_storage_type(), indices) : indices

@kernel function _reset_bessi_columns_kernel!(
    N,
    mass,
    mass_w,
    density,
    temperature,
    mass_base,
    smb_ice,
    runoff,
    melt,
    refreezing,
    vapor_mass,
    sublimation,
    latent_heat_flux_sum,
    Tsrf,
    snow_cover,
    albedo_dynamic,
    inactive_indices,
    Ntot::Int,
    density_init,
    temperature_init,
    surface_temperature_init,
    albedo_init,
)
    active_idx = @index(Global)
    if active_idx <= length(inactive_indices)
        idx = inactive_indices[active_idx]
        N[idx] = 0
        for layer_index in 1:Ntot
            mass[layer_index, idx] = zero(density_init)
            mass_w[layer_index, idx] = zero(density_init)
            density[layer_index, idx] = density_init
            temperature[layer_index, idx] = temperature_init
        end
        mass_base[idx] = zero(density_init)
        smb_ice[idx] = zero(density_init)
        runoff[idx] = zero(density_init)
        melt[idx] = zero(density_init)
        refreezing[idx] = zero(density_init)
        vapor_mass[idx] = zero(density_init)
        sublimation[idx] = zero(density_init)
        latent_heat_flux_sum[idx] = zero(density_init)
        Tsrf[idx] = surface_temperature_init
        snow_cover[idx] = zero(density_init)
        albedo_dynamic[idx] = albedo_init
    end
end

function _reset_model_columns!(model::BESSIModel, ::CurrentState, runtime, inactive_indices::Vector{Int})
    isempty(inactive_indices) && return nothing
    backend_indices = _backend_active_indices(inactive_indices, runtime)
    kernel! = _reset_bessi_columns_kernel!(_ka_backend(runtime.state.mass))
    event = kernel!(
        runtime.state.N,
        runtime.state.mass,
        runtime.state.mass_w,
        runtime.state.density,
        runtime.state.temperature,
        runtime.state.mass_base,
        runtime.state.smb_ice,
        runtime.state.runoff,
        runtime.state.melt,
        runtime.state.refreezing,
        runtime.state.vapor_mass,
        runtime.state.sublimation,
        runtime.state.latent_heat_flux_sum,
        runtime.state.Tsrf,
        runtime.state.snow_cover,
        runtime.state.albedo_dynamic,
        backend_indices,
        runtime.state.Ntot,
        convert(eltype(runtime.state.mass), model.density_init),
        convert(eltype(runtime.state.mass), model.temperature_init),
        runtime.state.c.T0,
        runtime.state.c.alpha_dry;
        ndrange=length(inactive_indices),
    )
    _wait_kernel(event)
    return nothing
end

@kernel function _reset_pdd_columns_kernel!(
    snowpack_swe,
    smb_ice,
    runoff,
    pdd_sum,
    inactive_indices,
)
    active_idx = @index(Global)
    if active_idx <= length(inactive_indices)
        idx = inactive_indices[active_idx]
        snowpack_swe[idx] = zero(eltype(snowpack_swe))
        smb_ice[idx] = zero(eltype(smb_ice))
        runoff[idx] = zero(eltype(runoff))
        pdd_sum[idx] = zero(eltype(pdd_sum))
    end
end

function _reset_model_columns!(::PDDModel, ::PDDState, runtime, inactive_indices::Vector{Int})
    isempty(inactive_indices) && return nothing
    backend_indices = _backend_active_indices(inactive_indices, runtime)
    kernel! = _reset_pdd_columns_kernel!(_ka_backend(runtime.snowpack_swe))
    event = kernel!(
        runtime.snowpack_swe,
        runtime.smb_ice,
        runtime.runoff,
        runtime.pdd_sum,
        backend_indices;
        ndrange=length(inactive_indices),
    )
    _wait_kernel(event)
    return nothing
end

_reset_model_columns!(::AbstractSnowModel, ::AbstractState, runtime, inactive_indices::Vector{Int}) = nothing

function _set_model_runtime_active_indices!(model_runtime::ModelRuntime, active::AbstractVector{Bool})
    length(active) == model_runtime.ncol || error("Active mask length must match the model column count.")
    active_v = Vector{Bool}(active)
    active_indices = findall(active_v)
    isempty(active_indices) && error("Active mask kept no Chion columns.")
    model_runtime.active = active_v
    model_runtime.active_indices = _backend_active_indices(active_indices, model_runtime.backend)
    return model_runtime
end

function init_model_runtime!(sim, options::RunOptions, timings::StepTimingStats)
    backend = prepare_runtime!(sim.model, sim.now, sim.forcing, options, timings)
    ncol = _model_column_count(sim.model)
    active = trues(ncol)
    active_indices = _backend_active_indices(collect(1:ncol), backend)
    return ModelRuntime(backend, ncol, sim.domain, active, active_indices)
end

function step_model!(model::BESSIModel, ::CurrentState, model_runtime::ModelRuntime, forcing::SnowpackForcing, time_index::Int)
    runtime = model_runtime.backend
    step!(
        runtime.state,
        forcing,
        time_index,
        runtime.workspace,
        model_runtime.active_indices;
        diurnal_shortwave_substeps=model.diurnal_shortwave_substeps,
        diurnal_shortwave_threshold=model.diurnal_shortwave_threshold,
        diurnal_shortwave_max_substeps=model.diurnal_shortwave_max_substeps,
        diurnal_shortwave_min_air_temperature=model.diurnal_shortwave_min_air_temperature,
        diurnal_temperature_cycle=model.diurnal_temperature_cycle,
        diurnal_temperature_amplitude=model.diurnal_temperature_amplitude,
    )
    return nothing
end

function step_model!(model::PDDModel, ::PDDState, model_runtime::ModelRuntime, forcing::SnowpackForcing, time_index::Int)
    runtime = model_runtime.backend
    if runtime.snowpack_swe isa Vector{Float64}
        if forcing_step_kind(forcing, time_index) === :monthly
            pdd_monthly_step!(
                runtime.snowpack_swe,
                runtime.smb_ice,
                runtime.runoff,
                runtime.pdd_sum,
                forcing,
                time_index,
                model,
                runtime.scratch,
            )
        else
            pdd_step!(
                runtime.snowpack_swe,
                runtime.smb_ice,
                runtime.runoff,
                runtime.pdd_sum,
                forcing,
                time_index,
                model.ddf_snow,
                model.ddf_ice,
                model.refreezing_fraction,
                runtime.scratch,
            )
        end
    else
        if forcing_step_kind(forcing, time_index) === :monthly
            pdd_monthly_step!(
                runtime.snowpack_swe,
                runtime.smb_ice,
                runtime.runoff,
                runtime.pdd_sum,
                forcing,
                time_index,
                model,
                model_runtime.active_indices,
            )
        else
            pdd_step!(
                runtime.snowpack_swe,
                runtime.smb_ice,
                runtime.runoff,
                runtime.pdd_sum,
                forcing,
                time_index,
                model.ddf_snow,
                model.ddf_ice,
                model.refreezing_fraction,
                model_runtime.active_indices,
            )
        end
    end
    return nothing
end

finalize_state!(::AbstractSnowModel, ::AbstractState, runtime, ::RunOptions, ::StepTimingStats) = nothing

function finalize_state!(::BESSIModel, state::CurrentState, runtime, options::RunOptions, timings::StepTimingStats)
    if runtime.is_gpu
        time_block!(timings, :gpu_transfer) do
            _copy_current_state!(state, cpu_state(runtime.state))
        end
    end
    update_diagnostics!(state)
    return nothing
end

function finalize_state!(::PDDModel, state::PDDState, runtime, options::RunOptions, timings::StepTimingStats)
    runtime.is_gpu || return nothing
    time_block!(timings, :gpu_transfer) do
        copyto!(state.snowpack_swe, Array(runtime.snowpack_swe))
        copyto!(state.smb_ice, Array(runtime.smb_ice))
        copyto!(state.runoff, Array(runtime.runoff))
        copyto!(state.pdd_sum, Array(runtime.pdd_sum))
    end
    return nothing
end
