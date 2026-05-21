"""Model runtime preparation, backend transfer, and state finalization."""

_model_grid(model::AbstractSnowModel) = model.grid
_model_column_count(model::AbstractSnowModel) = ncols(_model_grid(model))
model_layer_count(::AbstractSnowModel, ::AbstractSnowModelState, runtime) = 0
model_layer_count(::BESSIModel, ::BESSIState, runtime) = runtime.domain.Ntot

model_output_groups(::BESSIModel) = (:final, :layers, :history, :monthly, :step)
model_output_groups(::PDDModel) = (:final, :history, :step)
model_output_groups(::ITMModel) = ()

function _validate_model_outputs!(model::AbstractSnowModel, options::RunOptions)
    options.write_netcdf || return nothing
    supported = model_output_groups(model)
    allowed = Symbol[]
    for group in supported
        append!(allowed, getproperty(OUTPUT_GROUPS, group))
    end
    unsupported = setdiff(options.netcdf_variables, allowed)
    isempty(unsupported) || error("$(typeof(model)) output supports only $(join(string.(supported), ", ")) NetCDF groups; unsupported: $(join(string.(unsupported), ", ")).")
    return nothing
end

function validate_integrator_setup!(sim, options::RunOptions)
    model = sim.model
    forcing = sim.forcing
    grid = _model_grid(model)
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

function _prepare_backend!(timings::StepTimingStats, domain::SnowpackDomain, forcing::SnowpackForcing; is_gpu::Bool)
    step_fields = forcing
    if is_gpu
        cuda_available() || error("`backend=gpu` requested, but CUDA is not functional in the current environment.")
        domain = time_block!(timings, :gpu_transfer) do
            gpu_domain(domain)
        end
        step_fields = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), forcing)
        end
        workspace = time_block!(timings, :gpu_transfer) do
            ColumnarStepWorkspace(domain)
        end
        return (domain=domain, step_fields=step_fields, workspace=workspace, is_gpu=true)
    end
    workspace = time_block!(timings, :create_workspaces) do
        ColumnarStepWorkspace(domain)
    end
    return (domain=domain, step_fields=step_fields, workspace=workspace, is_gpu=false)
end

function prepare_runtime!(model::BESSIModel, state::BESSIState, forcing::SnowpackForcing, options::RunOptions, timings::StepTimingStats)
    return _prepare_backend!(timings, state.domain, forcing; is_gpu=options.backend == :gpu)
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

function prepare_runtime!(::ITMModel, ::AbstractSnowModelState, ::SnowpackForcing, ::RunOptions, ::StepTimingStats)
    error("ITMModel is not yet implemented. Physics coming soon.")
end

function init_model_runtime!(sim, options::RunOptions, timings::StepTimingStats)
    backend = prepare_runtime!(sim.model, sim.now, sim.forcing, options, timings)
    return ModelRuntime(backend, _model_column_count(sim.model), _model_grid(sim.model))
end

function step_model!(model::BESSIModel, ::BESSIState, runtime, forcing::SnowpackForcing, time_index::Int)
    step!(
        runtime.domain,
        forcing,
        time_index,
        runtime.workspace;
        diurnal_shortwave_substeps=model.diurnal_shortwave_substeps,
        diurnal_shortwave_threshold=model.diurnal_shortwave_threshold,
        diurnal_shortwave_max_substeps=model.diurnal_shortwave_max_substeps,
        diurnal_shortwave_min_air_temperature=model.diurnal_shortwave_min_air_temperature,
        diurnal_temperature_cycle=model.diurnal_temperature_cycle,
        diurnal_temperature_amplitude=model.diurnal_temperature_amplitude,
    )
    return nothing
end

function step_model!(model::PDDModel, ::PDDState, runtime, forcing::SnowpackForcing, time_index::Int)
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
            )
        end
    end
    return nothing
end

finalize_state!(::AbstractSnowModel, ::AbstractSnowModelState, runtime, ::RunOptions, ::StepTimingStats) = nothing

function finalize_state!(::BESSIModel, state::BESSIState, runtime, options::RunOptions, timings::StepTimingStats)
    runtime.is_gpu || return nothing
    time_block!(timings, :gpu_transfer) do
        _copy_domain_state!(state.domain, cpu_domain(runtime.domain))
    end
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
