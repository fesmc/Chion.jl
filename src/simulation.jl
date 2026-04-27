"""
Internal runtime types and helpers used by execute_model_run!.
These are not part of the public API.
"""

struct RunOptions
    name::String
    input_label::String
    output_dir::String
    netcdf_path::String
    write_outputs::Bool
    write_netcdf::Bool
    netcdf_variables::Vector{Symbol}
    cycles::Int
    backend::Symbol
    history_stride::Int
end

@inline function normalize_backend(backend)
    value = lowercase(strip(String(backend)))
    value == "cpu" && return :threads
    value in ("threads", "gpu") || error("Unsupported backend '$backend'. Use `threads`, `cpu`, or `gpu`.")
    return Symbol(value)
end

@inline normalize_history_stride(stride::Integer) =
    Int(stride) >= 0 ? Int(stride) : error("`history_stride` must be >= 0.")

@inline should_record_cycle_metrics(cycle::Int, cycles::Int, stride::Int) =
    cycle == cycles || (stride > 0 && mod(cycle, stride) == 0)

function RunOptions(;
    name::AbstractString="chion_run",
    input_label::AbstractString="",
    output_dir::AbstractString="",
    netcdf_path::AbstractString="",
    write_outputs::Bool=true,
    write_netcdf::Bool=true,
    netcdf_variables=copy(NETCDF_VARIABLES),
    cycles::Integer=10,
    backend=:threads,
    history_stride::Integer=1,
)
    resolved_name = String(name)
    resolved_output_dir = isempty(output_dir) ? _default_output_dir(resolved_name) : String(output_dir)
    return RunOptions(
        resolved_name,
        String(input_label),
        resolved_output_dir,
        String(netcdf_path),
        write_outputs,
        write_netcdf,
        normalize_netcdf_variables(netcdf_variables),
        Int(cycles),
        normalize_backend(backend),
        normalize_history_stride(history_stride),
    )
end

include("reporting.jl")

function _prepare_backend!(timings::StepTimingStats, domain::SnowpackDomain, forcing::SnowpackForcing; is_gpu::Bool)
    step_fields = forcing
    if is_gpu
        cuda_available() || error("`backend=gpu` requested, but CUDA is not functional in the current environment.")
        domain = time_block!(timings, :gpu_transfer) do
            gpu_domain(domain)
        end
        step_fields = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), step_fields)
        end
        workspace = time_block!(timings, :gpu_transfer) do
            ColumnarStepWorkspace(domain)
        end
        return domain, step_fields, workspace
    end
    workspace = time_block!(timings, :create_workspaces) do
        ColumnarStepWorkspace(domain)
    end
    return domain, step_fields, workspace
end

_model_grid(model::AbstractSnowModel) = model.grid
_model_column_count(model::AbstractSnowModel) = ncols(_model_grid(model))
_model_layer_count(::AbstractSnowModel, runtime) = 0
_model_layer_count(::BESSIModel, runtime) = runtime.domain.Ntot

_model_supported_output_groups(::BESSIModel) = (:final, :layers, :history, :monthly, :step)
_model_supported_output_groups(::PDDModel) = (:final, :history, :step)
_model_supported_output_groups(::ITMModel) = ()

function _validate_model_outputs!(model::AbstractSnowModel, options::RunOptions)
    options.write_netcdf || return nothing
    supported = _model_supported_output_groups(model)
    allowed = Symbol[]
    for group in supported
        append!(allowed, getproperty(OUTPUT_GROUPS, group))
    end
    unsupported = setdiff(options.netcdf_variables, allowed)
    isempty(unsupported) || error("$(typeof(model)) output supports only $(join(string.(supported), ", ")) NetCDF groups; unsupported: $(join(string.(unsupported), ", ")).")
    return nothing
end

function _prepare_model_runtime!(model::BESSIModel, forcing::SnowpackForcing, options::RunOptions, timings::StepTimingStats)
    is_gpu = options.backend == :gpu
    domain, step_fields, workspace = _prepare_backend!(timings, model.domain, forcing; is_gpu=is_gpu)
    return (domain=domain, step_fields=step_fields, workspace=workspace, is_gpu=is_gpu)
end

function _prepare_model_runtime!(model::PDDModel, forcing::SnowpackForcing, options::RunOptions, timings::StepTimingStats)
    is_gpu = options.backend == :gpu
    if is_gpu
        cuda_available() || error("`backend=gpu` requested, but CUDA is not functional in the current environment.")
        snowpack_swe = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), model.snowpack_swe)
        end
        smb_ice = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), model.smb_ice)
        end
        runoff = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), model.runoff)
        end
        pdd_sum = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), model.pdd_sum)
        end
        step_fields = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), forcing)
        end
        return (snowpack_swe=snowpack_swe, smb_ice=smb_ice, runoff=runoff, pdd_sum=pdd_sum, step_fields=step_fields, is_gpu=true)
    end
    ncol = ncols(model.grid)
    scratch = (a=Vector{Float64}(undef, ncol), b=Vector{Float64}(undef, ncol), c=Vector{Float64}(undef, ncol), d=Vector{Float64}(undef, ncol), e=Vector{Float64}(undef, ncol), f=Vector{Float64}(undef, ncol))
    return (snowpack_swe=model.snowpack_swe, smb_ice=model.smb_ice, runoff=model.runoff, pdd_sum=model.pdd_sum, step_fields=forcing, scratch=scratch, is_gpu=false)
end

function _prepare_model_runtime!(::ITMModel, ::SnowpackForcing, ::RunOptions, ::StepTimingStats)
    error("ITMModel is not yet implemented. Physics coming soon.")
end

function _allocate_cycle_backend_buffers(model::BESSIModel, runtime, ncol::Int)
    return allocate_cycle_summary_buffers(runtime.domain, ncol)
end
function _allocate_cycle_backend_buffers(::AbstractSnowModel, runtime, ncol::Int)
    return allocate_cycle_summary_buffers(ncol)
end
function _allocate_cycle_backend_buffers(::PDDModel, runtime, ncol::Int)
    return _named_buffers(CYCLE_BUFFER_NAMES, () -> similar(runtime.snowpack_swe, Float64, ncol))
end
function _allocate_step_backend_buffers(model::BESSIModel, runtime, ncol::Int)
    return allocate_summary_buffers(runtime.domain, ncol)
end
function _allocate_step_backend_buffers(::AbstractSnowModel, runtime, ncol::Int)
    return allocate_summary_buffers(ncol)
end
function _allocate_step_backend_buffers(::PDDModel, runtime, ncol::Int)
    return _named_buffers(SUMMARY_BUFFER_NAMES, () -> similar(runtime.snowpack_swe, Float64, ncol))
end

function _summarize_cycle_state!(summary, ::BESSIModel, runtime)
    summarize_cycle_state!(
        summary.thickness,
        summary.wet_mass,
        summary.bulk_density,
        summary.base_mass,
        runtime.domain,
    )
    return summary
end

function _summarize_cycle_state!(summary, model::PDDModel, runtime)
    summary.thickness .= runtime.snowpack_swe ./ 1000.0
    summary.wet_mass .= runtime.snowpack_swe
    fill!(summary.bulk_density, NaN)
    summary.base_mass .= runtime.smb_ice
    return summary
end

function _summarize_step_state!(summary, ::PDDModel, runtime)
    summary.thickness .= runtime.snowpack_swe ./ 1000.0
    summary.wet_mass .= runtime.snowpack_swe
    fill!(summary.bulk_density, NaN)
    summary.base_mass .= runtime.smb_ice
    summary.smb_ice .= runtime.smb_ice
    fill!(summary.liquid_water, 0.0)
    summary.runoff .= runtime.runoff
    summary.pdd .= runtime.pdd_sum
    return summary
end

function _summarize_step_state!(summary, ::BESSIModel, runtime)
    summarize_domain_state!(
        summary.thickness,
        summary.wet_mass,
        summary.bulk_density,
        summary.base_mass,
        summary.smb_ice,
        summary.liquid_water,
        summary.runoff,
        runtime.domain,
    )
    fill!(summary.pdd, 0.0)
    return summary
end

function _model_step!(::BESSIModel, runtime, forcing::SnowpackForcing, time_index::Int)
    step!(runtime.domain, runtime.step_fields, time_index, runtime.workspace)
    return nothing
end

function _model_step!(model::PDDModel, runtime, forcing::SnowpackForcing, time_index::Int)
    if runtime.snowpack_swe isa Vector{Float64}
        pdd_step!(
            runtime.snowpack_swe,
            runtime.smb_ice,
            runtime.runoff,
            runtime.pdd_sum,
            runtime.step_fields,
            time_index,
            model.ddf_snow,
            model.ddf_ice,
            model.refreezing_fraction,
            runtime.scratch,
        )
    else
        pdd_step!(
            runtime.snowpack_swe,
            runtime.smb_ice,
            runtime.runoff,
            runtime.pdd_sum,
            runtime.step_fields,
            time_index,
            model.ddf_snow,
            model.ddf_ice,
            model.refreezing_fraction,
        )
    end
    return nothing
end

_model_supports_cycle_step(::AbstractSnowModel, runtime) = false
_model_supports_cycle_step(::PDDModel, runtime) = runtime.snowpack_swe isa Vector{Float64}

function _model_cycle_step!(model::PDDModel, runtime, forcing::SnowpackForcing)
    pdd_step!(
        runtime.snowpack_swe,
        runtime.smb_ice,
        runtime.runoff,
        runtime.pdd_sum,
        runtime.step_fields,
        model.ddf_snow,
        model.ddf_ice,
        model.refreezing_fraction,
        runtime.scratch,
    )
    return nothing
end

_model_smb_ice_vector(::BESSIModel, runtime) = _host_vector(runtime.domain.smb_ice; copy_array=true)
_model_smb_ice_vector(model::PDDModel, runtime) = _host_vector(runtime.smb_ice; copy_array=true)
_model_runoff_vector(::BESSIModel, runtime) = _host_vector(runtime.domain.runoff; copy_array=true)
_model_runoff_vector(model::PDDModel, runtime) = _host_vector(runtime.runoff; copy_array=true)
#_model_pdd_vector(model::AbstractSnowModel, runtime) = zeros(Float64, ncols(model.grid))
_model_pdd_vector(model::PDDModel, runtime) = _host_vector(runtime.pdd_sum; copy_array=true)

function _update_cycle_smb_delta!(last_delta::Vector{Float64}, previous_cycle_smb_ice::Vector{Float64}, model::AbstractSnowModel, runtime)
    current = _model_smb_ice_vector(model, runtime)
    last_delta .= current .- previous_cycle_smb_ice
    previous_cycle_smb_ice .= current
    return nothing
end

function _generic_final_output_vector(model::AbstractSnowModel, runtime, spec, final_state, deltas)
    source = spec.source
    values = source === :final_state ? getfield(final_state, spec.field) :
        source === :deltas ? getfield(deltas, spec.field) :
        source === :domain && spec.field === :smb_ice ? _model_smb_ice_vector(model, runtime) :
        source === :domain && spec.field === :runoff ? _model_runoff_vector(model, runtime) :
        error("Unsupported final output source `$(source)`.")
    return _host_vector(values; copy_array=true)
end

function _model_final_output_vector(model::AbstractSnowModel, runtime, spec, final_state, deltas)
    return _generic_final_output_vector(model, runtime, spec, final_state, deltas)
end

function _model_final_output_vector(model::PDDModel, runtime, spec, final_state, deltas)
    spec.key === :final_base_mass && return zeros(Float64, ncols(model.grid))
    return _generic_final_output_vector(model, runtime, spec, final_state, deltas)
end

function _scatter_model_final_grids(model::AbstractSnowModel, runtime, final_state, deltas, layout)
    return NamedTuple{FINAL_GRID_KEYS}(ntuple(i -> begin
        values = _model_final_output_vector(model, runtime, FINAL_OUTPUT_SPECS[i], final_state, deltas)
        scatter_to_grid(values, layout.js, layout.is, _grid_shape(layout))
    end, length(FINAL_OUTPUT_SPECS)))
end

function _model_layer_grids(model::AbstractSnowModel, runtime, layout, need_layer_outputs::Bool, timings::StepTimingStats)
    return _empty_layer_grids()
end

function _model_layer_grids(model::BESSIModel, runtime, layout, need_layer_outputs::Bool, timings::StepTimingStats)
    need_layer_outputs || return _empty_layer_grids()
    final_domain = runtime.is_gpu ? time_block!(timings, :gpu_transfer) do
        cpu_domain(runtime.domain)
    end : runtime.domain
    return time_block!(timings, :collect_final_layer_grids) do
        collect_final_layer_grids(final_domain, layout.js, layout.is, _grid_shape(layout), final_domain.Ntot)
    end
end

function _finalize_model_runtime!(model::AbstractSnowModel, runtime, options::RunOptions, timings::StepTimingStats)
    return nothing
end

function _finalize_model_runtime!(model::BESSIModel, runtime, options::RunOptions, timings::StepTimingStats)
    if options.backend == :gpu
        _copy_domain_state!(model.domain, cpu_domain(runtime.domain))
    end
    return nothing
end

function _finalize_model_runtime!(model::PDDModel, runtime, options::RunOptions, timings::StepTimingStats)
    runtime.is_gpu || return nothing
    time_block!(timings, :gpu_transfer) do
        copyto!(model.snowpack_swe, Array(runtime.snowpack_swe))
        copyto!(model.smb_ice, Array(runtime.smb_ice))
        copyto!(model.runoff, Array(runtime.runoff))
        copyto!(model.pdd_sum, Array(runtime.pdd_sum))
    end
    return nothing
end

function execute_model_run!(
    model::AbstractSnowModel,
    forcing::SnowpackForcing;
    options::RunOptions=RunOptions(),
    io::IO=stdout,
    timings::StepTimingStats=StepTimingStats(),
    run_wall_t0::Integer=time_ns(),
)
    grid = _model_grid(model)
    ncol = _model_column_count(model)
    size(forcing.air_temperature, 1) == ncol || error("Forcing column count must match the model column count.")
    spatial_grid = has_spatial_coords(grid)
    options.write_netcdf && !spatial_grid && error("NetCDF output requires a grid with spatial coordinates.")
    spatial_grid && length(grid.js) != ncol && error("Grid point count must match the domain column count.")
    _validate_model_outputs!(model, options)

    selected = Set(options.netcdf_variables)
    write_final_fields = options.write_netcdf
    need_step_outputs = options.write_netcdf && any(var -> var in selected, OUTPUT_GROUPS.step)
    need_monthly_outputs = options.write_netcdf && any(var -> var in selected, OUTPUT_GROUPS.monthly)
    need_layer_outputs = options.write_netcdf && any(var -> var in selected, OUTPUT_GROUPS.layers)
    need_last_cycle_smb_delta = options.write_netcdf && (:last_cycle_delta_ice_sheet_smb in selected)
    need_step_diagnostics = need_step_outputs || need_monthly_outputs

    runtime = _prepare_model_runtime!(model, forcing, options, timings)
    prev = allocate_cycle_summary_buffers(ncol)
    final = allocate_cycle_summary_buffers(ncol)
    backend_cycle_summary = _allocate_cycle_backend_buffers(model, runtime, ncol)
    time_block!(timings, :summarize_columns_initial) do
        _summarize_cycle_state!(backend_cycle_summary, model, runtime)
        _copy_summary_fields!(prev, backend_cycle_summary, CYCLE_BUFFER_NAMES)
    end
    initial_thickness_vec = write_final_fields ? copy(prev.thickness) : Float64[]

    schedule = nothing
    writer = nothing
    nc_path = ""
    if options.write_netcdf
        schedule = time_block!(timings, :prepare_output_schedule) do
            _prepare_output_schedule(forcing.time_values, options.cycles)
        end
        initial_thickness = time_block!(timings, :prepare_initial_output_fields_grid) do
            scatter_to_grid(initial_thickness_vec, grid.js, grid.is, _grid_shape(grid))
        end
        nc_path = resolve_netcdf_path(options)
        writer = time_block!(timings, :init_netcdf) do
            init_netcdf(
                nc_path,
                options,
                forcing.time_values,
                _model_layer_count(model, runtime),
                grid,
                initial_thickness,
                schedule.month_cycle,
                schedule.month_of_year,
                schedule.source_month_code,
                schedule.annual_output.source_indices,
                schedule.annual_output.source_codes,
            )
        end
    end
    write_step_fields = writer !== nothing && need_step_outputs

    step_summary = need_step_diagnostics ? allocate_summary_buffers(ncol) : nothing
    backend_step_summary = need_step_diagnostics ? _allocate_step_backend_buffers(model, runtime, ncol) : nothing
    deltas = (
        thickness=fill(NaN, ncol),
        wet_mass=fill(NaN, ncol),
        base_mass=fill(NaN, ncol),
        ice_sheet_smb=need_last_cycle_smb_delta ? fill(NaN, ncol) : Float64[],
    )
    monthly_total = need_monthly_outputs && schedule !== nothing ? schedule.nmonth_total : 0
    monthly_sums = _allocate_monthly_sums(need_monthly_outputs, monthly_total, ncol)
    monthly_count = need_monthly_outputs ? zeros(Int32, monthly_total) : Int32[]
    step_vectors = _allocate_step_vectors(need_step_outputs, ncol)

    previous = (
        base_mass=_host_vector(prev.base_mass; copy_array=true),
        smb_ice=_model_smb_ice_vector(model, runtime),
        runoff=_model_runoff_vector(model, runtime),
        pdd=_model_pdd_vector(model, runtime),
    )
    previous_cycle_smb_ice = need_last_cycle_smb_delta ? _model_smb_ice_vector(model, runtime) : Float64[]
    history = NamedTuple[]
    steps_written = 0

    use_cycle_step = _model_supports_cycle_step(model, runtime) && !need_step_diagnostics
    simulation_wall_t0 = time_ns()
    progress = Progress(options.cycles; desc="Running cycles: ", output=io, showspeed=true)
    for cycle in 1:options.cycles
        if use_cycle_step
            time_counted_block!(timings, :model_step_wall, ncol) do
                _model_cycle_step!(model, runtime, forcing)
            end
        else
        for t in eachindex(forcing.time_values)
            month_idx = need_monthly_outputs ? (cycle - 1) * schedule.nmonth_per_cycle + schedule.step_month[t] : 0
            time_counted_block!(timings, :model_step_wall, ncol) do
                _model_step!(model, runtime, forcing, t)
            end
            if need_step_diagnostics
                time_counted_block!(timings, :step_diagnostics, ncol) do
                    backend_step_summary === nothing && error("Missing backend summary buffers.")
                    _summarize_step_state!(backend_step_summary, model, runtime)
                    _copy_summary_fields!(step_summary, backend_step_summary, SUMMARY_BUFFER_NAMES)
                    _accumulate_step_diagnostics!(step_summary, previous, monthly_sums, step_vectors, month_idx, need_monthly_outputs, need_step_outputs)
                end
                need_monthly_outputs && (monthly_count[month_idx] += 1)
            end
            if need_step_outputs && schedule.annual_output.write_output[t]
                steps_written += 1
                if write_step_fields
                    step_grids = time_block!(timings, :step_output_prepare) do
                        _step_output_grids(step_vectors, grid)
                    end
                    time_block!(timings, :step_output_write) do
                        for key in OUTPUT_GROUPS.step
                            maybe_write_step_output!(writer, steps_written, key, getfield(step_grids, key))
                        end
                    end
                    _reset_step_vectors!(step_vectors)
                end
            end
        end
        end

        time_block!(timings, :summarize_columns_cycle) do
            _summarize_cycle_state!(backend_cycle_summary, model, runtime)
            _copy_summary_fields!(final, backend_cycle_summary, CYCLE_BUFFER_NAMES)
        end
        if should_record_cycle_metrics(cycle, options.cycles, options.history_stride)
            record = time_block!(timings, :cycle_metrics) do
                make_cycle_record_and_deltas!(cycle, deltas.thickness, deltas.wet_mass, deltas.base_mass, final.thickness, final.wet_mass, final.bulk_density, final.base_mass, prev.thickness, prev.wet_mass, prev.base_mass)
            end
            push!(history, record)
            time_block!(timings, :cycle_logging) do
                println(io, cycle_log_line(record))
            end
        end
        if need_last_cycle_smb_delta
            time_block!(timings, :cycle_state_deltas) do
                _update_cycle_smb_delta!(deltas.ice_sheet_smb, previous_cycle_smb_ice, model, runtime)
            end
        end
        prev, final = final, prev
        next!(progress)
    end
    simulation_wall_sec = (time_ns() - simulation_wall_t0) * 1.0e-9

    summary_path = ""
    history_csv_path = ""
    if options.write_outputs
        mkpath(options.output_dir)
        summary_path = joinpath(options.output_dir, "$(options.name)_summary.txt")
        history_csv_path = joinpath(options.output_dir, "$(options.name)_history.csv")
        time_block!(timings, :write_summary_text) do
            write_run_summary(summary_path, options, forcing.time_values, ncol, history, :cycles, timings)
        end
        time_block!(timings, :write_history_csv) do
            write_run_history_csv(history_csv_path, history)
        end
    end
    if writer !== nothing
        final_grids = time_block!(timings, :scatter_final_outputs) do
            _scatter_model_final_grids(model, runtime, prev, deltas, grid)
        end
        layer_grids = _model_layer_grids(model, runtime, grid, need_layer_outputs, timings)
        monthly_grids = need_monthly_outputs ? time_block!(timings, :aggregate_monthly_outputs) do
            _finalize_monthly_grids(monthly_sums, monthly_count, grid)
        end : empty_monthly_grids()
        time_block!(timings, :write_netcdf) do
            finalize_netcdf!(writer, final_grids, layer_grids, history, monthly_grids, :cycles, length(history), steps_written)
        end
    end

    _finalize_model_runtime!(model, runtime, options, timings)
    run_wall_sec = (time_ns() - run_wall_t0) * 1.0e-9
    print_run_report(io, options, forcing.time_values, history, :cycles, simulation_wall_sec, run_wall_sec, timings; nc_path=nc_path, summary_path=summary_path, history_csv_path=history_csv_path)
    return (
        history=history,
        status=:cycles,
        timings=timings,
        simulation_wall_sec=simulation_wall_sec,
        run_wall_sec=run_wall_sec,
        netcdf_path=nc_path,
        summary_path=summary_path,
        history_csv_path=history_csv_path,
    )
end

const SUMMARY_BUFFER_NAMES = (:thickness, :wet_mass, :bulk_density, :base_mass, :smb_ice, :liquid_water, :runoff, :pdd)
const CYCLE_BUFFER_NAMES = (:thickness, :wet_mass, :bulk_density, :base_mass)

_named_buffers(names::NTuple{N, Symbol}, build::F) where {N, F <: Function} = NamedTuple{names}(ntuple(_ -> build(), N))

allocate_summary_buffers(n::Int) = _named_buffers(SUMMARY_BUFFER_NAMES, () -> Vector{Float64}(undef, n))
allocate_summary_buffers(domain::AbstractSnowpackDomain, n::Int) = _named_buffers(SUMMARY_BUFFER_NAMES, () -> similar(domain.mass, Float64, n))
allocate_cycle_summary_buffers(n::Int) = _named_buffers(CYCLE_BUFFER_NAMES, () -> Vector{Float64}(undef, n))
allocate_cycle_summary_buffers(domain::AbstractSnowpackDomain, n::Int) = _named_buffers(CYCLE_BUFFER_NAMES, () -> similar(domain.mass, Float64, n))

"""
Simulation-first public API for Chion.
"""

"""
    SimulationOptions(; name="chion_run", input_label="", cycles=1, backend=:threads, history_stride=1)

Execution controls for a `Simulation`.
"""
struct SimulationOptions
    name::String
    input_label::String
    cycles::Int
    backend::Symbol
    history_stride::Int
end

function SimulationOptions(;
    name::AbstractString="chion_run",
    input_label::AbstractString="",
    cycles::Integer=1,
    backend=:threads,
    history_stride::Integer=1,
)
    return SimulationOptions(
        String(name),
        String(input_label),
        Int(cycles),
        normalize_backend(backend),
        normalize_history_stride(history_stride),
    )
end

"""
    OutputOptions(; save=Symbol[], output_dir="", netcdf_path="", write_outputs=false)

Output controls for text, CSV, and NetCDF products written by `run!`.
"""
struct OutputOptions
    variables::Vector{Symbol}
    output_dir::String
    netcdf_path::String
    write_outputs::Bool
end

function OutputOptions(;
    save=Symbol[],
    output_dir::AbstractString="",
    netcdf_path::AbstractString="",
    write_outputs::Bool=false,
)
    return OutputOptions(
        _normalize_run_save(save),
        String(output_dir),
        String(netcdf_path),
        write_outputs,
    )
end

"""
    SimulationResult

Summary returned by `run!`, including cycle history, timing diagnostics, and
paths to any files written during the run.
"""
struct SimulationResult
    history::Vector{NamedTuple}
    status::Symbol
    timings::StepTimingStats
    simulation_wall_sec::Float64
    run_wall_sec::Float64
    netcdf_path::String
    summary_path::String
    history_csv_path::String
end

function _copy_domain_state!(dest::SnowpackDomain, src::SnowpackDomain)
    dest.Ntot == src.Ntot || error("Cannot copy domain state with different `Ntot`.")
    dest.ncol == src.ncol || error("Cannot copy domain state with different column count.")
    dest.N .= src.N
    dest.mass .= src.mass
    dest.mass_w .= src.mass_w
    dest.density .= src.density
    dest.temperature .= src.temperature
    dest.mass_base .= src.mass_base
    dest.smb_ice .= src.smb_ice
    dest.runoff .= src.runoff
    dest.Tsrf .= src.Tsrf
    dest.snow_cover .= src.snow_cover
    dest.albedo_dynamic .= src.albedo_dynamic
    return dest
end

function _run_options(options::SimulationOptions, output::OutputOptions)
    return RunOptions(
        name=options.name,
        input_label=options.input_label,
        cycles=options.cycles,
        backend=options.backend,
        history_stride=options.history_stride,
        write_outputs=output.write_outputs,
        output_dir=output.output_dir,
        netcdf_path=output.netcdf_path,
        write_netcdf=!isempty(output.variables),
        netcdf_variables=output.variables,
    )
end

function _simulation_result(result)
    return SimulationResult(
        result.history,
        result.status,
        result.timings,
        result.simulation_wall_sec,
        result.run_wall_sec,
        result.netcdf_path,
        result.summary_path,
        result.history_csv_path,
    )
end

"""
    Simulation(model; forcing, options=SimulationOptions(), output=OutputOptions(), ...)

Couple a snow model with time-varying forcing and execution/output options.
"""
mutable struct Simulation{M <: AbstractSnowModel}
    model::M
    forcing::SnowpackForcing
    Δt::Float64
    stop_time::Union{Nothing, Float64}
    options::SimulationOptions
    output::OutputOptions
end

function Simulation(
    model::AbstractSnowModel;
    forcing::SnowpackForcing,
    Δt::Real=1.0,
    stop_time=nothing,
    options::SimulationOptions=SimulationOptions(),
    output::OutputOptions=OutputOptions(),
    cycles::Union{Nothing, Integer}=nothing,
    backend=nothing,
    history_stride::Union{Nothing, Integer}=nothing,
    save=nothing,
    output_dir::Union{Nothing, AbstractString}=nothing,
    netcdf_path::Union{Nothing, AbstractString}=nothing,
    write_outputs::Union{Nothing, Bool}=nothing,
    name::Union{Nothing, AbstractString}=nothing,
    input_label::Union{Nothing, AbstractString}=nothing,
)
    resolved_options = options
    if !isnothing(cycles) || !isnothing(backend) || !isnothing(history_stride) || !isnothing(name) || !isnothing(input_label)
        resolved_options = SimulationOptions(
            name=isnothing(name) ? options.name : name,
            input_label=isnothing(input_label) ? options.input_label : input_label,
            cycles=isnothing(cycles) ? options.cycles : cycles,
            backend=isnothing(backend) ? options.backend : backend,
            history_stride=isnothing(history_stride) ? options.history_stride : history_stride,
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
    return Simulation(
        model,
        forcing,
        Float64(Δt),
        isnothing(stop_time) ? nothing : Float64(stop_time),
        resolved_options,
        resolved_output,
    )
end

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
    println(io, "  columns: ", ncols(sim.model.grid))
    println(io, "  forcing steps: ", length(sim.forcing.time_values))
    println(io, "  backend: ", sim.options.backend)
    println(io, "  cycles: ", sim.options.cycles)
    println(io, "  netcdf vars: ", isempty(sim.output.variables) ? "(none)" : join(string.(sim.output.variables), ", "))
end

"""
    run!(simulation; options=simulation.options, output=simulation.output, io=stdout)

Advance a `Simulation` in place and return a `SimulationResult`.
"""
function run!(
    sim::Simulation;
    options::SimulationOptions=sim.options,
    output::OutputOptions=sim.output,
    io::IO=stdout,
)
    run_options = _run_options(options, output)
    result = execute_model_run!(sim.model, sim.forcing; options=run_options, io=io)
    return _simulation_result(result)
end
