function time_block!(stats, key::Symbol, f; synchronize=nothing)
    synchronize === nothing || synchronize()
    t0 = time_ns()
    value = f()
    synchronize === nothing || synchronize()
    add_timing!(stats, key, (time_ns() - t0) * 1.0e-9)
    return value
end

time_block!(f, stats, key::Symbol; kwargs...) = time_block!(stats, key, f; kwargs...)

function time_counted_block!(stats, key::Symbol, count::Int, f; synchronize=nothing)
    synchronize === nothing || synchronize()
    t0 = time_ns()
    value = f()
    synchronize === nothing || synchronize()
    add_timing!(stats, key, (time_ns() - t0) * 1.0e-9, count)
    return value
end

time_counted_block!(f, stats, key::Symbol, count::Int; kwargs...) =
    time_counted_block!(stats, key, count, f; kwargs...)

using Dates
using Base.Threads: @threads
using NCDatasets
import CUDA
import Libdl

include("runtime_core.jl")
include("netcdf.jl")
include("reporting.jl")

function _prepare_backend!(timings::StepTimingStats, domain::SnowpackDomain, forcing::ForcingData; is_gpu::Bool)
    step_fields = SnowpackStepFields(forcing)
    if is_gpu
        cuda_available() || error("`backend=gpu` requested, but CUDA is not functional in the current environment.")
        domain = time_block!(timings, :gpu_transfer) do
            gpu_domain(domain)
        end
        step_fields = time_block!(timings, :gpu_transfer) do
            adapt(CUDA.CuArray, step_fields)
        end
        workspaces = time_block!(timings, :gpu_transfer) do
            ColumnarStepWorkspace(domain)
        end
        return domain, step_fields, workspaces
    end
    workspaces = time_block!(timings, :create_workspaces) do
        threaded_workspaces(domain)
    end
    return domain, step_fields, workspaces
end

function execute_run!(
    domain::SnowpackDomain,
    forcing::ForcingData;
    layout::Union{Nothing, GridLayout}=nothing,
    options::RunOptions=RunOptions(),
    io::IO=stdout,
    timings::StepTimingStats=StepTimingStats(),
    run_wall_t0::Integer=time_ns(),
)
    ncol = column_count(domain)
    size(forcing.air_temperature, 1) == ncol || error("Forcing column count must match the domain column count.")
    options.write_netcdf && isnothing(layout) && error("NetCDF output requires a grid layout.")
    isnothing(layout) || length(layout.js) == ncol || error("Grid-layout point count must match the domain column count.")

    selected = Set(options.netcdf_variables)
    is_gpu = options.backend == :gpu
    summary_backend = is_gpu ? :kernelabstractions : :threads
    synchronize = is_gpu ? CUDA.synchronize : nothing
    write_final_fields = options.write_netcdf && !isnothing(layout)
    need_step_outputs = options.write_netcdf && any(var -> var in selected, OUTPUT_GROUPS.step)
    need_monthly_outputs = options.write_netcdf && any(var -> var in selected, OUTPUT_GROUPS.monthly)
    need_layer_outputs = options.write_netcdf && any(var -> var in selected, OUTPUT_GROUPS.layers)
    need_last_cycle_smb_delta = options.write_netcdf && (:last_cycle_delta_ice_sheet_smb in selected)
    need_step_diagnostics = need_step_outputs || need_monthly_outputs

    initial_thickness_vec = if write_final_fields
        summary = allocate_cycle_summary_buffers(ncol)
        time_block!(timings, :prepare_initial_output_fields) do
            summarize_cycle_state!(summary.thickness, summary.wet_mass, summary.bulk_density, summary.base_mass, domain)
        end
        copy(summary.thickness)
    else
        Float64[]
    end

    domain, step_fields, workspaces = _prepare_backend!(timings, domain, forcing; is_gpu=is_gpu)
    schedule = nothing
    writer = nothing
    nc_path = ""
    if options.write_netcdf
        schedule = time_block!(timings, :prepare_output_schedule) do
            _prepare_output_schedule(forcing.time_values, options.cycles)
        end
        initial_thickness = time_block!(timings, :prepare_initial_output_fields_grid) do
            scatter_to_grid(initial_thickness_vec, layout.js, layout.is, _grid_shape(layout))
        end
        nc_path = resolve_netcdf_path(options)
        writer = time_block!(timings, :init_netcdf) do
            init_netcdf(
                nc_path,
                options,
                forcing.time_values,
                domain.Ntot,
                layout,
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

    prev = allocate_cycle_summary_buffers(ncol)
    final = allocate_cycle_summary_buffers(ncol)
    device_cycle_summary = is_gpu ? allocate_cycle_summary_buffers(domain, ncol) : nothing
    step_summary = need_step_diagnostics ? allocate_summary_buffers(ncol) : nothing
    device_step_summary = is_gpu && need_step_diagnostics ? allocate_summary_buffers(domain, ncol) : nothing
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

    time_block!(timings, :summarize_columns_initial; synchronize=synchronize) do
        if summary_backend == :kernelabstractions
            device_cycle_summary === nothing && error("`device_summary` must be provided for `backend=:kernelabstractions`.")
            summarize_cycle_state!(
                device_cycle_summary.thickness,
                device_cycle_summary.wet_mass,
                device_cycle_summary.bulk_density,
                device_cycle_summary.base_mass,
                domain;
                backend=summary_backend,
            )
            _copy_summary_fields!(prev, device_cycle_summary, CYCLE_BUFFER_NAMES)
        else
            summarize_cycle_state!(
                prev.thickness,
                prev.wet_mass,
                prev.bulk_density,
                prev.base_mass,
                domain;
                backend=summary_backend,
            )
        end
    end
    previous = (
        base_mass=_host_vector(prev.base_mass; copy_array=true),
        smb_ice=_host_vector(domain.smb_ice; copy_array=true),
        runoff=_host_vector(domain.runoff; copy_array=true),
    )
    previous_cycle_smb_ice = need_last_cycle_smb_delta ? _host_vector(domain.smb_ice; copy_array=true) : Float64[]
    history = NamedTuple[]
    steps_written = 0

    simulation_wall_t0 = time_ns()
    for cycle in 1:options.cycles
        for t in eachindex(forcing.time_values)
            month_idx = need_monthly_outputs ? (cycle - 1) * schedule.nmonth_per_cycle + schedule.step_month[t] : 0
            time_counted_block!(timings, :model_step_wall, ncol; synchronize=synchronize) do
                step!(domain, step_fields, t, workspaces)
            end
            if need_step_diagnostics
                time_counted_block!(timings, :step_diagnostics, ncol; synchronize=synchronize) do
                    if summary_backend == :kernelabstractions
                        device_step_summary === nothing && error("`device_summary` must be provided for `backend=:kernelabstractions`.")
                        summarize_domain_state!(
                            device_step_summary.thickness,
                            device_step_summary.wet_mass,
                            device_step_summary.bulk_density,
                            device_step_summary.base_mass,
                            device_step_summary.smb_ice,
                            device_step_summary.liquid_water,
                            device_step_summary.runoff,
                            domain;
                            backend=summary_backend,
                        )
                        _copy_summary_fields!(step_summary, device_step_summary, SUMMARY_BUFFER_NAMES)
                    else
                        summarize_domain_state!(
                            step_summary.thickness,
                            step_summary.wet_mass,
                            step_summary.bulk_density,
                            step_summary.base_mass,
                            step_summary.smb_ice,
                            step_summary.liquid_water,
                            step_summary.runoff,
                            domain;
                            backend=summary_backend,
                        )
                    end
                    _accumulate_step_diagnostics!(step_summary, previous, monthly_sums, step_vectors, month_idx, need_monthly_outputs, need_step_outputs)
                end
                need_monthly_outputs && (monthly_count[month_idx] += 1)
            end
        if need_step_outputs && schedule.annual_output.write_output[t]
            steps_written += 1
            if write_step_fields
                step_grids = time_block!(timings, :step_output_prepare) do
                    _step_output_grids(step_vectors, layout)
                end
                    time_block!(timings, :step_output_write) do
                        for key in OUTPUT_GROUPS.step
                            maybe_write_step_output!(writer, steps_written, key, getfield(step_grids, key))
                        end
                    end
                end
                _reset_step_vectors!(step_vectors)
            end
        end

        time_block!(timings, :summarize_columns_cycle; synchronize=synchronize) do
            if summary_backend == :kernelabstractions
                device_cycle_summary === nothing && error("`device_summary` must be provided for `backend=:kernelabstractions`.")
                summarize_cycle_state!(
                    device_cycle_summary.thickness,
                    device_cycle_summary.wet_mass,
                    device_cycle_summary.bulk_density,
                    device_cycle_summary.base_mass,
                    domain;
                    backend=summary_backend,
                )
                _copy_summary_fields!(final, device_cycle_summary, CYCLE_BUFFER_NAMES)
            else
                summarize_cycle_state!(
                    final.thickness,
                    final.wet_mass,
                    final.bulk_density,
                    final.base_mass,
                    domain;
                    backend=summary_backend,
                )
            end
        end
        if should_record_cycle_metrics(cycle, options.cycles, options.history_stride)
            record = time_block!(timings, :cycle_metrics; synchronize=synchronize) do
                make_cycle_record_and_deltas!(cycle, deltas.thickness, deltas.wet_mass, deltas.base_mass, final.thickness, final.wet_mass, final.bulk_density, final.base_mass, prev.thickness, prev.wet_mass, prev.base_mass)
            end
            push!(history, record)
            time_block!(timings, :cycle_logging) do
                println(io, cycle_log_line(record))
            end
        end
        if need_last_cycle_smb_delta
            time_block!(timings, :cycle_state_deltas; synchronize=synchronize) do
                _update_cycle_smb_delta!(deltas.ice_sheet_smb, previous_cycle_smb_ice, domain)
            end
        end
        prev, final = final, prev
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
            _scatter_final_grids(prev, domain, deltas, layout)
        end
        layer_grids = if need_layer_outputs
            final_domain = is_gpu ? time_block!(timings, :gpu_transfer) do
                cpu_domain(domain)
            end : domain
            time_block!(timings, :collect_final_layer_grids) do
                collect_final_layer_grids(final_domain, layout.js, layout.is, _grid_shape(layout), final_domain.Ntot)
            end
        else
            _empty_layer_grids()
        end
        monthly_grids = need_monthly_outputs ? time_block!(timings, :aggregate_monthly_outputs) do
            _finalize_monthly_grids(monthly_sums, monthly_count, layout)
        end : empty_monthly_grids()
        time_block!(timings, :write_netcdf) do
            finalize_netcdf!(writer, final_grids, layer_grids, history, monthly_grids, :cycles, length(history), steps_written)
        end
    end

    run_wall_sec = (time_ns() - run_wall_t0) * 1.0e-9
    print_run_report(io, options, forcing.time_values, history, :cycles, simulation_wall_sec, run_wall_sec, timings; nc_path=nc_path, summary_path=summary_path, history_csv_path=history_csv_path)
    return RunResult(history, :cycles, timings, simulation_wall_sec, run_wall_sec, nc_path, summary_path, history_csv_path, domain, options)
end
const SUMMARY_BUFFER_NAMES = (:thickness, :wet_mass, :bulk_density, :base_mass, :smb_ice, :liquid_water, :runoff)
const CYCLE_BUFFER_NAMES = (:thickness, :wet_mass, :bulk_density, :base_mass)

_named_buffers(names::NTuple{N, Symbol}, build::F) where {N, F <: Function} = NamedTuple{names}(ntuple(_ -> build(), N))

allocate_summary_buffers(n::Int) = _named_buffers(SUMMARY_BUFFER_NAMES, () -> Vector{Float64}(undef, n))
allocate_summary_buffers(domain::AbstractSnowpackDomain, n::Int) = _named_buffers(SUMMARY_BUFFER_NAMES, () -> similar(domain.mass, Float64, n))
allocate_cycle_summary_buffers(n::Int) = _named_buffers(CYCLE_BUFFER_NAMES, () -> Vector{Float64}(undef, n))
allocate_cycle_summary_buffers(domain::AbstractSnowpackDomain, n::Int) = _named_buffers(CYCLE_BUFFER_NAMES, () -> similar(domain.mass, Float64, n))

@inline _host_vector(data::Vector{Float64}; copy_array::Bool=false) = copy_array ? copy(data) : data
@inline _host_vector(data; copy_array::Bool=false) = Float64.(Array(data))

function _prepare_output_schedule(time_values::Vector{DateTime}, cycles::Int)
    write_output = falses(length(time_values))
    output_slot = zeros(Int, length(time_values))
    source_indices = Int32[]
    source_codes = Int32[]
    years = unique(year.(time_values))
    for (slot, yr) in enumerate(years)
        last_t = findlast(t -> year(time_values[t]) == yr, eachindex(time_values))
        isnothing(last_t) && error("Could not determine the last timestep for source year $yr.")
        write_output[last_t] = true
        output_slot[last_t] = slot
        push!(source_indices, Int32(last_t))
        ts = time_values[last_t]
        push!(source_codes, Int32(year(ts) * 1000000 + month(ts) * 10000 + day(ts) * 100 + hour(ts)))
    end
    month_keys = unique((year(t), month(t)) for t in time_values)
    month_lookup = Dict(key => idx for (idx, key) in enumerate(month_keys))
    return (
        annual_output=(write_output=write_output, output_slot=output_slot, source_indices=source_indices, source_codes=source_codes, years=years),
        step_month=[month_lookup[(year(t), month(t))] for t in time_values],
        nmonth_per_cycle=length(month_keys),
        nmonth_total=cycles * length(month_keys),
        month_cycle=Int32[cyc for cyc in 1:cycles for _ in month_keys],
        month_of_year=Int32[key[2] for _ in 1:cycles for key in month_keys],
        source_month_code=Int32[key[1] * 100 + key[2] for _ in 1:cycles for key in month_keys],
    )
end
