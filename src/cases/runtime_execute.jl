@inline _backend_info(options::RunConfig) = (
    is_gpu=options.backend == :gpu,
    summary_backend=options.backend == :gpu ? :kernelabstractions : :threads,
    sync=options.backend == :gpu ? CUDA.synchronize : nothing,
)

@inline _step_cycle!(domain, step_fields, workspaces, options::RunConfig) =
    options.backend == :gpu ? SM.step!(domain, step_fields, workspaces; update_snow_cover=false) : SM.step!(domain, step_fields, workspaces)

@inline _step_timestep!(domain, step_fields, t::Int, workspaces, ::RunConfig) = SM.step!(domain, step_fields, t, workspaces)

function _prepare_backend!(timings::TimingStats, options::RunConfig, domain::SM.SnowpackDomain, forcing::ForcingData)
    step_fields = SM.SnowpackStepFields(forcing)
    if options.backend == :gpu
        SM.cuda_available() || error("`backend=gpu` requested, but CUDA is not functional in the current environment.")
        domain = time_block!(timings, :gpu_transfer) do
            SM.gpu_domain(domain)
        end
        step_fields = time_block!(timings, :gpu_transfer) do
            SM.adapt(CUDA.CuArray, step_fields)
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

function _case_flags(options::RunConfig, layout::Union{Nothing, GridLayout})
    selected = Set(options.netcdf_variables)
    need_step_outputs = options.write_netcdf && case_selected(options, :step)
    need_monthly_outputs = options.write_netcdf && case_selected(options, :monthly)
    return (
        write_final_fields=options.write_netcdf && !isnothing(layout),
        need_step_outputs=need_step_outputs,
        need_monthly_outputs=need_monthly_outputs,
        need_layer_outputs=options.write_netcdf && case_selected(options, :layers),
        need_last_cycle_smb_delta=options.write_netcdf && (:last_cycle_delta_ice_sheet_smb in selected),
        need_step_diagnostics=need_step_outputs || need_monthly_outputs,
    )
end

function _prepare_output_schedule(time_values::Vector{DateTime}, cycles::Int)
    annual_output = build_annual_output_schedule(time_values)
    month_keys = unique((year(t), month(t)) for t in time_values)
    month_lookup = Dict(key => idx for (idx, key) in enumerate(month_keys))
    return (
        annual_output=annual_output,
        step_month=[month_lookup[(year(t), month(t))] for t in time_values],
        nmonth_per_cycle=length(month_keys),
        nmonth_total=cycles * length(month_keys),
        month_cycle=Int32[cyc for cyc in 1:cycles for _ in month_keys],
        month_of_year=Int32[key[2] for _ in 1:cycles for key in month_keys],
        source_month_code=Int32[key[1] * 100 + key[2] for _ in 1:cycles for key in month_keys],
    )
end

_allocate_step_vectors(active::Bool, ncol::Int) = NamedTuple{CASE_OUTPUT_GROUPS.step}(ntuple(_ -> active ? zeros(Float64, ncol) : Float64[], length(CASE_OUTPUT_GROUPS.step)))
_allocate_monthly_sums(active::Bool, nmonth_total::Int, ncol::Int) = NamedTuple{MONTHLY_GRID_KEYS}(ntuple(_ -> active ? zeros(Float64, nmonth_total, ncol) : Matrix{Float64}(undef, 0, 0), length(MONTHLY_GRID_KEYS)))

function _step_output_grids(step_vectors, layout::GridLayout)
    return NamedTuple{CASE_OUTPUT_GROUPS.step}(ntuple(i -> scatter_to_grid(getfield(step_vectors, CASE_OUTPUT_GROUPS.step[i]), layout.js, layout.is, _grid_shape(layout)), length(CASE_OUTPUT_GROUPS.step)))
end

function _reset_step_vectors!(step_vectors)
    for key in CASE_OUTPUT_GROUPS.step
        isempty(getfield(step_vectors, key)) || fill!(getfield(step_vectors, key), 0.0)
    end
    return
end

function _accumulate_step_diagnostics!(summary, previous, monthly_sums, step_vectors, month_idx::Int, flags)
    current_base = summary.base_mass
    current_smb = summary.smb_ice
    current_runoff = summary.runoff
    delta_base = current_base .- previous.base_mass
    delta_smb = current_smb .- previous.smb_ice
    delta_runoff = current_runoff .- previous.runoff

    if flags.need_monthly_outputs
        monthly_sums.monthly_mean_thickness[month_idx, :] .+= summary.thickness
        monthly_sums.monthly_mean_wet_mass[month_idx, :] .+= summary.wet_mass
        monthly_sums.monthly_mean_bulk_density[month_idx, :] .+= summary.bulk_density
        monthly_sums.monthly_mean_base_mass[month_idx, :] .+= current_base
        monthly_sums.monthly_mean_ice_sheet_smb[month_idx, :] .+= delta_smb
        monthly_sums.monthly_export_to_ice[month_idx, :] .+= delta_base
        monthly_sums.monthly_net_ice_sheet_forcing[month_idx, :] .+= delta_smb
        monthly_sums.monthly_runoff[month_idx, :] .+= delta_runoff
    end
    if flags.need_step_outputs
        step_vectors.step_export_to_ice .+= delta_base
        step_vectors.step_ice_sheet_smb .+= delta_smb
    end
    previous.base_mass .= current_base
    previous.smb_ice .= current_smb
    previous.runoff .= current_runoff
    return
end

function _update_cycle_smb_delta!(last_delta::Vector{Float64}, previous_cycle_smb_ice::Vector{Float64}, domain)
    current = _host_vector(domain.smb_ice; copy_array=true)
    last_delta .= current .- previous_cycle_smb_ice
    previous_cycle_smb_ice .= current
    return
end

function _scatter_final_grids(final_state, domain, deltas, layout::GridLayout)
    final_smb_ice = _host_vector(domain.smb_ice; copy_array=true)
    final_runoff = _host_vector(domain.runoff; copy_array=true)
    return (
        final_thickness=scatter_to_grid(final_state.thickness, layout.js, layout.is, _grid_shape(layout)),
        final_wet_mass=scatter_to_grid(final_state.wet_mass, layout.js, layout.is, _grid_shape(layout)),
        final_bulk_density=scatter_to_grid(final_state.bulk_density, layout.js, layout.is, _grid_shape(layout)),
        final_base_mass=scatter_to_grid(final_state.base_mass, layout.js, layout.is, _grid_shape(layout)),
        final_ice_sheet_smb=scatter_to_grid(final_smb_ice, layout.js, layout.is, _grid_shape(layout)),
        final_runoff=scatter_to_grid(final_runoff, layout.js, layout.is, _grid_shape(layout)),
        last_cycle_delta_thickness=scatter_to_grid(deltas.thickness, layout.js, layout.is, _grid_shape(layout)),
        last_cycle_delta_wet_mass=scatter_to_grid(deltas.wet_mass, layout.js, layout.is, _grid_shape(layout)),
        last_cycle_delta_base_mass=scatter_to_grid(deltas.base_mass, layout.js, layout.is, _grid_shape(layout)),
        last_cycle_delta_ice_sheet_smb=scatter_to_grid(deltas.ice_sheet_smb, layout.js, layout.is, _grid_shape(layout)),
    )
end

const MONTHLY_MEAN_KEYS = (:monthly_mean_thickness, :monthly_mean_wet_mass, :monthly_mean_bulk_density)

function _finalize_monthly_grids(monthly_sums, monthly_count::Vector{Int32}, layout::GridLayout)
    vectors = NamedTuple{MONTHLY_GRID_KEYS}(ntuple(i -> begin
        key = MONTHLY_GRID_KEYS[i]
        data = copy(getfield(monthly_sums, key))
        if key in MONTHLY_MEAN_KEYS
            @inbounds for m in axes(data, 1)
                data[m, :] ./= max(monthly_count[m], 1)
            end
        end
        data
    end, length(MONTHLY_GRID_KEYS)))
    return NamedTuple{MONTHLY_GRID_KEYS}(ntuple(i -> monthly_vectors_to_grids(getfield(vectors, MONTHLY_GRID_KEYS[i]), layout.js, layout.is, _grid_shape(layout)), length(MONTHLY_GRID_KEYS)))
end

function _write_optional_outputs!(
    timings::TimingStats,
    options::RunConfig,
    time_values::Vector{DateTime},
    ncol::Int,
    history::Vector{NamedTuple},
    status::Symbol,
)
    options.write_outputs || return "", ""
    mkpath(options.output_dir)
    summary_path = joinpath(options.output_dir, "$(options.name)_summary.txt")
    history_csv_path = joinpath(options.output_dir, "$(options.name)_history.csv")
    time_block!(timings, :write_summary_text) do
        write_case_summary(summary_path, options, time_values, ncol, history, status, timings)
    end
    time_block!(timings, :write_history_csv) do
        write_case_history_csv(history_csv_path, history)
    end
    return summary_path, history_csv_path
end

function run_case_cycles_no_netcdf!(
    timings::TimingStats,
    options::RunConfig,
    domain::SM.SnowpackDomain,
    workspaces,
    step_fields::SM.SnowpackStepFields,
    io::IO,
)
    backend = _backend_info(options)
    ncol = SM.column_count(domain)
    prev = allocate_cycle_summary_buffers(ncol)
    final = allocate_cycle_summary_buffers(ncol)
    device_cycle_summary = backend.is_gpu ? allocate_cycle_summary_buffers(domain, ncol) : nothing
    deltas = (
        thickness=fill(NaN, ncol),
        wet_mass=fill(NaN, ncol),
        base_mass=fill(NaN, ncol),
    )
    time_block!(timings, :summarize_columns_initial; synchronize=backend.sync) do
        summarize_cycle_columns!(prev, domain; backend=backend.summary_backend, device_summary=device_cycle_summary)
    end

    history = NamedTuple[]
    simulation_wall_t0 = time_ns()
    for cycle in 1:options.cycles
        time_counted_block!(timings, :model_step_wall, ncol * size(step_fields.air_temperature, 2); synchronize=backend.sync) do
            _step_cycle!(domain, step_fields, workspaces, options)
        end
        time_block!(timings, :summarize_columns_cycle; synchronize=backend.sync) do
            summarize_cycle_columns!(final, domain; backend=backend.summary_backend, device_summary=device_cycle_summary)
        end
        if should_record_cycle_metrics(cycle, options.cycles, options.history_stride)
            record = time_block!(timings, :cycle_metrics; synchronize=backend.sync) do
                make_cycle_record_and_deltas!(cycle, deltas.thickness, deltas.wet_mass, deltas.base_mass, final.thickness, final.wet_mass, final.bulk_density, final.base_mass, prev.thickness, prev.wet_mass, prev.base_mass)
            end
            push!(history, record)
            time_block!(timings, :cycle_logging) do
                println(io, cycle_log_line(record))
            end
        end
        prev, final = final, prev
    end
    return (history=history, status=:cycles, simulation_wall_sec=(time_ns() - simulation_wall_t0) * 1.0e-9)
end

function _print_run_report(
    io::IO,
    options::RunConfig,
    time_values::Vector{DateTime},
    history::Vector{NamedTuple},
    status::Symbol,
    simulation_wall_sec::Float64,
    run_wall_sec::Float64,
    timings::TimingStats;
    nc_path::AbstractString="",
    summary_path::AbstractString="",
    history_csv_path::AbstractString="",
)
    println(io, "$(options.name) complete.")
    println(io, "Input label     : ", isempty(options.input_label) ? "(not provided)" : options.input_label)
    println(io, "Forcing start   : $(first(time_values))")
    println(io, "Forcing end     : $(last(time_values))")
    println(io, "Backend         : $(String(options.backend))")
    println(io, "Cycles          : $(completed_cycle_count(history, status, options.cycles))")
    println(io, "Status          : $(string(status))")
    println(io, "Cycle metrics   : $(cycle_metrics_schedule_label(options.history_stride))")
    println(io, @sprintf("Simulation wall : %.3f s", simulation_wall_sec))
    println(io, @sprintf("Run wall total  : %.3f s", run_wall_sec))
    haskey(timings.totals, :model_step_wall) && println(io, @sprintf("Model step wall : %.3f s", timings.totals[:model_step_wall]))
    println(io, "Output NetCDF   : ", options.write_netcdf ? abspath(nc_path) : "skipped (--no-nc)")
    if options.write_outputs
        println(io, "History CSV     : $(abspath(history_csv_path))")
        println(io, "Summary         : $(abspath(summary_path))")
    else
        println(io, "File outputs    : skipped (--no-output)")
    end
    print_timing_summary(io, timings; total_wall_sec=run_wall_sec)
end

function execute_case!(
    domain::SM.SnowpackDomain,
    forcing::ForcingData;
    layout::Union{Nothing, GridLayout}=nothing,
    options::RunConfig=RunConfig(),
    io::IO=stdout,
    timings::TimingStats=TimingStats(),
    run_wall_t0::Integer=time_ns(),
)
    ncol = SM.column_count(domain)
    size(forcing.air_temperature, 1) == ncol || error("Forcing column count must match the domain column count.")
    options.write_netcdf && isnothing(layout) && error("NetCDF output requires a grid layout.")
    isnothing(layout) || length(layout.js) == ncol || error("Grid-layout point count must match the domain column count.")

    flags = _case_flags(options, layout)
    initial_thickness_vec = if flags.write_final_fields
        summary = allocate_cycle_summary_buffers(ncol)
        time_block!(timings, :prepare_initial_output_fields) do
            summarize_cycle_columns!(summary, domain)
        end
        copy(summary.thickness)
    else
        Float64[]
    end

    domain, step_fields, workspaces = _prepare_backend!(timings, options, domain, forcing)
    if !options.write_netcdf
        sim = run_case_cycles_no_netcdf!(timings, options, domain, workspaces, step_fields, io)
        summary_path, history_csv_path = _write_optional_outputs!(timings, options, forcing.time_values, ncol, sim.history, sim.status)
        run_wall_sec = (time_ns() - run_wall_t0) * 1.0e-9
        _print_run_report(io, options, forcing.time_values, sim.history, sim.status, sim.simulation_wall_sec, run_wall_sec, timings; summary_path=summary_path, history_csv_path=history_csv_path)
        return RunResult(sim.history, sim.status, timings, sim.simulation_wall_sec, run_wall_sec, "", summary_path, history_csv_path, domain, options)
    end

    backend = _backend_info(options)
    schedule = time_block!(timings, :prepare_output_schedule) do
        _prepare_output_schedule(forcing.time_values, options.cycles)
    end
    initial_thickness = time_block!(timings, :prepare_initial_output_fields_grid) do
        scatter_to_grid(initial_thickness_vec, layout.js, layout.is, _grid_shape(layout))
    end
    nc_path = resolve_case_netcdf_path(options)
    writer = time_block!(timings, :init_netcdf) do
        init_case_netcdf(
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
    write_step_fields = any(key -> haskey(writer.vars, key), CASE_OUTPUT_GROUPS.step)

    prev = allocate_cycle_summary_buffers(ncol)
    final = allocate_cycle_summary_buffers(ncol)
    device_cycle_summary = backend.is_gpu ? allocate_cycle_summary_buffers(domain, ncol) : nothing
    step_summary = flags.need_step_diagnostics ? allocate_summary_buffers(ncol) : nothing
    device_step_summary = backend.is_gpu && flags.need_step_diagnostics ? allocate_summary_buffers(domain, ncol) : nothing
    deltas = (
        thickness=fill(NaN, ncol),
        wet_mass=fill(NaN, ncol),
        base_mass=fill(NaN, ncol),
        ice_sheet_smb=flags.need_last_cycle_smb_delta ? fill(NaN, ncol) : Float64[],
    )
    monthly_sums = _allocate_monthly_sums(flags.need_monthly_outputs, schedule.nmonth_total, ncol)
    monthly_count = flags.need_monthly_outputs ? zeros(Int32, schedule.nmonth_total) : Int32[]
    step_vectors = _allocate_step_vectors(flags.need_step_outputs, ncol)

    time_block!(timings, :summarize_columns_initial; synchronize=backend.sync) do
        summarize_cycle_columns!(prev, domain; backend=backend.summary_backend, device_summary=device_cycle_summary)
    end
    previous = (
        base_mass=_host_vector(prev.base_mass; copy_array=true),
        smb_ice=_host_vector(domain.smb_ice; copy_array=true),
        runoff=_host_vector(domain.runoff; copy_array=true),
    )
    previous_cycle_smb_ice = flags.need_last_cycle_smb_delta ? _host_vector(domain.smb_ice; copy_array=true) : Float64[]
    history = NamedTuple[]
    steps_written = 0

    simulation_wall_t0 = time_ns()
    for cycle in 1:options.cycles
        for t in eachindex(forcing.time_values)
            month_idx = flags.need_monthly_outputs ? (cycle - 1) * schedule.nmonth_per_cycle + schedule.step_month[t] : 0
            time_counted_block!(timings, :model_step_wall, ncol; synchronize=backend.sync) do
                _step_timestep!(domain, step_fields, t, workspaces, options)
            end
            if flags.need_step_diagnostics
                time_counted_block!(timings, :step_diagnostics, ncol; synchronize=backend.sync) do
                    summarize_columns!(step_summary, domain; backend=backend.summary_backend, device_summary=device_step_summary)
                    _accumulate_step_diagnostics!(step_summary, previous, monthly_sums, step_vectors, month_idx, flags)
                end
                flags.need_monthly_outputs && (monthly_count[month_idx] += 1)
            end
            if flags.need_step_outputs && schedule.annual_output.write_output[t]
                steps_written += 1
                if write_step_fields
                    step_grids = time_block!(timings, :step_output_prepare) do
                        _step_output_grids(step_vectors, layout)
                    end
                    time_block!(timings, :step_output_write) do
                        for key in CASE_OUTPUT_GROUPS.step
                            maybe_write_step_output!(writer, steps_written, key, getfield(step_grids, key))
                        end
                    end
                end
                _reset_step_vectors!(step_vectors)
            end
        end

        time_block!(timings, :summarize_columns_cycle; synchronize=backend.sync) do
            summarize_cycle_columns!(final, domain; backend=backend.summary_backend, device_summary=device_cycle_summary)
        end
        if should_record_cycle_metrics(cycle, options.cycles, options.history_stride)
            record = time_block!(timings, :cycle_metrics; synchronize=backend.sync) do
                make_cycle_record_and_deltas!(cycle, deltas.thickness, deltas.wet_mass, deltas.base_mass, final.thickness, final.wet_mass, final.bulk_density, final.base_mass, prev.thickness, prev.wet_mass, prev.base_mass)
            end
            push!(history, record)
            time_block!(timings, :cycle_logging) do
                println(io, cycle_log_line(record))
            end
        end
        if flags.need_last_cycle_smb_delta
            time_block!(timings, :cycle_state_deltas; synchronize=backend.sync) do
                _update_cycle_smb_delta!(deltas.ice_sheet_smb, previous_cycle_smb_ice, domain)
            end
        end
        prev, final = final, prev
    end
    simulation_wall_sec = (time_ns() - simulation_wall_t0) * 1.0e-9

    final_grids = time_block!(timings, :scatter_final_outputs) do
        _scatter_final_grids(prev, domain, deltas, layout)
    end
    layer_grids = if flags.need_layer_outputs
        final_domain = backend.is_gpu ? time_block!(timings, :gpu_transfer) do
            SM.cpu_domain(domain)
        end : domain
        time_block!(timings, :collect_final_layer_grids) do
            collect_final_layer_grids(final_domain, layout.js, layout.is, _grid_shape(layout), final_domain.Ntot)
        end
    else
        _empty_layer_grids()
    end
    monthly_grids = flags.need_monthly_outputs ? time_block!(timings, :aggregate_monthly_outputs) do
        _finalize_monthly_grids(monthly_sums, monthly_count, layout)
    end : empty_monthly_grids()

    summary_path, history_csv_path = _write_optional_outputs!(timings, options, forcing.time_values, ncol, history, :cycles)
    time_block!(timings, :write_netcdf) do
        finalize_case_netcdf!(writer, final_grids, layer_grids, history, monthly_grids, :cycles, length(history), steps_written)
    end

    run_wall_sec = (time_ns() - run_wall_t0) * 1.0e-9
    _print_run_report(io, options, forcing.time_values, history, :cycles, simulation_wall_sec, run_wall_sec, timings; nc_path=nc_path, summary_path=summary_path, history_csv_path=history_csv_path)
    return RunResult(history, :cycles, timings, simulation_wall_sec, run_wall_sec, nc_path, summary_path, history_csv_path, domain, options)
end
