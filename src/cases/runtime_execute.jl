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
            SM.summarize_cycle_state!(summary.thickness, summary.wet_mass, summary.bulk_density, summary.base_mass, domain)
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
        if backend.summary_backend == :kernelabstractions
            device_cycle_summary === nothing && error("`device_summary` must be provided for `backend=:kernelabstractions`.")
            SM.summarize_cycle_state!(
                device_cycle_summary.thickness,
                device_cycle_summary.wet_mass,
                device_cycle_summary.bulk_density,
                device_cycle_summary.base_mass,
                domain;
                backend=backend.summary_backend,
            )
            _copy_summary_fields!(prev, device_cycle_summary, CYCLE_BUFFER_NAMES)
        else
            SM.summarize_cycle_state!(
                prev.thickness,
                prev.wet_mass,
                prev.bulk_density,
                prev.base_mass,
                domain;
                backend=backend.summary_backend,
            )
        end
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
                    if backend.summary_backend == :kernelabstractions
                        device_step_summary === nothing && error("`device_summary` must be provided for `backend=:kernelabstractions`.")
                        SM.summarize_domain_state!(
                            device_step_summary.thickness,
                            device_step_summary.wet_mass,
                            device_step_summary.bulk_density,
                            device_step_summary.base_mass,
                            device_step_summary.smb_ice,
                            device_step_summary.liquid_water,
                            device_step_summary.runoff,
                            domain;
                            backend=backend.summary_backend,
                        )
                        _copy_summary_fields!(step_summary, device_step_summary, SUMMARY_BUFFER_NAMES)
                    else
                        SM.summarize_domain_state!(
                            step_summary.thickness,
                            step_summary.wet_mass,
                            step_summary.bulk_density,
                            step_summary.base_mass,
                            step_summary.smb_ice,
                            step_summary.liquid_water,
                            step_summary.runoff,
                            domain;
                            backend=backend.summary_backend,
                        )
                    end
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
            if backend.summary_backend == :kernelabstractions
                device_cycle_summary === nothing && error("`device_summary` must be provided for `backend=:kernelabstractions`.")
                SM.summarize_cycle_state!(
                    device_cycle_summary.thickness,
                    device_cycle_summary.wet_mass,
                    device_cycle_summary.bulk_density,
                    device_cycle_summary.base_mass,
                    domain;
                    backend=backend.summary_backend,
                )
                _copy_summary_fields!(final, device_cycle_summary, CYCLE_BUFFER_NAMES)
            else
                SM.summarize_cycle_state!(
                    final.thickness,
                    final.wet_mass,
                    final.bulk_density,
                    final.base_mass,
                    domain;
                    backend=backend.summary_backend,
                )
            end
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
