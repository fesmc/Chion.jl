function run_case_cycles_threads_no_netcdf!(
    timings::TimingStats,
    options::RunConfig,
    domain::SM.SnowpackDomain,
    workspaces::AbstractVector{<:SM.StepWorkspace},
    step_fields::SM.SnowpackStepFields,
)
    options.backend == :threads || error("Threaded fast path requires `backend=:threads`.")
    !options.write_netcdf || error("Threaded fast path only applies when NetCDF output is disabled.")
    ntime = size(step_fields.air_temperature, 2)
    ncol = SM.column_count(domain)
    prev = allocate_cycle_summary_buffers(ncol)
    final = allocate_cycle_summary_buffers(ncol)
    time_block!(timings, :summarize_columns_initial) do
        summarize_cycle_columns!(prev, domain; backend=:threads)
    end
    history = NamedTuple[]
    status = :cycles
    last_delta_thickness_vec = fill(NaN, ncol)
    last_delta_wet_mass_vec = fill(NaN, ncol)
    last_delta_base_mass_vec = fill(NaN, ncol)
    simulation_wall_t0 = time_ns()
    for cycle in 1:options.cycles
        time_counted_block!(timings, :model_step_wall, ncol * ntime) do
            SM.step!(domain, step_fields, workspaces)
        end
        time_block!(timings, :summarize_columns_cycle) do
            summarize_cycle_columns!(final, domain; backend=:threads)
        end
        if should_record_cycle_metrics(cycle, options.cycles, options.history_stride)
            record = time_block!(timings, :cycle_metrics) do
                make_cycle_record_and_deltas!(
                    cycle,
                    last_delta_thickness_vec,
                    last_delta_wet_mass_vec,
                    last_delta_base_mass_vec,
                    final.thickness,
                    final.wet_mass,
                    final.bulk_density,
                    final.base_mass,
                    prev.thickness,
                    prev.wet_mass,
                    prev.base_mass,
                )
            end
            push!(history, record)
            time_block!(timings, :cycle_logging) do
                println(cycle_log_line(record))
            end
        end
        prev, final = final, prev
    end
    simulation_wall_sec = (time_ns() - simulation_wall_t0) * 1.0e-9
    return (history=history, status=status, simulation_wall_sec=simulation_wall_sec)
end

function run_case_cycles_gpu_no_netcdf!(
    timings::TimingStats,
    options::RunConfig,
    domain::SM.SnowpackDomain,
    workspace::SM.ColumnarStepWorkspace,
    step_fields::SM.SnowpackStepFields,
)
    options.backend == :gpu || error("GPU fast path requires `backend=:gpu`.")
    !options.write_netcdf || error("GPU fast path only applies when NetCDF output is disabled.")
    ntime = size(step_fields.air_temperature, 2)
    ncol = SM.column_count(domain)
    prev = allocate_cycle_summary_buffers(domain, ncol)
    final = allocate_cycle_summary_buffers(domain, ncol)
    time_block!(timings, :summarize_columns_initial; synchronize=CUDA.synchronize) do
        SM.summarize_cycle_state!(
            prev.thickness,
            prev.wet_mass,
            prev.bulk_density,
            prev.base_mass,
            domain;
            backend=:kernelabstractions,
        )
    end
    history = NamedTuple[]
    status = :cycles
    last_delta_thickness_vec = similar(domain.mass, Float64, ncol)
    last_delta_wet_mass_vec = similar(domain.mass, Float64, ncol)
    last_delta_base_mass_vec = similar(domain.mass, Float64, ncol)
    cycle_metrics_workspace = CycleMetricsWorkspace(domain)
    simulation_wall_t0 = time_ns()
    for cycle in 1:options.cycles
        time_counted_block!(timings, :model_step_wall, ncol * ntime; synchronize=CUDA.synchronize) do
            SM.step!(domain, step_fields, workspace; update_snow_cover=false)
        end
        time_block!(timings, :summarize_columns_cycle; synchronize=CUDA.synchronize) do
            SM.summarize_cycle_state!(
                final.thickness,
                final.wet_mass,
                final.bulk_density,
                final.base_mass,
                domain;
                backend=:kernelabstractions,
            )
        end
        if should_record_cycle_metrics(cycle, options.cycles, options.history_stride)
            record = time_block!(timings, :cycle_metrics; synchronize=CUDA.synchronize) do
                make_cycle_record_and_deltas!(
                    cycle_metrics_workspace,
                    cycle,
                    last_delta_thickness_vec,
                    last_delta_wet_mass_vec,
                    last_delta_base_mass_vec,
                    final.thickness,
                    final.wet_mass,
                    final.bulk_density,
                    final.base_mass,
                    prev.thickness,
                    prev.wet_mass,
                    prev.base_mass,
                )
            end
            push!(history, record)
            time_block!(timings, :cycle_logging) do
                println(cycle_log_line(record))
            end
        end
        prev, final = final, prev
    end
    simulation_wall_sec = (time_ns() - simulation_wall_t0) * 1.0e-9
    return (history=history, status=status, simulation_wall_sec=simulation_wall_sec)
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
    cycles_completed = completed_cycle_count(history, status, options.cycles)
    println(io, "$(options.name) complete.")
    println(io, "Input label     : ", isempty(options.input_label) ? "(not provided)" : options.input_label)
    println(io, "Forcing start   : $(first(time_values))")
    println(io, "Forcing end     : $(last(time_values))")
    println(io, "Backend         : $(String(options.backend))")
    println(io, "Cycles          : $(cycles_completed)")
    println(io, "Status          : $(string(status))")
    println(io, "Cycle metrics   : $(cycle_metrics_schedule_label(options.history_stride))")
    println(io, @sprintf("Simulation wall : %.3f s", simulation_wall_sec))
    println(io, @sprintf("Run wall total  : %.3f s", run_wall_sec))
    if haskey(timings.totals, :model_step_wall)
        println(io, @sprintf("Model step wall : %.3f s", timings.totals[:model_step_wall]))
    end
    if options.write_netcdf
        println(io, "Output NetCDF   : $(abspath(nc_path))")
    else
        println(io, "Output NetCDF   : skipped (--no-nc)")
    end
    if options.write_outputs
        println(io, "History CSV     : $(abspath(history_csv_path))")
        println(io, "Summary         : $(abspath(summary_path))")
    else
        println(io, "File outputs    : skipped (--no-output)")
    end
    print_timing_summary(io, timings; total_wall_sec=run_wall_sec)
    return
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
    !isnothing(layout) && length(layout.js) == ncol || isnothing(layout) || error("Grid-layout point count must match the domain column count.")

    write_final_fields = options.write_netcdf && !isnothing(layout)
    selected = Set(options.netcdf_variables)
    need_step_outputs = options.write_netcdf && case_selected(options, :step)
    need_monthly_outputs = options.write_netcdf && case_selected(options, :monthly)
    need_layer_outputs = options.write_netcdf && case_selected(options, :layers)
    need_last_cycle_smb_delta = options.write_netcdf && :last_cycle_delta_ice_sheet_smb in selected
    need_netcdf_step_diagnostics = need_step_outputs || need_monthly_outputs

    annual_output, unique_month_keys, step_month, nmonth_per_cycle, nmonth_total, month_cycle, month_of_year, source_month_code =
        time_block!(timings, :prepare_output_schedule) do
            annual_output_local = options.write_netcdf ? build_annual_output_schedule(forcing.time_values) : nothing
            unique_month_keys_local = options.write_netcdf ? unique((year(t), month(t)) for t in forcing.time_values) : Tuple{Int, Int}[]
            month_lookup = Dict{Tuple{Int, Int}, Int}()
            if options.write_netcdf
                for (idx, key) in enumerate(unique_month_keys_local)
                    month_lookup[key] = idx
                end
            end
            step_month_local = options.write_netcdf ? [month_lookup[(year(t), month(t))] for t in forcing.time_values] : Int[]
            nmonth_per_cycle_local = options.write_netcdf ? length(unique_month_keys_local) : 0
            nmonth_total_local = options.write_netcdf ? options.cycles * nmonth_per_cycle_local : 0
            month_cycle_local = Int32[]
            month_of_year_local = Int32[]
            source_month_code_local = Int32[]
            if options.write_netcdf
                for cyc in 1:options.cycles, key in unique_month_keys_local
                    push!(month_cycle_local, Int32(cyc))
                    push!(month_of_year_local, Int32(key[2]))
                    push!(source_month_code_local, Int32(key[1] * 100 + key[2]))
                end
            end
            return (
                annual_output_local,
                unique_month_keys_local,
                step_month_local,
                nmonth_per_cycle_local,
                nmonth_total_local,
                month_cycle_local,
                month_of_year_local,
                source_month_code_local,
            )
        end

    initial_thickness_vec = write_final_fields ? Vector{Float64}(undef, ncol) : Float64[]
    if write_final_fields
        time_block!(timings, :prepare_initial_output_fields) do
            @threads :static for idx in 1:ncol
                summary = summarize_column_state(domain, idx)
                initial_thickness_vec[idx] = summary.thickness
            end
        end
    end

    if options.backend == :threads && !options.write_outputs && !options.write_netcdf
        step_fields = SM.SnowpackStepFields(forcing)
        workspaces = time_block!(timings, :create_workspaces) do
            threaded_workspaces(domain)
        end
        sim = redirect_stdout(io) do
            run_case_cycles_threads_no_netcdf!(timings, options, domain, workspaces, step_fields)
        end
        run_wall_sec = (time_ns() - run_wall_t0) * 1.0e-9
        _print_run_report(io, options, forcing.time_values, sim.history, sim.status, sim.simulation_wall_sec, run_wall_sec, timings)
        return RunResult(sim.history, sim.status, timings, sim.simulation_wall_sec, run_wall_sec, "", "", "", domain, options)
    end

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
        if !options.write_outputs && !options.write_netcdf
            sim = redirect_stdout(io) do
                run_case_cycles_gpu_no_netcdf!(
                    timings,
                    options,
                    domain,
                    workspaces,
                    step_fields,
                )
            end
            run_wall_sec = (time_ns() - run_wall_t0) * 1.0e-9
            _print_run_report(io, options, forcing.time_values, sim.history, sim.status, sim.simulation_wall_sec, run_wall_sec, timings)
            return RunResult(sim.history, sim.status, timings, sim.simulation_wall_sec, run_wall_sec, "", "", "", domain, options)
        end
    else
        workspaces = time_block!(timings, :create_workspaces) do
            threaded_workspaces(domain)
        end
    end

    gpu_netcdf_diagnostics = options.backend == :gpu && need_netcdf_step_diagnostics
    gpu_stage_sync = options.backend == :gpu ? CUDA.synchronize : nothing
    prev, final, last_delta_thickness_vec, last_delta_wet_mass_vec, last_delta_base_mass_vec, last_delta_ice_sheet_smb_vec =
        time_block!(timings, :allocate_cycle_buffers) do
            prev_local = options.backend == :gpu ? allocate_cycle_summary_buffers(domain, ncol) : allocate_cycle_summary_buffers(ncol)
            final_local = options.backend == :gpu ? allocate_cycle_summary_buffers(domain, ncol) : allocate_cycle_summary_buffers(ncol)
            last_delta_thickness_vec_local = options.backend == :gpu ? similar(domain.mass, Float64, ncol) : fill(NaN, ncol)
            last_delta_wet_mass_vec_local = options.backend == :gpu ? similar(domain.mass, Float64, ncol) : fill(NaN, ncol)
            last_delta_base_mass_vec_local = options.backend == :gpu ? similar(domain.mass, Float64, ncol) : fill(NaN, ncol)
            last_delta_ice_sheet_smb_vec_local = options.backend == :gpu ? similar(domain.mass, Float64, ncol) : fill(NaN, ncol)
            return (
                prev_local,
                final_local,
                last_delta_thickness_vec_local,
                last_delta_wet_mass_vec_local,
                last_delta_base_mass_vec_local,
                last_delta_ice_sheet_smb_vec_local,
            )
        end
    cycle_metrics_workspace = options.backend == :gpu ? CycleMetricsWorkspace(domain) : nothing
    time_block!(timings, :summarize_columns_initial; synchronize=gpu_stage_sync) do
        if options.backend == :gpu
            SM.summarize_cycle_state!(
                prev.thickness,
                prev.wet_mass,
                prev.bulk_density,
                prev.base_mass,
                domain;
                backend=:kernelabstractions,
            )
        else
            summarize_cycle_columns!(prev, domain; backend=:threads)
        end
    end
    previous_base_mass_vec, previous_smb_ice_vec, previous_runoff_vec, previous_cycle_smb_ice_vec =
        time_block!(timings, :initialize_cycle_tracking) do
            previous_base_mass_vec_local = if gpu_netcdf_diagnostics
                copy(prev.base_mass)
            elseif need_netcdf_step_diagnostics
                copy(prev.base_mass)
            else
                Float64[]
            end
            previous_smb_ice_vec_local = if gpu_netcdf_diagnostics
                copy(domain.smb_ice)
            elseif need_netcdf_step_diagnostics
                copy(domain.smb_ice)
            else
                Float64[]
            end
            previous_runoff_vec_local = if gpu_netcdf_diagnostics
                copy(domain.runoff)
            elseif need_netcdf_step_diagnostics
                copy(domain.runoff)
            else
                Float64[]
            end
            previous_cycle_smb_ice_vec_local = need_last_cycle_smb_delta ? copy(domain.smb_ice) : Float64[]
            return (
                previous_base_mass_vec_local,
                previous_smb_ice_vec_local,
                previous_runoff_vec_local,
                previous_cycle_smb_ice_vec_local,
            )
        end

    initial_thickness = write_final_fields ? time_block!(timings, :prepare_initial_output_fields_grid) do
        scatter_to_grid(initial_thickness_vec, layout.js, layout.is, _grid_shape(layout))
    end : Matrix{Float64}(undef, 0, 0)
    nc_path = resolve_case_netcdf_path(options)
    writer = options.write_netcdf ? time_block!(timings, :init_netcdf) do
        init_case_netcdf(
            nc_path,
            options,
            forcing.time_values,
            domain.Ntot,
            layout,
            initial_thickness,
            month_cycle,
            month_of_year,
            source_month_code,
            annual_output.source_indices,
            annual_output.source_codes,
        )
    end : nothing
    write_step_fields = options.write_netcdf && !isnothing(writer) &&
        (haskey(writer.vars, :step_export_to_ice) || haskey(writer.vars, :step_ice_sheet_smb))
    step_summary, device_step_summary, annual_export_to_ice, annual_ice_sheet_smb, steps_written, history, status,
    monthly_sum_thickness, monthly_sum_wet_mass, monthly_sum_bulk_density, monthly_sum_base_mass,
    monthly_sum_ice_sheet_smb, monthly_sum_export, monthly_sum_net_ice_sheet_forcing, monthly_sum_runoff, monthly_count =
        time_block!(timings, :allocate_output_buffers) do
            step_summary_local = need_netcdf_step_diagnostics && !gpu_netcdf_diagnostics ? allocate_summary_buffers(ncol) : nothing
            device_step_summary_local = gpu_netcdf_diagnostics ? allocate_summary_buffers(domain, ncol) : nothing
            annual_export_to_ice_local = if gpu_netcdf_diagnostics && need_step_outputs
                CUDA.zeros(Float64, ncol)
            elseif need_step_outputs
                zeros(Float64, ncol)
            else
                Float64[]
            end
            annual_ice_sheet_smb_local = if gpu_netcdf_diagnostics && need_step_outputs
                CUDA.zeros(Float64, ncol)
            elseif need_step_outputs
                zeros(Float64, ncol)
            else
                Float64[]
            end
            history_local = NamedTuple[]
            status_local = :cycles
            monthly_sum_thickness_local = if gpu_netcdf_diagnostics && need_monthly_outputs
                CUDA.zeros(Float64, nmonth_total, ncol)
            elseif need_monthly_outputs
                zeros(Float64, nmonth_total, ncol)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_wet_mass_local = if gpu_netcdf_diagnostics && need_monthly_outputs
                CUDA.zeros(Float64, nmonth_total, ncol)
            elseif need_monthly_outputs
                zeros(Float64, nmonth_total, ncol)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_bulk_density_local = if gpu_netcdf_diagnostics && need_monthly_outputs
                CUDA.zeros(Float64, nmonth_total, ncol)
            elseif need_monthly_outputs
                zeros(Float64, nmonth_total, ncol)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_base_mass_local = if gpu_netcdf_diagnostics && need_monthly_outputs
                CUDA.zeros(Float64, nmonth_total, ncol)
            elseif need_monthly_outputs
                zeros(Float64, nmonth_total, ncol)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_ice_sheet_smb_local = if gpu_netcdf_diagnostics && need_monthly_outputs
                CUDA.zeros(Float64, nmonth_total, ncol)
            elseif need_monthly_outputs
                zeros(Float64, nmonth_total, ncol)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_export_local = if gpu_netcdf_diagnostics && need_monthly_outputs
                CUDA.zeros(Float64, nmonth_total, ncol)
            elseif need_monthly_outputs
                zeros(Float64, nmonth_total, ncol)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_net_ice_sheet_forcing_local = if gpu_netcdf_diagnostics && need_monthly_outputs
                CUDA.zeros(Float64, nmonth_total, ncol)
            elseif need_monthly_outputs
                zeros(Float64, nmonth_total, ncol)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_runoff_local = if gpu_netcdf_diagnostics && need_monthly_outputs
                CUDA.zeros(Float64, nmonth_total, ncol)
            elseif need_monthly_outputs
                zeros(Float64, nmonth_total, ncol)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_count_local = need_monthly_outputs ? zeros(Int32, nmonth_total) : Int32[]
            return (
                step_summary_local,
                device_step_summary_local,
                annual_export_to_ice_local,
                annual_ice_sheet_smb_local,
                0,
                history_local,
                status_local,
                monthly_sum_thickness_local,
                monthly_sum_wet_mass_local,
                monthly_sum_bulk_density_local,
                monthly_sum_base_mass_local,
                monthly_sum_ice_sheet_smb_local,
                monthly_sum_export_local,
                monthly_sum_net_ice_sheet_forcing_local,
                monthly_sum_runoff_local,
                monthly_count_local,
            )
        end

    simulation_wall_t0 = time_ns()
    ntime = length(forcing.time_values)
    for cycle in 1:options.cycles
        for t in 1:ntime
            month_idx = need_monthly_outputs ? (cycle - 1) * nmonth_per_cycle + step_month[t] : 0
            if options.backend == :gpu
                time_counted_block!(timings, :model_step_wall, ncol; synchronize=gpu_stage_sync) do
                    SM.step!(domain, step_fields, t, workspaces)
                end
                if need_netcdf_step_diagnostics
                    time_counted_block!(timings, :step_diagnostics, ncol; synchronize=gpu_stage_sync) do
                        if gpu_netcdf_diagnostics
                            SM.summarize_domain_state!(
                                device_step_summary.thickness,
                                device_step_summary.wet_mass,
                                device_step_summary.bulk_density,
                                device_step_summary.base_mass,
                                device_step_summary.smb_ice,
                                device_step_summary.liquid_water,
                                device_step_summary.runoff,
                                domain;
                                backend=:kernelabstractions,
                            )
                            current_base = device_step_summary.base_mass
                            current_smb_ice = device_step_summary.smb_ice
                            current_runoff = device_step_summary.runoff
                            if need_monthly_outputs
                                monthly_sum_thickness[month_idx, :] .+= device_step_summary.thickness
                                monthly_sum_wet_mass[month_idx, :] .+= device_step_summary.wet_mass
                                monthly_sum_bulk_density[month_idx, :] .+= device_step_summary.bulk_density
                                monthly_sum_base_mass[month_idx, :] .+= current_base
                                monthly_sum_ice_sheet_smb[month_idx, :] .+= current_smb_ice .- previous_smb_ice_vec
                                monthly_sum_export[month_idx, :] .+= current_base .- previous_base_mass_vec
                                monthly_sum_net_ice_sheet_forcing[month_idx, :] .+= current_smb_ice .- previous_smb_ice_vec
                                monthly_sum_runoff[month_idx, :] .+= current_runoff .- previous_runoff_vec
                            end
                            if need_step_outputs
                                annual_export_to_ice .+= current_base .- previous_base_mass_vec
                                annual_ice_sheet_smb .+= current_smb_ice .- previous_smb_ice_vec
                            end
                            previous_base_mass_vec .= current_base
                            previous_smb_ice_vec .= current_smb_ice
                            previous_runoff_vec .= current_runoff
                        else
                            summarize_columns!(step_summary, domain; backend=:kernelabstractions, device_summary=device_step_summary)
                            current_base = step_summary.base_mass
                            current_smb_ice = step_summary.smb_ice
                            current_runoff = step_summary.runoff
                            if need_monthly_outputs
                                monthly_sum_thickness[month_idx, :] .+= step_summary.thickness
                                monthly_sum_wet_mass[month_idx, :] .+= step_summary.wet_mass
                                monthly_sum_bulk_density[month_idx, :] .+= step_summary.bulk_density
                                monthly_sum_base_mass[month_idx, :] .+= current_base
                                monthly_sum_ice_sheet_smb[month_idx, :] .+= current_smb_ice .- previous_smb_ice_vec
                                monthly_sum_export[month_idx, :] .+= current_base .- previous_base_mass_vec
                                monthly_sum_net_ice_sheet_forcing[month_idx, :] .+= current_smb_ice .- previous_smb_ice_vec
                                monthly_sum_runoff[month_idx, :] .+= current_runoff .- previous_runoff_vec
                            end
                            if need_step_outputs
                                annual_export_to_ice .+= current_base .- previous_base_mass_vec
                                annual_ice_sheet_smb .+= current_smb_ice .- previous_smb_ice_vec
                            end
                            previous_base_mass_vec .= current_base
                            previous_smb_ice_vec .= current_smb_ice
                            previous_runoff_vec .= current_runoff
                        end
                    end
                end
            elseif need_netcdf_step_diagnostics
                t0 = time_ns()
                SM.step!(domain, step_fields, t, workspaces)
                add_timing!(timings, :model_step_wall, (time_ns() - t0) * 1.0e-9, ncol)

                diag_t0 = time_ns()
                summarize_columns!(step_summary, domain; backend=:threads)
                @threads :static for idx in 1:ncol
                    current_base = step_summary.base_mass[idx]
                    current_smb_ice = step_summary.smb_ice[idx]
                    current_runoff = step_summary.runoff[idx]
                    if need_monthly_outputs
                        monthly_sum_thickness[month_idx, idx] += step_summary.thickness[idx]
                        monthly_sum_wet_mass[month_idx, idx] += step_summary.wet_mass[idx]
                        monthly_sum_bulk_density[month_idx, idx] += step_summary.bulk_density[idx]
                        monthly_sum_base_mass[month_idx, idx] += current_base
                        monthly_sum_ice_sheet_smb[month_idx, idx] += current_smb_ice - previous_smb_ice_vec[idx]
                        monthly_sum_export[month_idx, idx] += current_base - previous_base_mass_vec[idx]
                        monthly_sum_net_ice_sheet_forcing[month_idx, idx] += current_smb_ice - previous_smb_ice_vec[idx]
                        monthly_sum_runoff[month_idx, idx] += current_runoff - previous_runoff_vec[idx]
                    end
                    if need_step_outputs
                        annual_export_to_ice[idx] += current_base - previous_base_mass_vec[idx]
                        annual_ice_sheet_smb[idx] += current_smb_ice - previous_smb_ice_vec[idx]
                    end
                    previous_base_mass_vec[idx] = current_base
                    previous_smb_ice_vec[idx] = current_smb_ice
                    previous_runoff_vec[idx] = current_runoff
                end
                add_timing!(timings, :step_diagnostics, (time_ns() - diag_t0) * 1.0e-9, ncol)
            else
                t0 = time_ns()
                SM.step!(domain, step_fields, t, workspaces)
                add_timing!(timings, :model_step_wall, (time_ns() - t0) * 1.0e-9, ncol)
            end
            if need_monthly_outputs
                monthly_count[month_idx] += 1
            end
            if need_step_outputs && annual_output.write_output[t]
                steps_written += 1
                if write_step_fields
                    step_export_to_ice_vec, step_ice_sheet_smb_vec, step_export_to_ice, step_ice_sheet_smb =
                        time_block!(timings, :step_output_prepare) do
                            step_export_to_ice_vec_local = gpu_netcdf_diagnostics ? Array(annual_export_to_ice) : annual_export_to_ice
                            step_ice_sheet_smb_vec_local = gpu_netcdf_diagnostics ? Array(annual_ice_sheet_smb) : annual_ice_sheet_smb
                            step_export_to_ice_local = scatter_to_grid(step_export_to_ice_vec_local, layout.js, layout.is, _grid_shape(layout))
                            step_ice_sheet_smb_local = scatter_to_grid(step_ice_sheet_smb_vec_local, layout.js, layout.is, _grid_shape(layout))
                            return (
                                step_export_to_ice_vec_local,
                                step_ice_sheet_smb_vec_local,
                                step_export_to_ice_local,
                                step_ice_sheet_smb_local,
                            )
                        end
                    time_block!(timings, :step_output_write) do
                        maybe_write_step_output!(writer, steps_written, :step_export_to_ice, step_export_to_ice)
                        maybe_write_step_output!(writer, steps_written, :step_ice_sheet_smb, step_ice_sheet_smb)
                    end
                end
                fill!(annual_export_to_ice, 0.0)
                fill!(annual_ice_sheet_smb, 0.0)
            end
        end

        if options.backend == :gpu
            time_block!(timings, :summarize_columns_cycle; synchronize=gpu_stage_sync) do
                SM.summarize_cycle_state!(
                    final.thickness,
                    final.wet_mass,
                    final.bulk_density,
                    final.base_mass,
                    domain;
                    backend=:kernelabstractions,
                )
            end
        else
            time_block!(timings, :summarize_columns_cycle) do
                summarize_cycle_columns!(final, domain; backend=:threads)
            end
        end
        if should_record_cycle_metrics(cycle, options.cycles, options.history_stride)
            record = if options.backend == :gpu
                time_block!(timings, :cycle_metrics; synchronize=gpu_stage_sync) do
                    record_local = make_cycle_record_and_deltas!(
                        cycle_metrics_workspace,
                        cycle,
                        last_delta_thickness_vec,
                        last_delta_wet_mass_vec,
                        last_delta_base_mass_vec,
                        final.thickness,
                        final.wet_mass,
                        final.bulk_density,
                        final.base_mass,
                        prev.thickness,
                        prev.wet_mass,
                        prev.base_mass,
                    )
                    if need_last_cycle_smb_delta
                        last_delta_ice_sheet_smb_vec .= domain.smb_ice .- previous_cycle_smb_ice_vec
                        previous_cycle_smb_ice_vec .= domain.smb_ice
                    end
                    record_local
                end
            else
                time_block!(timings, :cycle_metrics) do
                    record_local = make_cycle_record_and_deltas!(
                        cycle,
                        last_delta_thickness_vec,
                        last_delta_wet_mass_vec,
                        last_delta_base_mass_vec,
                        final.thickness,
                        final.wet_mass,
                        final.bulk_density,
                        final.base_mass,
                        prev.thickness,
                        prev.wet_mass,
                        prev.base_mass,
                    )
                    if need_last_cycle_smb_delta
                        current_smb_ice_vec = copy(domain.smb_ice)
                        last_delta_ice_sheet_smb_vec .= current_smb_ice_vec .- previous_cycle_smb_ice_vec
                        previous_cycle_smb_ice_vec .= current_smb_ice_vec
                    end
                    record_local
                end
            end
            push!(history, record)
            time_block!(timings, :cycle_logging) do
                println(io, cycle_log_line(record))
            end
        elseif options.backend == :gpu
            time_block!(timings, :cycle_state_deltas; synchronize=gpu_stage_sync) do
                last_delta_thickness_vec .= final.thickness .- prev.thickness
                last_delta_wet_mass_vec .= final.wet_mass .- prev.wet_mass
                last_delta_base_mass_vec .= final.base_mass .- prev.base_mass
                if need_last_cycle_smb_delta
                    last_delta_ice_sheet_smb_vec .= domain.smb_ice .- previous_cycle_smb_ice_vec
                    previous_cycle_smb_ice_vec .= domain.smb_ice
                end
            end
        else
            time_block!(timings, :cycle_state_deltas) do
                last_delta_thickness_vec .= final.thickness .- prev.thickness
                last_delta_wet_mass_vec .= final.wet_mass .- prev.wet_mass
                last_delta_base_mass_vec .= final.base_mass .- prev.base_mass
                if need_last_cycle_smb_delta
                    current_smb_ice_vec = copy(domain.smb_ice)
                    last_delta_ice_sheet_smb_vec .= current_smb_ice_vec .- previous_cycle_smb_ice_vec
                    previous_cycle_smb_ice_vec .= current_smb_ice_vec
                end
            end
        end
        prev, final = final, prev
    end
    simulation_wall_sec = (time_ns() - simulation_wall_t0) * 1.0e-9
    final_state = if status == :cycles
        history[end].cycle == options.cycles ? final : prev
    else
        final
    end

    final_thickness = Matrix{Float64}(undef, 0, 0)
    final_wet_mass = Matrix{Float64}(undef, 0, 0)
    final_bulk_density = Matrix{Float64}(undef, 0, 0)
    final_base_mass = Matrix{Float64}(undef, 0, 0)
    final_ice_sheet_smb = Matrix{Float64}(undef, 0, 0)
    final_runoff = Matrix{Float64}(undef, 0, 0)
    last_delta_thickness = Matrix{Float64}(undef, 0, 0)
    last_delta_wet_mass = Matrix{Float64}(undef, 0, 0)
    last_delta_base_mass = Matrix{Float64}(undef, 0, 0)
    last_delta_ice_sheet_smb = Matrix{Float64}(undef, 0, 0)
    if write_final_fields
        final_state_host, final_smb_ice_vec, final_runoff_vec, last_delta_thickness_host, last_delta_wet_mass_host,
        last_delta_base_mass_host, last_delta_ice_sheet_smb_host = time_block!(timings, :finalize_state_transfer) do
            final_state_host_local = options.backend == :gpu ? cpu_cycle_summary(final_state) : final_state
            final_smb_ice_vec_local = options.backend == :gpu ? Array(domain.smb_ice) : copy(domain.smb_ice)
            final_runoff_vec_local = options.backend == :gpu ? Array(domain.runoff) : copy(domain.runoff)
            last_delta_thickness_host_local = options.backend == :gpu ? Array(last_delta_thickness_vec) : last_delta_thickness_vec
            last_delta_wet_mass_host_local = options.backend == :gpu ? Array(last_delta_wet_mass_vec) : last_delta_wet_mass_vec
            last_delta_base_mass_host_local = options.backend == :gpu ? Array(last_delta_base_mass_vec) : last_delta_base_mass_vec
            last_delta_ice_sheet_smb_host_local = if need_last_cycle_smb_delta
                options.backend == :gpu ? Array(last_delta_ice_sheet_smb_vec) : last_delta_ice_sheet_smb_vec
            else
                fill(NaN, ncol)
            end
            return (
                final_state_host_local,
                final_smb_ice_vec_local,
                final_runoff_vec_local,
                last_delta_thickness_host_local,
                last_delta_wet_mass_host_local,
                last_delta_base_mass_host_local,
                last_delta_ice_sheet_smb_host_local,
            )
        end
        final_thickness, final_wet_mass, final_bulk_density, final_base_mass, final_ice_sheet_smb, final_runoff,
        last_delta_thickness, last_delta_wet_mass, last_delta_base_mass, last_delta_ice_sheet_smb =
            time_block!(timings, :scatter_final_outputs) do
                return (
                    scatter_to_grid(final_state_host.thickness, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(final_state_host.wet_mass, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(final_state_host.bulk_density, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(final_state_host.base_mass, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(final_smb_ice_vec, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(final_runoff_vec, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(last_delta_thickness_host, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(last_delta_wet_mass_host, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(last_delta_base_mass_host, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(last_delta_ice_sheet_smb_host, layout.js, layout.is, _grid_shape(layout)),
                )
            end
    end

    layer_grids = _empty_layer_grids()
    monthly_mean_thickness_grid = _empty_monthly_grid()
    monthly_mean_wet_mass_grid = _empty_monthly_grid()
    monthly_mean_bulk_density_grid = _empty_monthly_grid()
    monthly_mean_base_mass_grid = _empty_monthly_grid()
    monthly_mean_ice_sheet_smb_grid = _empty_monthly_grid()
    monthly_export_to_ice_grid = _empty_monthly_grid()
    monthly_net_ice_sheet_forcing_grid = _empty_monthly_grid()
    monthly_runoff_grid = _empty_monthly_grid()
    if options.write_netcdf
        if need_layer_outputs
            final_domain = options.backend == :gpu ? time_block!(timings, :gpu_transfer) do
                SM.cpu_domain(domain)
            end : domain
            layer_grids = time_block!(timings, :collect_final_layer_grids) do
                collect_final_layer_grids(final_domain, layout.js, layout.is, _grid_shape(layout), final_domain.Ntot)
            end
        end
        if need_monthly_outputs
            monthly_sum_thickness_host = gpu_netcdf_diagnostics ? time_block!(timings, :gpu_transfer) do
                Array(monthly_sum_thickness)
            end : monthly_sum_thickness
            monthly_sum_wet_mass_host = gpu_netcdf_diagnostics ? time_block!(timings, :gpu_transfer) do
                Array(monthly_sum_wet_mass)
            end : monthly_sum_wet_mass
            monthly_sum_bulk_density_host = gpu_netcdf_diagnostics ? time_block!(timings, :gpu_transfer) do
                Array(monthly_sum_bulk_density)
            end : monthly_sum_bulk_density
            monthly_sum_base_mass_host = gpu_netcdf_diagnostics ? time_block!(timings, :gpu_transfer) do
                Array(monthly_sum_base_mass)
            end : monthly_sum_base_mass
            monthly_sum_ice_sheet_smb_host = gpu_netcdf_diagnostics ? time_block!(timings, :gpu_transfer) do
                Array(monthly_sum_ice_sheet_smb)
            end : monthly_sum_ice_sheet_smb
            monthly_sum_export_host = gpu_netcdf_diagnostics ? time_block!(timings, :gpu_transfer) do
                Array(monthly_sum_export)
            end : monthly_sum_export
            monthly_sum_net_ice_sheet_forcing_host = gpu_netcdf_diagnostics ? time_block!(timings, :gpu_transfer) do
                Array(monthly_sum_net_ice_sheet_forcing)
            end : monthly_sum_net_ice_sheet_forcing
            monthly_sum_runoff_host = gpu_netcdf_diagnostics ? time_block!(timings, :gpu_transfer) do
                Array(monthly_sum_runoff)
            end : monthly_sum_runoff
            monthly_mean_thickness = similar(monthly_sum_thickness_host)
            monthly_mean_wet_mass = similar(monthly_sum_wet_mass_host)
            monthly_mean_bulk_density = similar(monthly_sum_bulk_density_host)
            monthly_mean_base_mass = similar(monthly_sum_base_mass_host)
            monthly_mean_ice_sheet_smb = similar(monthly_sum_ice_sheet_smb_host)
            monthly_export_to_ice = similar(monthly_sum_export_host)
            monthly_net_ice_sheet_forcing = similar(monthly_sum_net_ice_sheet_forcing_host)
            monthly_runoff = similar(monthly_sum_runoff_host)
            monthly_mean_thickness_grid, monthly_mean_wet_mass_grid, monthly_mean_bulk_density_grid,
            monthly_mean_base_mass_grid, monthly_mean_ice_sheet_smb_grid, monthly_export_to_ice_grid,
            monthly_net_ice_sheet_forcing_grid, monthly_runoff_grid = time_block!(timings, :aggregate_monthly_outputs) do
                @inbounds for m in 1:nmonth_total
                    c = max(monthly_count[m], 1)
                    monthly_mean_thickness[m, :] .= monthly_sum_thickness_host[m, :] ./ c
                    monthly_mean_wet_mass[m, :] .= monthly_sum_wet_mass_host[m, :] ./ c
                    monthly_mean_bulk_density[m, :] .= monthly_sum_bulk_density_host[m, :] ./ c
                    monthly_mean_base_mass[m, :] .= monthly_sum_base_mass_host[m, :]
                    monthly_mean_ice_sheet_smb[m, :] .= monthly_sum_ice_sheet_smb_host[m, :]
                    monthly_export_to_ice[m, :] .= monthly_sum_export_host[m, :]
                    monthly_net_ice_sheet_forcing[m, :] .= monthly_sum_net_ice_sheet_forcing_host[m, :]
                    monthly_runoff[m, :] .= monthly_sum_runoff_host[m, :]
                end
                return (
                    monthly_vectors_to_grids(monthly_mean_thickness, layout.js, layout.is, _grid_shape(layout)),
                    monthly_vectors_to_grids(monthly_mean_wet_mass, layout.js, layout.is, _grid_shape(layout)),
                    monthly_vectors_to_grids(monthly_mean_bulk_density, layout.js, layout.is, _grid_shape(layout)),
                    monthly_vectors_to_grids(monthly_mean_base_mass, layout.js, layout.is, _grid_shape(layout)),
                    monthly_vectors_to_grids(monthly_mean_ice_sheet_smb, layout.js, layout.is, _grid_shape(layout)),
                    monthly_vectors_to_grids(monthly_export_to_ice, layout.js, layout.is, _grid_shape(layout)),
                    monthly_vectors_to_grids(monthly_net_ice_sheet_forcing, layout.js, layout.is, _grid_shape(layout)),
                    monthly_vectors_to_grids(monthly_runoff, layout.js, layout.is, _grid_shape(layout)),
                )
            end
        end
    end

    summary_path = ""
    history_csv_path = ""
    if options.write_outputs
        mkpath(options.output_dir)
        summary_path = joinpath(options.output_dir, "$(options.name)_summary.txt")
        history_csv_path = joinpath(options.output_dir, "$(options.name)_history.csv")
        time_block!(timings, :write_summary_text) do
            write_case_summary(summary_path, options, forcing.time_values, ncol, history, status, timings)
        end
        time_block!(timings, :write_history_csv) do
            write_case_history_csv(history_csv_path, history)
        end
    end
    if options.write_netcdf
        time_block!(timings, :write_netcdf) do
            finalize_case_netcdf!(
                writer,
                final_thickness,
                final_wet_mass,
                final_bulk_density,
                final_base_mass,
                final_ice_sheet_smb,
                last_delta_thickness,
                last_delta_wet_mass,
                last_delta_base_mass,
                last_delta_ice_sheet_smb,
                final_runoff,
                layer_grids,
                history,
                monthly_mean_thickness_grid,
                monthly_mean_wet_mass_grid,
                monthly_mean_bulk_density_grid,
                monthly_mean_base_mass_grid,
                monthly_mean_ice_sheet_smb_grid,
                monthly_export_to_ice_grid,
                monthly_net_ice_sheet_forcing_grid,
                monthly_runoff_grid,
                status,
                length(history),
                steps_written,
            )
        end
    end

    run_wall_sec = (time_ns() - run_wall_t0) * 1.0e-9
    _print_run_report(
        io,
        options,
        forcing.time_values,
        history,
        status,
        simulation_wall_sec,
        run_wall_sec,
        timings;
        nc_path=nc_path,
        summary_path=summary_path,
        history_csv_path=history_csv_path,
    )
    return RunResult(history, status, timings, simulation_wall_sec, run_wall_sec, options.write_netcdf ? nc_path : "", summary_path, history_csv_path, domain, options)
end
