@inline _backend_info(options::RunConfig) = (
    is_gpu=options.backend == :gpu,
    summary_backend=options.backend == :gpu ? :kernelabstractions : :threads,
    sync=options.backend == :gpu ? CUDA.synchronize : nothing,
)

@inline _step_cycle!(domain, step_fields, workspaces, options::RunConfig) =
    options.backend == :gpu ? step!(domain, step_fields, workspaces; update_snow_cover=false) : step!(domain, step_fields, workspaces)

@inline _step_timestep!(domain, step_fields, t::Int, workspaces, ::RunConfig) = step!(domain, step_fields, t, workspaces)

function _prepare_backend!(timings::StepTimingStats, options::RunConfig, domain::SnowpackDomain, forcing::ForcingData)
    step_fields = SnowpackStepFields(forcing)
    if options.backend == :gpu
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

function _write_optional_outputs!(
    timings::StepTimingStats,
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

function _print_run_report(
    io::IO,
    options::RunConfig,
    time_values::Vector{DateTime},
    history::Vector{NamedTuple},
    status::Symbol,
    simulation_wall_sec::Float64,
    run_wall_sec::Float64,
    timings::StepTimingStats;
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
    domain::SnowpackDomain,
    forcing::ForcingData;
    layout::Union{Nothing, GridLayout}=nothing,
    options::RunConfig=RunConfig(),
    io::IO=stdout,
    timings::StepTimingStats=StepTimingStats(),
    run_wall_t0::Integer=time_ns(),
)
    ncol = column_count(domain)
    size(forcing.air_temperature, 1) == ncol || error("Forcing column count must match the domain column count.")
    options.write_netcdf && isnothing(layout) && error("NetCDF output requires a grid layout.")
    isnothing(layout) || length(layout.js) == ncol || error("Grid-layout point count must match the domain column count.")

    flags = _case_flags(options, layout)
    initial_thickness_vec = if flags.write_final_fields
        summary = allocate_cycle_summary_buffers(ncol)
        time_block!(timings, :prepare_initial_output_fields) do
            summarize_cycle_state!(summary.thickness, summary.wet_mass, summary.bulk_density, summary.base_mass, domain)
        end
        copy(summary.thickness)
    else
        Float64[]
    end

    domain, step_fields, workspaces = _prepare_backend!(timings, options, domain, forcing)
    backend = _backend_info(options)
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
    end
    write_step_fields = writer !== nothing && any(key -> haskey(writer.vars, key), CASE_OUTPUT_GROUPS.step)

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
    monthly_total = flags.need_monthly_outputs && schedule !== nothing ? schedule.nmonth_total : 0
    monthly_sums = _allocate_monthly_sums(flags.need_monthly_outputs, monthly_total, ncol)
    monthly_count = flags.need_monthly_outputs ? zeros(Int32, monthly_total) : Int32[]
    step_vectors = _allocate_step_vectors(flags.need_step_outputs, ncol)

    time_block!(timings, :summarize_columns_initial; synchronize=backend.sync) do
        if backend.summary_backend == :kernelabstractions
            device_cycle_summary === nothing && error("`device_summary` must be provided for `backend=:kernelabstractions`.")
            summarize_cycle_state!(
                device_cycle_summary.thickness,
                device_cycle_summary.wet_mass,
                device_cycle_summary.bulk_density,
                device_cycle_summary.base_mass,
                domain;
                backend=backend.summary_backend,
            )
            _copy_summary_fields!(prev, device_cycle_summary, CYCLE_BUFFER_NAMES)
        else
            summarize_cycle_state!(
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
                        summarize_domain_state!(
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
                        summarize_domain_state!(
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
                summarize_cycle_state!(
                    device_cycle_summary.thickness,
                    device_cycle_summary.wet_mass,
                    device_cycle_summary.bulk_density,
                    device_cycle_summary.base_mass,
                    domain;
                    backend=backend.summary_backend,
                )
                _copy_summary_fields!(final, device_cycle_summary, CYCLE_BUFFER_NAMES)
            else
                summarize_cycle_state!(
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

    summary_path, history_csv_path = _write_optional_outputs!(timings, options, forcing.time_values, ncol, history, :cycles)
    if writer !== nothing
        final_grids = time_block!(timings, :scatter_final_outputs) do
            _scatter_final_grids(prev, domain, deltas, layout)
        end
        layer_grids = if flags.need_layer_outputs
            final_domain = backend.is_gpu ? time_block!(timings, :gpu_transfer) do
                cpu_domain(domain)
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
        time_block!(timings, :write_netcdf) do
            finalize_case_netcdf!(writer, final_grids, layer_grids, history, monthly_grids, :cycles, length(history), steps_written)
        end
    end

    run_wall_sec = (time_ns() - run_wall_t0) * 1.0e-9
    _print_run_report(io, options, forcing.time_values, history, :cycles, simulation_wall_sec, run_wall_sec, timings; nc_path=nc_path, summary_path=summary_path, history_csv_path=history_csv_path)
    return RunResult(history, :cycles, timings, simulation_wall_sec, run_wall_sec, nc_path, summary_path, history_csv_path, domain, options)
end
function scatter_to_grid(values::Vector{Float64}, js::Vector{Int}, is::Vector{Int}, grid_shape::Tuple{Int, Int})
    out = fill(NaN, grid_shape)
    @inbounds for idx in eachindex(values)
        out[js[idx], is[idx]] = values[idx]
    end
    return out
end

function monthly_vectors_to_grids(values::Matrix{Float64}, js::Vector{Int}, is::Vector{Int}, grid_shape::Tuple{Int, Int})
    nmonth, nvalid = size(values)
    ny, nx = grid_shape
    out = fill(NaN, nmonth, ny, nx)
    @inbounds for m in 1:nmonth, idx in 1:nvalid
        out[m, js[idx], is[idx]] = values[m, idx]
    end
    return out
end

@inline function case_selected(options::RunConfig, group::Symbol)
    group_vars = getproperty(CASE_OUTPUT_GROUPS, group)
    return any(var -> var in group_vars, options.netcdf_variables)
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

function _finite_mean(data)
    total = 0.0
    count = 0
    for value in data
        isfinite(value) || continue
        total += value
        count += 1
    end
    return count == 0 ? NaN : total / count
end

function _delta_stats(data)
    signed_total = 0.0
    abs_total = 0.0
    max_abs = 0.0
    count = 0
    for value in data
        isfinite(value) || continue
        abs_value = abs(value)
        signed_total += value
        abs_total += abs_value
        max_abs = max(max_abs, abs_value)
        count += 1
    end
    return (
        mean_signed=count == 0 ? NaN : signed_total / count,
        mean_abs=count == 0 ? NaN : abs_total / count,
        max_abs=count == 0 ? NaN : max_abs,
    )
end

function make_cycle_record_and_deltas!(
    cycle::Int,
    delta_thickness,
    delta_wet_mass,
    delta_base_mass,
    thickness,
    wet_mass,
    bulk_density,
    base_mass,
    prev_thickness,
    prev_wet_mass,
    prev_base_mass,
)
    delta_thickness .= thickness .- prev_thickness
    delta_wet_mass .= wet_mass .- prev_wet_mass
    delta_base_mass .= base_mass .- prev_base_mass
    dth = _delta_stats(_host_vector(delta_thickness))
    dwet = _delta_stats(_host_vector(delta_wet_mass))
    dbase = _delta_stats(_host_vector(delta_base_mass))
    return (
        cycle=cycle,
        mean_thickness=_finite_mean(_host_vector(thickness)),
        mean_wet_mass=_finite_mean(_host_vector(wet_mass)),
        mean_bulk_density=_finite_mean(_host_vector(bulk_density)),
        mean_base_mass=_finite_mean(_host_vector(base_mass)),
        mean_signed_delta_thickness=dth.mean_signed,
        mean_abs_delta_thickness=dth.mean_abs,
        max_abs_delta_thickness=dth.max_abs,
        mean_signed_delta_wet_mass=dwet.mean_signed,
        mean_abs_delta_wet_mass=dwet.mean_abs,
        max_abs_delta_wet_mass=dwet.max_abs,
        mean_signed_delta_base_mass=dbase.mean_signed,
        mean_abs_delta_base_mass=dbase.mean_abs,
        max_abs_delta_base_mass=dbase.max_abs,
    )
end

function _copy_summary_fields!(summary, device_summary, fields)
    for field in fields
        copyto!(getfield(summary, field), getfield(device_summary, field))
    end
    return summary
end

function build_annual_output_schedule(time_values::Vector{DateTime})
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
    return (write_output=write_output, output_slot=output_slot, source_indices=source_indices, source_codes=source_codes, years=years)
end

function cycle_log_line(record)
    return @sprintf(
        "cycle=%d mean_th=%.5f m mean_swe=%.5f mmWE mean_base=%.5f mmWE mean_abs_dth=%.5f m mean_abs_dswe=%.5f mmWE mean_abs_dbase=%.5f mmWE",
        record.cycle,
        record.mean_thickness,
        record.mean_wet_mass,
        record.mean_base_mass,
        record.mean_abs_delta_thickness,
        record.mean_abs_delta_wet_mass,
        record.mean_abs_delta_base_mass,
    )
end

const HISTORY_CSV_SPECS = (
    (key=:cycle, label="cycle", integer=true),
    (key=:mean_thickness, label="mean_thickness_m", integer=false),
    (key=:mean_wet_mass, label="mean_wet_mass_mmwe", integer=false),
    (key=:mean_bulk_density, label="mean_bulk_density_kgm3", integer=false),
    (key=:mean_base_mass, label="mean_base_mass_mmwe", integer=false),
    (key=:mean_signed_delta_thickness, label="mean_signed_delta_thickness_m", integer=false),
    (key=:mean_abs_delta_thickness, label="mean_abs_delta_thickness_m", integer=false),
    (key=:max_abs_delta_thickness, label="max_abs_delta_thickness_m", integer=false),
    (key=:mean_signed_delta_wet_mass, label="mean_signed_delta_wet_mass_mmwe", integer=false),
    (key=:mean_abs_delta_wet_mass, label="mean_abs_delta_wet_mass_mmwe", integer=false),
    (key=:max_abs_delta_wet_mass, label="max_abs_delta_wet_mass_mmwe", integer=false),
    (key=:mean_signed_delta_base_mass, label="mean_signed_delta_base_mass_mmwe", integer=false),
    (key=:mean_abs_delta_base_mass, label="mean_abs_delta_base_mass_mmwe", integer=false),
    (key=:max_abs_delta_base_mass, label="max_abs_delta_base_mass_mmwe", integer=false),
)

function write_case_history_csv(out_path::AbstractString, history::Vector{NamedTuple})
    mkpath(dirname(out_path))
    open(out_path, "w") do io
        println(io, join((spec.label for spec in HISTORY_CSV_SPECS), ","))
        for rec in history
            values = String[]
            for spec in HISTORY_CSV_SPECS
                value = getfield(rec, spec.key)
                push!(values, spec.integer ? string(value) : @sprintf("%.10f", value))
            end
            println(io, join(values, ","))
        end
    end
end

const SUMMARY_REPORT_SPECS = (
    (title="Final domain means", fields=(
        ("Thickness (m)", :mean_thickness),
        ("Wet mass (mmWE)", :mean_wet_mass),
        ("Bulk density (kg m-3)", :mean_bulk_density),
        ("Firn-to-ice mass (mmWE)", :mean_base_mass),
    )),
    (title="Last cycle deltas", fields=(
        ("Mean signed dThickness (m)", :mean_signed_delta_thickness),
        ("Mean abs dThickness (m)", :mean_abs_delta_thickness),
        ("Max abs dThickness (m)", :max_abs_delta_thickness),
        ("Mean signed dSWE (mmWE)", :mean_signed_delta_wet_mass),
        ("Mean abs dSWE (mmWE)", :mean_abs_delta_wet_mass),
        ("Max abs dSWE (mmWE)", :max_abs_delta_wet_mass),
        ("Mean signed dBase (mmWE)", :mean_signed_delta_base_mass),
        ("Mean abs dBase (mmWE)", :mean_abs_delta_base_mass),
        ("Max abs dBase (mmWE)", :max_abs_delta_base_mass),
    )),
)

function write_case_summary(
    out_path::AbstractString,
    options::RunConfig,
    time_values::Vector{DateTime},
    ncol::Int,
    history::Vector{NamedTuple},
    status::Symbol,
    timings::StepTimingStats,
)
    last_record = history[end]
    mkpath(dirname(out_path))
    open(out_path, "w") do io
        println(io, options.name)
        println(io, "Input label        : ", isempty(options.input_label) ? "(not provided)" : options.input_label)
        println(io, "Forcing start      : ", first(time_values))
        println(io, "Forcing end        : ", last(time_values))
        println(io, "Forcing steps      : ", length(time_values))
        println(io, "Columns            : ", ncol)
        println(io, "Backend            : ", String(options.backend))
        println(io, "Threads            : ", nthreads())
        println(io, "File output        : ", options.write_outputs ? "enabled" : "disabled (--no-output)")
        println(io, "NetCDF output      : ", options.write_netcdf ? "enabled" : "disabled (--no-nc)")
        println(io, "Cycle metrics      : ", cycle_metrics_schedule_label(options.history_stride))
        println(io, "Status             : ", string(status))
        println(io, "Cycles completed   : ", completed_cycle_count(history, status, options.cycles))
        for section in SUMMARY_REPORT_SPECS
            println(io)
            println(io, section.title)
            for (label, key) in section.fields
                println(io, @sprintf("%-28s : %.6f", label, getfield(last_record, key)))
            end
        end
        println(io)
        println(io, "Interpretation     : Requested cycles completed.")
        println(io)
        print_timing_summary(io, timings)
    end
end

function collect_final_layer_grids(
    domain::SnowpackDomain,
    js::Vector{Int},
    is::Vector{Int},
    grid_shape::Tuple{Int, Int},
    nlayer::Int,
)
    ny, nx = grid_shape
    n_active = fill(Int32(-1), ny, nx)
    layer_density = fill(NaN, nlayer, ny, nx)
    layer_thickness = fill(NaN, nlayer, ny, nx)
    layer_snow_mass = fill(NaN, nlayer, ny, nx)
    layer_liquid_mass = fill(NaN, nlayer, ny, nx)
    layer_temperature_c = fill(NaN, nlayer, ny, nx)
    @inbounds for idx in 1:column_count(domain)
        j, i = js[idx], is[idx]
        n_active[j, i] = Int32(domain.N[idx])
        for k in 1:domain.N[idx]
            rho = domain.density[k, idx]
            m = domain.mass[k, idx]
            mw = domain.mass_w[k, idx]
            layer_density[k, j, i] = rho
            layer_snow_mass[k, j, i] = m
            layer_liquid_mass[k, j, i] = mw
            layer_temperature_c[k, j, i] = domain.temperature[k, idx] - domain.c.T0
            if isfinite(rho) && rho > 0.0 && isfinite(m)
                layer_thickness[k, j, i] = m / rho
            end
        end
    end
    return (
        n_active=n_active,
        layer_density=layer_density,
        layer_thickness=layer_thickness,
        layer_snow_mass=layer_snow_mass,
        layer_liquid_mass=layer_liquid_mass,
        layer_temperature_c=layer_temperature_c,
    )
end

empty_final_grids() = NamedTuple{FINAL_GRID_KEYS}(ntuple(_ -> Matrix{Float64}(undef, 0, 0), length(FINAL_GRID_KEYS)))
empty_monthly_grids() = NamedTuple{MONTHLY_GRID_KEYS}(ntuple(_ -> Array{Float64}(undef, 0, 0, 0), length(MONTHLY_GRID_KEYS)))

function _empty_layer_grids()
    return (
        n_active=Matrix{Int32}(undef, 0, 0),
        layer_density=Array{Float64}(undef, 0, 0, 0),
        layer_thickness=Array{Float64}(undef, 0, 0, 0),
        layer_snow_mass=Array{Float64}(undef, 0, 0, 0),
        layer_liquid_mass=Array{Float64}(undef, 0, 0, 0),
        layer_temperature_c=Array{Float64}(undef, 0, 0, 0),
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
