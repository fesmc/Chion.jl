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

@inline _grid_shape(layout::GridLayout) = size(layout.mask)

const SUMMARY_BUFFER_NAMES = (:thickness, :wet_mass, :bulk_density, :base_mass, :smb_ice, :liquid_water, :runoff)
const CYCLE_BUFFER_NAMES = (:thickness, :wet_mass, :bulk_density, :base_mass)

_named_buffers(names::NTuple{N, Symbol}, build::F) where {N, F <: Function} = NamedTuple{names}(ntuple(_ -> build(), N))

allocate_summary_buffers(n::Int) = _named_buffers(SUMMARY_BUFFER_NAMES, () -> Vector{Float64}(undef, n))
allocate_summary_buffers(domain::SM.AbstractSnowpackDomain, n::Int) = _named_buffers(SUMMARY_BUFFER_NAMES, () -> similar(domain.mass, Float64, n))
allocate_cycle_summary_buffers(n::Int) = _named_buffers(CYCLE_BUFFER_NAMES, () -> Vector{Float64}(undef, n))
allocate_cycle_summary_buffers(domain::SM.AbstractSnowpackDomain, n::Int) = _named_buffers(CYCLE_BUFFER_NAMES, () -> similar(domain.mass, Float64, n))

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
    timings::TimingStats,
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
    domain::SM.SnowpackDomain,
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
    @inbounds for idx in 1:SM.column_count(domain)
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

    history = NamedTuple[]
    simulation_wall_t0 = time_ns()
    for cycle in 1:options.cycles
        time_counted_block!(timings, :model_step_wall, ncol * size(step_fields.air_temperature, 2); synchronize=backend.sync) do
            _step_cycle!(domain, step_fields, workspaces, options)
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

@inline _nc_float_attrib(long_name::AbstractString, units::AbstractString) = Dict("long_name" => String(long_name), "units" => String(units))
@inline _nc_text_attrib(long_name::AbstractString) = Dict("long_name" => String(long_name))
_define_nc_output_variable(ds::NCDataset, name::AbstractString, dims; long_name::AbstractString, units::AbstractString) =
    defVar(ds, String(name), Float32, dims; fillvalue=NaN32, attrib=_nc_float_attrib(long_name, units))
_define_nc_int_variable(ds::NCDataset, name::AbstractString, dims; long_name::AbstractString) =
    defVar(ds, String(name), Int32, dims; attrib=_nc_text_attrib(long_name))
_define_nc_double_variable(ds::NCDataset, name::AbstractString, dims; attrib::Dict{String, String}=Dict{String, String}()) =
    defVar(ds, String(name), Float64, dims; attrib=attrib)

const CASE_NC_SPECS = (
    (key=:final_thickness, dims=("x", "y"), name="final_thickness", long_name="Final snow thickness", units="m", integer=false),
    (key=:final_wet_mass, dims=("x", "y"), name="final_wet_mass", long_name="Final snow wet mass", units="mmWE", integer=false),
    (key=:final_bulk_density, dims=("x", "y"), name="final_bulk_density", long_name="Final bulk snow density", units="kg m-3", integer=false),
    (key=:final_base_mass, dims=("x", "y"), name="final_base_mass", long_name="Cumulative firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:final_ice_sheet_smb, dims=("x", "y"), name="final_ice_sheet_smb", long_name="Cumulative net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:final_runoff, dims=("x", "y"), name="final_runoff", long_name="Final cumulative runoff", units="mmWE", integer=false),
    (key=:last_cycle_delta_thickness, dims=("x", "y"), name="last_cycle_delta_thickness", long_name="Last cycle snow-thickness change", units="m", integer=false),
    (key=:last_cycle_delta_wet_mass, dims=("x", "y"), name="last_cycle_delta_wet_mass", long_name="Last cycle wet-mass change", units="mmWE", integer=false),
    (key=:last_cycle_delta_base_mass, dims=("x", "y"), name="last_cycle_delta_base_mass", long_name="Last cycle firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:last_cycle_delta_ice_sheet_smb, dims=("x", "y"), name="last_cycle_delta_ice_sheet_smb", long_name="Last cycle net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:n_active, dims=("x", "y"), name="n_active", long_name="Number of active Chion layers", units="", integer=true),
    (key=:layer_density, dims=("layer", "x", "y"), name="layer_density", long_name="Final Chion layer density", units="kg m-3", integer=false),
    (key=:layer_thickness, dims=("layer", "x", "y"), name="layer_thickness", long_name="Final Chion layer thickness", units="m", integer=false),
    (key=:layer_snow_mass, dims=("layer", "x", "y"), name="layer_snow_mass", long_name="Final Chion layer snow mass", units="kg m-2", integer=false),
    (key=:layer_liquid_mass, dims=("layer", "x", "y"), name="layer_liquid_mass", long_name="Final Chion layer liquid-water mass", units="kg m-2", integer=false),
    (key=:layer_temperature_c, dims=("layer", "x", "y"), name="layer_temperature_c", long_name="Final Chion layer temperature", units="C", integer=false),
    (key=:history_mean_thickness, dims=("cycle",), name="history_mean_thickness", long_name="Cycle-mean snow thickness", units="m", integer=false),
    (key=:history_mean_wet_mass, dims=("cycle",), name="history_mean_wet_mass", long_name="Cycle-mean snow wet mass", units="mmWE", integer=false),
    (key=:history_mean_bulk_density, dims=("cycle",), name="history_mean_bulk_density", long_name="Cycle-mean bulk snow density", units="kg m-3", integer=false),
    (key=:history_mean_base_mass, dims=("cycle",), name="history_mean_base_mass", long_name="Cycle-mean firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:history_mean_abs_delta_thickness, dims=("cycle",), name="history_mean_abs_delta_thickness", long_name="Cycle mean absolute snow-thickness change", units="m", integer=false),
    (key=:history_mean_abs_delta_wet_mass, dims=("cycle",), name="history_mean_abs_delta_wet_mass", long_name="Cycle mean absolute wet-mass change", units="mmWE", integer=false),
    (key=:history_mean_abs_delta_base_mass, dims=("cycle",), name="history_mean_abs_delta_base_mass", long_name="Cycle mean absolute firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:monthly_mean_thickness, dims=("month", "x", "y"), name="monthly_mean_thickness", long_name="Monthly mean snow thickness", units="m", integer=false),
    (key=:monthly_mean_wet_mass, dims=("month", "x", "y"), name="monthly_mean_wet_mass", long_name="Monthly mean snow wet mass", units="mmWE", integer=false),
    (key=:monthly_mean_bulk_density, dims=("month", "x", "y"), name="monthly_mean_bulk_density", long_name="Monthly mean bulk snow density", units="kg m-3", integer=false),
    (key=:monthly_mean_base_mass, dims=("month", "x", "y"), name="monthly_mean_base_mass", long_name="Monthly mean cumulative firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:monthly_mean_ice_sheet_smb, dims=("month", "x", "y"), name="monthly_mean_ice_sheet_smb", long_name="Monthly net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:monthly_export_to_ice, dims=("month", "x", "y"), name="monthly_export_to_ice", long_name="Monthly firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:monthly_net_ice_sheet_forcing, dims=("month", "x", "y"), name="monthly_net_ice_sheet_forcing", long_name="Monthly net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:monthly_runoff, dims=("month", "x", "y"), name="monthly_runoff", long_name="Monthly runoff production", units="mmWE", integer=false),
    (key=:step_export_to_ice, dims=("step", "x", "y"), name="step_export_to_ice", long_name="Annual firn mass exported to the ice model for each written output interval", units="mmWE", integer=false),
    (key=:step_ice_sheet_smb, dims=("step", "x", "y"), name="step_ice_sheet_smb", long_name="Annual net mass forcing to the ice sheet for each written output interval", units="mmWE", integer=false),
)

function _define_selected_nc_variables!(ds::NCDataset, selected::Set{Symbol})
    vars = Dict{Symbol, Any}()
    for spec in CASE_NC_SPECS
        spec.key in selected || continue
        vars[spec.key] = spec.integer ?
            _define_nc_int_variable(ds, spec.name, spec.dims; long_name=spec.long_name) :
            _define_nc_output_variable(ds, spec.name, spec.dims; long_name=spec.long_name, units=spec.units)
    end
    return vars
end

function init_case_netcdf(
    netcdf_path::AbstractString,
    options::RunConfig,
    time_values::Vector{DateTime},
    nlayer::Int,
    layout::GridLayout,
    initial_thickness::Matrix{Float64},
    month_cycle::Vector{Int32},
    month_of_year::Vector{Int32},
    source_month_code::Vector{Int32},
    annual_output_source_indices::Vector{Int32},
    annual_output_source_codes::Vector{Int32},
)
    mkpath(dirname(netcdf_path))
    isdir(netcdf_path) && error("NetCDF output path '$(abspath(netcdf_path))' is a directory; pass a file path ending in `.nc`.")
    ny, nx = _grid_shape(layout)
    max_steps = options.cycles * length(annual_output_source_indices)
    selected = Set(options.netcdf_variables)
    step_cycle = Int32[cyc for cyc in 1:options.cycles for _ in annual_output_source_indices]
    step_source_index = Int32[idx for _ in 1:options.cycles for idx in annual_output_source_indices]
    step_source_code = Int32[code for _ in 1:options.cycles for code in annual_output_source_codes]

    ds = NCDataset(netcdf_path, "c")
    for (name, len) in (("x", nx), ("y", ny), ("layer", max(nlayer, 1)), ("cycle", options.cycles), ("month", length(month_cycle)), ("point", length(layout.js)), ("step", max_steps))
        defDim(ds, name, len)
    end

    _define_nc_double_variable(ds, "x", ("x",); attrib=Dict("units" => "km", "axis" => "X"))[:] = layout.x
    _define_nc_double_variable(ds, "y", ("y",); attrib=Dict("units" => "km", "axis" => "Y"))[:] = layout.y
    _define_nc_int_variable(ds, "layer", ("layer",); long_name="Chion internal layer index from surface downward")[:] = Int32.(collect(1:max(nlayer, 1)))
    _define_nc_int_variable(ds, "cycle", ("cycle",); long_name="Repeated annual forcing cycle index")[:] = Int32.(collect(1:options.cycles))
    _define_nc_int_variable(ds, "month", ("month",); long_name="Sequential monthly output index")[:] = Int32.(collect(1:length(month_cycle)))
    _define_nc_int_variable(ds, "month_cycle", ("month",); long_name="Forcing cycle associated with monthly output")[:] = month_cycle
    _define_nc_int_variable(ds, "month_of_year", ("month",); long_name="Calendar month of the repeated forcing")[:] = month_of_year
    _define_nc_int_variable(ds, "source_month_code", ("month",); long_name="Source forcing month code YYYYMM")[:] = source_month_code
    _define_nc_int_variable(ds, "point", ("point",); long_name="Compact valid cell index")[:] = Int32.(collect(1:length(layout.js)))
    _define_nc_int_variable(ds, "point_j", ("point",); long_name="1-based y-index for each compact valid cell")[:] = Int32.(layout.js)
    _define_nc_int_variable(ds, "point_i", ("point",); long_name="1-based x-index for each compact valid cell")[:] = Int32.(layout.is)
    _define_nc_double_variable(ds, "point_y_km", ("point",); attrib=Dict("long_name" => "Y coordinate for each compact valid cell", "units" => "km"))[:] = layout.y[layout.js]
    _define_nc_double_variable(ds, "point_x_km", ("point",); attrib=Dict("long_name" => "X coordinate for each compact valid cell", "units" => "km"))[:] = layout.x[layout.is]
    _define_nc_int_variable(ds, "step", ("step",); long_name="Sequential yearly output index across repeated annual cycles")[:] = Int32.(collect(1:max_steps))
    _define_nc_int_variable(ds, "step_cycle", ("step",); long_name="Repeated annual forcing cycle index for each yearly output")[:] = step_cycle
    _define_nc_int_variable(ds, "step_source_index", ("step",); long_name="1-based index of the last forcing step included in each yearly output")[:] = step_source_index
    _define_nc_int_variable(ds, "step_source_code", ("step",); long_name="Source forcing timestamp code YYYYMMDDHH for the final step included in each yearly output")[:] = step_source_code
    _define_nc_output_variable(ds, "domain_mask", ("x", "y"); long_name="Domain mask", units="1")[:, :] = Float32.(permutedims(layout.mask, (2, 1)))
    _define_nc_output_variable(ds, "initial_thickness", ("x", "y"); long_name="Initial snow thickness", units="m")[:, :] = Float32.(permutedims(initial_thickness, (2, 1)))

    vars = _define_selected_nc_variables!(ds, selected)
    if any(key -> key in selected, CASE_OUTPUT_GROUPS.step)
        vars[:step_valid] = _define_nc_int_variable(ds, "step_valid", ("step",); long_name="1 where a yearly output record was completed and written, 0 for unused trailing slots")
    end

    ds.attrib["title"] = options.name
    ds.attrib["source_model"] = "Chion"
    ds.attrib["input_label"] = isempty(options.input_label) ? "not provided" : options.input_label
    ds.attrib["forcing_start"] = string(first(time_values))
    ds.attrib["forcing_end"] = string(last(time_values))
    ds.attrib["cycles_completed"] = "pending"
    ds.attrib["status"] = "pending"
    ds.attrib["created"] = string(now())
    return CaseNetCDFWriter(ds, vars, max_steps, options.cycles)
end

function maybe_write_step_output!(writer::CaseNetCDFWriter, step_index::Int, key::Symbol, data::AbstractMatrix{<:Real})
    haskey(writer.vars, key) && (writer.vars[key][step_index, :, :] = Float32.(permutedims(data, (2, 1))))
    return
end

function _write_dataset_var!(var, data::Vector{Float64})
    var[:] = Float32.(data)
    return
end

function _write_dataset_var!(var, data::Vector{Int32})
    var[:] = data
    return
end

function _write_dataset_var!(var, data::AbstractMatrix{Int32})
    var[:, :] = permutedims(data, (2, 1))
    return
end

function _write_dataset_var!(var, data::AbstractMatrix{<:Real})
    var[:, :] = Float32.(permutedims(data, (2, 1)))
    return
end

function _write_dataset_var!(var, data::Array{Float64, 3})
    var[:, :, :] = Float32.(permutedims(data, (1, 3, 2)))
    return
end

function maybe_write_output!(writer::CaseNetCDFWriter, key::Symbol, data)
    haskey(writer.vars, key) && _write_dataset_var!(writer.vars[key], data)
    return
end

function _history_vectors(history::Vector{NamedTuple}, cycles::Int)
    out = Dict(spec.output => fill(NaN, cycles) for spec in HISTORY_OUTPUT_SPECS)
    for rec in history
        idx = rec.cycle
        for spec in HISTORY_OUTPUT_SPECS
            out[spec.output][idx] = getfield(rec, spec.record)
        end
    end
    return out
end

function finalize_case_netcdf!(
    writer::CaseNetCDFWriter,
    final_grids,
    layer_grids,
    history::Vector{NamedTuple},
    monthly_grids,
    status::Symbol,
    cycles_completed::Int,
    steps_written::Int,
)
    for key in FINAL_GRID_KEYS
        maybe_write_output!(writer, key, getfield(final_grids, key))
    end
    for key in LAYER_GRID_KEYS
        maybe_write_output!(writer, key, getfield(layer_grids, key))
    end
    for (key, values) in _history_vectors(history, writer.cycles)
        maybe_write_output!(writer, key, values)
    end
    for key in MONTHLY_GRID_KEYS
        maybe_write_output!(writer, key, getfield(monthly_grids, key))
    end
    if haskey(writer.vars, :step_valid)
        step_valid = zeros(Int32, writer.max_steps)
        step_valid[1:steps_written] .= 1
        writer.vars[:step_valid][:] = step_valid
    end
    writer.dataset.attrib["cycles_completed"] = string(cycles_completed)
    writer.dataset.attrib["status"] = string(status)
    writer.dataset.attrib["steps_written"] = string(steps_written)
    close(writer.dataset)
    return
end
