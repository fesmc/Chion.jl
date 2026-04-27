# Output grids, NetCDF schema, and NetCDF writer helpers for Simulation runs.

const FINAL_OUTPUT_SPECS = (
    (key=:final_thickness, source=:final_state, field=:thickness),
    (key=:final_wet_mass, source=:final_state, field=:wet_mass),
    (key=:final_bulk_density, source=:final_state, field=:bulk_density),
    (key=:final_base_mass, source=:final_state, field=:base_mass),
    (key=:final_ice_sheet_smb, source=:domain, field=:smb_ice),
    (key=:final_runoff, source=:domain, field=:runoff),
    (key=:last_cycle_delta_thickness, source=:deltas, field=:thickness),
    (key=:last_cycle_delta_wet_mass, source=:deltas, field=:wet_mass),
    (key=:last_cycle_delta_base_mass, source=:deltas, field=:base_mass),
    (key=:last_cycle_delta_ice_sheet_smb, source=:deltas, field=:ice_sheet_smb),
)

const HISTORY_OUTPUT_SPECS = (
    (key=:history_mean_thickness, record=:mean_thickness),
    (key=:history_mean_wet_mass, record=:mean_wet_mass),
    (key=:history_mean_bulk_density, record=:mean_bulk_density),
    (key=:history_mean_base_mass, record=:mean_base_mass),
    (key=:history_mean_abs_delta_thickness, record=:mean_abs_delta_thickness),
    (key=:history_mean_abs_delta_wet_mass, record=:mean_abs_delta_wet_mass),
    (key=:history_mean_abs_delta_base_mass, record=:mean_abs_delta_base_mass),
)

const MONTHLY_OUTPUT_SPECS = (
    (key=:monthly_mean_thickness, aggregate=:mean),
    (key=:monthly_mean_wet_mass, aggregate=:mean),
    (key=:monthly_mean_bulk_density, aggregate=:mean),
    (key=:monthly_mean_base_mass, aggregate=:sum),
    (key=:monthly_mean_ice_sheet_smb, aggregate=:sum),
    (key=:monthly_export_to_ice, aggregate=:sum),
    (key=:monthly_net_ice_sheet_forcing, aggregate=:sum),
    (key=:monthly_runoff, aggregate=:sum),
)

const OUTPUT_GROUPS = (
    final=map(spec -> spec.key, FINAL_OUTPUT_SPECS),
    layers=(
        :n_active,
        :layer_density,
        :layer_thickness,
        :layer_snow_mass,
        :layer_liquid_mass,
        :layer_temperature_c,
    ),
    history=map(spec -> spec.key, HISTORY_OUTPUT_SPECS),
    monthly=map(spec -> spec.key, MONTHLY_OUTPUT_SPECS),
    step=(:step_export_to_ice, :step_ice_sheet_smb, :step_pdd),
)
const NETCDF_VARIABLES = unique(Symbol[var for group in values(OUTPUT_GROUPS) for var in group])
const FINAL_GRID_KEYS = OUTPUT_GROUPS.final
const LAYER_GRID_KEYS = OUTPUT_GROUPS.layers
const MONTHLY_GRID_KEYS = OUTPUT_GROUPS.monthly
@inline _grid_shape(layout) = size(layout.mask)

@inline _host_vector(data::Vector{Float64}; copy_array::Bool=false) = copy_array ? copy(data) : data
@inline _host_vector(data; copy_array::Bool=false) = Float64.(Array(data))

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

empty_final_grids() = NamedTuple{FINAL_GRID_KEYS}(ntuple(_ -> Matrix{Float64}(undef, 0, 0), length(FINAL_GRID_KEYS)))
empty_monthly_grids() = NamedTuple{MONTHLY_GRID_KEYS}(ntuple(_ -> Array{Float64}(undef, 0, 0, 0), length(MONTHLY_GRID_KEYS)))

function collect_final_layer_grids(
    domain::SnowpackDomain,
    js::Vector{Int},
    is::Vector{Int},
    grid_shape::Tuple{Int, Int},
    nlayer::Int,
)
    ny, nx = grid_shape
    ncol = length(js)
    n_active = fill(Int32(0), ny, nx)
    layer_density     = fill(NaN, nlayer, ny, nx)
    layer_thickness   = fill(NaN, nlayer, ny, nx)
    layer_snow_mass   = fill(NaN, nlayer, ny, nx)
    layer_liquid_mass = fill(NaN, nlayer, ny, nx)
    layer_temperature_c = fill(NaN, nlayer, ny, nx)
    c = domain.c
    @inbounds for col in 1:ncol
        j, i = js[col], is[col]
        n_active[j, i] = Int32(domain.N[col])
        for k in 1:nlayer
            rho = domain.density[k, col]
            m   = domain.mass[k, col]
            layer_density[k, j, i]       = rho
            layer_snow_mass[k, j, i]     = m
            layer_liquid_mass[k, j, i]   = domain.mass_w[k, col]
            layer_temperature_c[k, j, i] = domain.temperature[k, col] - c.T0
            layer_thickness[k, j, i]     = rho > 0 ? m / rho : 0.0
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

_allocate_step_vectors(active::Bool, ncol::Int) = NamedTuple{OUTPUT_GROUPS.step}(ntuple(_ -> active ? zeros(Float64, ncol) : Float64[], length(OUTPUT_GROUPS.step)))
_allocate_monthly_sums(active::Bool, nmonth_total::Int, ncol::Int) = NamedTuple{MONTHLY_GRID_KEYS}(ntuple(_ -> active ? zeros(Float64, nmonth_total, ncol) : Matrix{Float64}(undef, 0, 0), length(MONTHLY_GRID_KEYS)))

function _step_output_grids(step_vectors, layout)
    return NamedTuple{OUTPUT_GROUPS.step}(ntuple(i -> scatter_to_grid(getfield(step_vectors, OUTPUT_GROUPS.step[i]), layout.js, layout.is, _grid_shape(layout)), length(OUTPUT_GROUPS.step)))
end

function _reset_step_vectors!(step_vectors)
    for key in OUTPUT_GROUPS.step
        isempty(getfield(step_vectors, key)) || fill!(getfield(step_vectors, key), 0.0)
    end
    return
end

function _accumulate_step_diagnostics!(
    summary,
    previous,
    monthly_sums,
    step_vectors,
    month_idx::Int,
    need_monthly_outputs::Bool,
    need_step_outputs::Bool,
)
    current_base = summary.base_mass
    current_smb = summary.smb_ice
    current_runoff = summary.runoff
    current_pdd = summary.pdd
    delta_base = current_base .- previous.base_mass
    delta_smb = current_smb .- previous.smb_ice
    delta_runoff = current_runoff .- previous.runoff
    delta_pdd = current_pdd .- previous.pdd

    if need_monthly_outputs
        monthly_sums.monthly_mean_thickness[month_idx, :] .+= summary.thickness
        monthly_sums.monthly_mean_wet_mass[month_idx, :] .+= summary.wet_mass
        monthly_sums.monthly_mean_bulk_density[month_idx, :] .+= summary.bulk_density
        monthly_sums.monthly_mean_base_mass[month_idx, :] .+= current_base
        monthly_sums.monthly_mean_ice_sheet_smb[month_idx, :] .+= delta_smb
        monthly_sums.monthly_export_to_ice[month_idx, :] .+= delta_base
        monthly_sums.monthly_net_ice_sheet_forcing[month_idx, :] .+= delta_smb
        monthly_sums.monthly_runoff[month_idx, :] .+= delta_runoff
    end
    if need_step_outputs
        step_vectors.step_export_to_ice .+= delta_base
        step_vectors.step_ice_sheet_smb .+= delta_smb
        step_vectors.step_pdd .+= delta_pdd
    end
    previous.base_mass .= current_base
    previous.smb_ice .= current_smb
    previous.runoff .= current_runoff
    previous.pdd .= current_pdd
    return
end

function _update_cycle_smb_delta!(last_delta::Vector{Float64}, previous_cycle_smb_ice::Vector{Float64}, domain)
    current = _host_vector(domain.smb_ice; copy_array=true)
    last_delta .= current .- previous_cycle_smb_ice
    previous_cycle_smb_ice .= current
    return
end

function _final_output_vector(spec, final_state, domain, deltas)
    source = spec.source
    values = source === :final_state ? getfield(final_state, spec.field) :
        source === :domain ? getfield(domain, spec.field) :
        source === :deltas ? getfield(deltas, spec.field) :
        error("Unsupported final output source `$(source)`.")
    return _host_vector(values; copy_array=true)
end

function _scatter_final_grids(final_state, domain, deltas, layout)
    return NamedTuple{FINAL_GRID_KEYS}(ntuple(i -> begin
        values = _final_output_vector(FINAL_OUTPUT_SPECS[i], final_state, domain, deltas)
        scatter_to_grid(values, layout.js, layout.is, _grid_shape(layout))
    end, length(FINAL_OUTPUT_SPECS)))
end

function _finalize_monthly_grids(monthly_sums, monthly_count::Vector{Int32}, layout)
    vectors = NamedTuple{MONTHLY_GRID_KEYS}(ntuple(i -> begin
        spec = MONTHLY_OUTPUT_SPECS[i]
        key = spec.key
        data = copy(getfield(monthly_sums, key))
        if spec.aggregate === :mean
            @inbounds for m in axes(data, 1)
                data[m, :] ./= max(monthly_count[m], 1)
            end
        end
        data
    end, length(MONTHLY_GRID_KEYS)))
    return NamedTuple{MONTHLY_GRID_KEYS}(ntuple(i -> monthly_vectors_to_grids(getfield(vectors, MONTHLY_GRID_KEYS[i]), layout.js, layout.is, _grid_shape(layout)), length(MONTHLY_GRID_KEYS)))
end
struct NetCDFWriter
    dataset::NCDataset
    vars::Dict{Symbol, Any}
    max_steps::Int
    cycles::Int
end

@inline _looks_like_directory_path(path::AbstractString) =
    !isempty(path) && (endswith(path, '/') || endswith(path, '\\'))

function _slug(name::AbstractString)
    slug = strip(replace(lowercase(strip(String(name))), r"[^a-z0-9]+" => "_"), '_')
    return isempty(slug) ? "run" : slug
end

_default_output_dir(name::AbstractString) = joinpath(pwd(), "run_output", _slug(name))

function resolve_netcdf_path(options)
    default_name = "$(options.name)_final_state.nc"
    isempty(options.netcdf_path) && return joinpath(options.output_dir, default_name)
    return isdir(options.netcdf_path) || _looks_like_directory_path(options.netcdf_path) ?
        joinpath(options.netcdf_path, default_name) :
        options.netcdf_path
end

function normalize_netcdf_variables(spec)
    if spec isa AbstractVector
        tokens = String[string(x) for x in spec]
    else
        text = lowercase(strip(String(spec)))
        isempty(text) && return copy(NETCDF_VARIABLES)
        tokens = split(text, ',')
    end
    selected = Symbol[]
    allowed_groups = String.(propertynames(OUTPUT_GROUPS))
    for token in tokens
        stripped = strip(token)
        isempty(stripped) && continue
        key = Symbol(lowercase(stripped))
        if key == :all
            append!(selected, NETCDF_VARIABLES)
        elseif key == :none
            continue
        elseif hasproperty(OUTPUT_GROUPS, key)
            append!(selected, getproperty(OUTPUT_GROUPS, key))
        elseif key in NETCDF_VARIABLES
            push!(selected, key)
        else
            error("Unsupported NetCDF variable selector '$token'. Use `all`, `none`, a group ($(join(sort!(allowed_groups), ", "))), or an explicit variable name.")
        end
    end
    return unique(selected)
end

@inline _nc_attrib(long_name::AbstractString, units::AbstractString="") =
    units == "" ? Dict("long_name" => String(long_name)) : Dict("long_name" => String(long_name), "units" => String(units))

@inline function _def_nc_var(ds::NCDataset, name::AbstractString, dims, spec)
    return spec.integer ?
        defVar(ds, String(name), Int32, dims; attrib=_nc_attrib(spec.long_name)) :
        defVar(ds, String(name), Float32, dims; fillvalue=NaN32, attrib=_nc_attrib(spec.long_name, spec.units))
end

@inline _write_nc_var!(ds::NCDataset, name, dims, spec, value) = (_def_nc_var(ds, name, dims, spec)[:] = value)

const NC_SPECS = (
    (key=:final_thickness,        name="final_thickness",        dims=("x", "y"),       long_name="Final snow thickness", units="m", integer=false),
    (key=:final_wet_mass,         name="final_wet_mass",         dims=("x", "y"),       long_name="Final snow wet mass", units="mmWE", integer=false),
    (key=:final_bulk_density,     name="final_bulk_density",     dims=("x", "y"),       long_name="Final bulk snow density", units="kg m-3", integer=false),
    (key=:final_base_mass,        name="final_base_mass",        dims=("x", "y"),       long_name="Cumulative firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:final_ice_sheet_smb,    name="final_ice_sheet_smb",    dims=("x", "y"),       long_name="Cumulative net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:final_runoff,           name="final_runoff",           dims=("x", "y"),       long_name="Final cumulative runoff", units="mmWE", integer=false),
    (key=:last_cycle_delta_thickness, name="last_cycle_delta_thickness", dims=("x", "y"), long_name="Last cycle snow-thickness change", units="m", integer=false),
    (key=:last_cycle_delta_wet_mass,  name="last_cycle_delta_wet_mass",  dims=("x", "y"), long_name="Last cycle wet-mass change", units="mmWE", integer=false),
    (key=:last_cycle_delta_base_mass, name="last_cycle_delta_base_mass", dims=("x", "y"), long_name="Last cycle firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:last_cycle_delta_ice_sheet_smb, name="last_cycle_delta_ice_sheet_smb", dims=("x", "y"), long_name="Last cycle net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:n_active,               name="n_active",               dims=("x", "y"),       long_name="Number of active Chion layers", units="", integer=true),
    (key=:layer_density,         name="layer_density",         dims=("layer", "x", "y"), long_name="Final Chion layer density", units="kg m-3", integer=false),
    (key=:layer_thickness,       name="layer_thickness",       dims=("layer", "x", "y"), long_name="Final Chion layer thickness", units="m", integer=false),
    (key=:layer_snow_mass,       name="layer_snow_mass",       dims=("layer", "x", "y"), long_name="Final Chion layer snow mass", units="kg m-2", integer=false),
    (key=:layer_liquid_mass,     name="layer_liquid_mass",     dims=("layer", "x", "y"), long_name="Final Chion layer liquid-water mass", units="kg m-2", integer=false),
    (key=:layer_temperature_c,   name="layer_temperature_c",   dims=("layer", "x", "y"), long_name="Final Chion layer temperature", units="C", integer=false),
    (key=:history_mean_thickness,      name="history_mean_thickness",      dims=("cycle",), long_name="Cycle-mean snow thickness", units="m", integer=false),
    (key=:history_mean_wet_mass,       name="history_mean_wet_mass",       dims=("cycle",), long_name="Cycle-mean snow wet mass", units="mmWE", integer=false),
    (key=:history_mean_bulk_density,   name="history_mean_bulk_density",   dims=("cycle",), long_name="Cycle-mean bulk snow density", units="kg m-3", integer=false),
    (key=:history_mean_base_mass,      name="history_mean_base_mass",      dims=("cycle",), long_name="Cycle-mean firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:history_mean_abs_delta_thickness, name="history_mean_abs_delta_thickness", dims=("cycle",), long_name="Cycle mean absolute snow-thickness change", units="m", integer=false),
    (key=:history_mean_abs_delta_wet_mass,  name="history_mean_abs_delta_wet_mass",  dims=("cycle",), long_name="Cycle mean absolute wet-mass change", units="mmWE", integer=false),
    (key=:history_mean_abs_delta_base_mass, name="history_mean_abs_delta_base_mass", dims=("cycle",), long_name="Cycle mean absolute firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:monthly_mean_thickness,      name="monthly_mean_thickness",      dims=("month", "x", "y"), long_name="Monthly mean snow thickness", units="m", integer=false),
    (key=:monthly_mean_wet_mass,       name="monthly_mean_wet_mass",       dims=("month", "x", "y"), long_name="Monthly mean snow wet mass", units="mmWE", integer=false),
    (key=:monthly_mean_bulk_density,   name="monthly_mean_bulk_density",   dims=("month", "x", "y"), long_name="Monthly mean bulk snow density", units="kg m-3", integer=false),
    (key=:monthly_mean_base_mass,      name="monthly_mean_base_mass",      dims=("month", "x", "y"), long_name="Monthly mean cumulative firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:monthly_mean_ice_sheet_smb,  name="monthly_mean_ice_sheet_smb",  dims=("month", "x", "y"), long_name="Monthly net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:monthly_export_to_ice,       name="monthly_export_to_ice",       dims=("month", "x", "y"), long_name="Monthly firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:monthly_net_ice_sheet_forcing, name="monthly_net_ice_sheet_forcing", dims=("month", "x", "y"), long_name="Monthly net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:monthly_runoff,              name="monthly_runoff",              dims=("month", "x", "y"), long_name="Monthly runoff production", units="mmWE", integer=false),
    (key=:step_export_to_ice,          name="step_export_to_ice",          dims=("step", "x", "y"), long_name="Annual firn mass exported to the ice model for each written output interval", units="mmWE", integer=false),
    (key=:step_ice_sheet_smb,          name="step_ice_sheet_smb",          dims=("step", "x", "y"), long_name="Annual net mass forcing to the ice sheet for each written output interval", units="mmWE", integer=false),
    (key=:step_pdd,                    name="step_pdd",                    dims=("step", "x", "y"), long_name="Annual positive degree days for each written output interval", units="degC d", integer=false),
)

function _define_selected_nc_variables!(ds::NCDataset, selected::Set{Symbol})
    vars = Dict{Symbol, Any}()
    for spec in NC_SPECS
        spec.key in selected || continue
        vars[spec.key] = _def_nc_var(ds, spec.name, spec.dims, spec)
    end
    return vars
end

function init_netcdf(
    netcdf_path::AbstractString,
    options,
    time_values::Vector{DateTime},
    nlayer::Int,
    layout::SnowpackGrid,
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

    for spec in (
        (name="x", dims=("x",), meta=(key=:x, long_name="X coordinate", units="km", integer=false), value=layout.x),
        (name="y", dims=("y",), meta=(key=:y, long_name="Y coordinate", units="km", integer=false), value=layout.y),
        (name="layer", dims=("layer",), meta=(key=:layer, long_name="Chion internal layer index from surface downward", units="", integer=true), value=Int32.(collect(1:max(nlayer, 1)))),
        (name="cycle", dims=("cycle",), meta=(key=:cycle, long_name="Repeated annual forcing cycle index", units="", integer=true), value=Int32.(collect(1:options.cycles))),
        (name="month", dims=("month",), meta=(key=:month, long_name="Sequential monthly output index", units="", integer=true), value=Int32.(collect(1:length(month_cycle)))),
        (name="month_cycle", dims=("month",), meta=(key=:month_cycle, long_name="Forcing cycle associated with monthly output", units="", integer=true), value=month_cycle),
        (name="month_of_year", dims=("month",), meta=(key=:month_of_year, long_name="Calendar month of the repeated forcing", units="", integer=true), value=month_of_year),
        (name="source_month_code", dims=("month",), meta=(key=:source_month_code, long_name="Source forcing month code YYYYMM", units="", integer=true), value=source_month_code),
        (name="point", dims=("point",), meta=(key=:point, long_name="Compact valid cell index", units="", integer=true), value=Int32.(collect(1:length(layout.js)))),
        (name="point_j", dims=("point",), meta=(key=:point_j, long_name="1-based y-index for each compact valid cell", units="", integer=true), value=Int32.(layout.js)),
        (name="point_i", dims=("point",), meta=(key=:point_i, long_name="1-based x-index for each compact valid cell", units="", integer=true), value=Int32.(layout.is)),
        (name="point_y_km", dims=("point",), meta=(key=:point_y_km, long_name="Y coordinate for each compact valid cell", units="km", integer=false), value=layout.y[layout.js]),
        (name="point_x_km", dims=("point",), meta=(key=:point_x_km, long_name="X coordinate for each compact valid cell", units="km", integer=false), value=layout.x[layout.is]),
        (name="step", dims=("step",), meta=(key=:step, long_name="Sequential yearly output index across repeated annual cycles", units="", integer=true), value=Int32.(collect(1:max_steps))),
        (name="step_cycle", dims=("step",), meta=(key=:step_cycle, long_name="Repeated annual forcing cycle index for each yearly output", units="", integer=true), value=step_cycle),
        (name="step_source_index", dims=("step",), meta=(key=:step_source_index, long_name="1-based index of the last forcing step included in each yearly output", units="", integer=true), value=step_source_index),
        (name="step_source_code", dims=("step",), meta=(key=:step_source_code, long_name="Source forcing timestamp code YYYYMMDDHH for the final step included in each yearly output", units="", integer=true), value=step_source_code),
        (name="domain_mask", dims=("x", "y"), meta=(key=:domain_mask, long_name="Domain mask", units="1", integer=false), value=Float32.(permutedims(layout.mask, (2, 1)))),
        (name="initial_thickness", dims=("x", "y"), meta=(key=:initial_thickness, long_name="Initial snow thickness", units="m", integer=false), value=Float32.(permutedims(initial_thickness, (2, 1)))),
    )
        _write_nc_var!(ds, spec.name, spec.dims, spec.meta, spec.value)
    end

    vars = _define_selected_nc_variables!(ds, selected)
    any(key -> key in selected, OUTPUT_GROUPS.step) && (vars[:step_valid] = _def_nc_var(ds, "step_valid", ("step",), (key=:step_valid, long_name="1 where a yearly output record was completed and written, 0 for unused trailing slots", units="", integer=true)))

    ds.attrib["title"] = options.name
    ds.attrib["source_model"] = "Chion"
    ds.attrib["input_label"] = isempty(options.input_label) ? "not provided" : options.input_label
    ds.attrib["forcing_start"] = string(first(time_values))
    ds.attrib["forcing_end"] = string(last(time_values))
    ds.attrib["cycles_completed"] = "pending"
    ds.attrib["status"] = "pending"
    ds.attrib["created"] = string(now())
    return NetCDFWriter(ds, vars, max_steps, options.cycles)
end

@inline _write_dataset_var!(var, data::AbstractVector) = (var[:] = eltype(var) <: Integer ? data : Float32.(data))
@inline _write_dataset_var!(var, data::AbstractMatrix) = (var[:, :] = eltype(var) <: Integer ? permutedims(data, (2, 1)) : Float32.(permutedims(data, (2, 1))))
@inline _write_dataset_var!(var, data::AbstractArray{<:Real, 3}) = (var[:, :, :] = Float32.(permutedims(data, (1, 3, 2))))

maybe_write_step_output!(writer::NetCDFWriter, step_index::Int, key::Symbol, data::AbstractMatrix{<:Real}) =
    (haskey(writer.vars, key) && (writer.vars[key][step_index, :, :] = Float32.(permutedims(data, (2, 1))); nothing))

maybe_write_output!(writer::NetCDFWriter, key::Symbol, data) = (haskey(writer.vars, key) && (_write_dataset_var!(writer.vars[key], data); nothing))

function _history_vectors(history::Vector{NamedTuple}, cycles::Int)
    out = Dict(spec.key => fill(NaN, cycles) for spec in HISTORY_OUTPUT_SPECS)
    for rec in history, spec in HISTORY_OUTPUT_SPECS
        out[spec.key][rec.cycle] = getfield(rec, spec.record)
    end
    return out
end

function _write_output_group!(writer::NetCDFWriter, keys, values)
    for key in keys
        maybe_write_output!(writer, key, getfield(values, key))
    end
    return nothing
end

function finalize_netcdf!(
    writer::NetCDFWriter,
    final_grids,
    layer_grids,
    history::Vector{NamedTuple},
    monthly_grids,
    status::Symbol,
    cycles_completed::Int,
    steps_written::Int,
)
    _write_output_group!(writer, FINAL_GRID_KEYS, final_grids)
    _write_output_group!(writer, LAYER_GRID_KEYS, layer_grids)
    for (key, values) in _history_vectors(history, writer.cycles)
        maybe_write_output!(writer, key, values)
    end
    _write_output_group!(writer, MONTHLY_GRID_KEYS, monthly_grids)
    if haskey(writer.vars, :step_valid)
        step_valid = zeros(Int32, writer.max_steps)
        step_valid[1:steps_written] .= 1
        writer.vars[:step_valid][:] = step_valid
    end
    writer.dataset.attrib["cycles_completed"] = string(cycles_completed)
    writer.dataset.attrib["status"] = string(status)
    writer.dataset.attrib["steps_written"] = string(steps_written)
    close(writer.dataset)
end

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
