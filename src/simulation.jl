const OUTPUT_GROUPS = (
    final=(
        :final_thickness,
        :final_wet_mass,
        :final_bulk_density,
        :final_base_mass,
        :final_ice_sheet_smb,
        :final_runoff,
        :last_cycle_delta_thickness,
        :last_cycle_delta_wet_mass,
        :last_cycle_delta_base_mass,
        :last_cycle_delta_ice_sheet_smb,
    ),
    layers=(
        :n_active,
        :layer_density,
        :layer_thickness,
        :layer_snow_mass,
        :layer_liquid_mass,
        :layer_temperature_c,
    ),
    history=(
        :history_mean_thickness,
        :history_mean_wet_mass,
        :history_mean_bulk_density,
        :history_mean_base_mass,
        :history_mean_abs_delta_thickness,
        :history_mean_abs_delta_wet_mass,
        :history_mean_abs_delta_base_mass,
    ),
    monthly=(
        :monthly_mean_thickness,
        :monthly_mean_wet_mass,
        :monthly_mean_bulk_density,
        :monthly_mean_base_mass,
        :monthly_mean_ice_sheet_smb,
        :monthly_export_to_ice,
        :monthly_net_ice_sheet_forcing,
        :monthly_runoff,
    ),
    step=(:step_export_to_ice, :step_ice_sheet_smb),
)
const NETCDF_VARIABLES = unique(Symbol[var for group in values(OUTPUT_GROUPS) for var in group])
const FINAL_GRID_KEYS = OUTPUT_GROUPS.final
const LAYER_GRID_KEYS = OUTPUT_GROUPS.layers
const MONTHLY_GRID_KEYS = OUTPUT_GROUPS.monthly
@inline _grid_shape(layout) = size(layout.mask)

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
    delta_base = current_base .- previous.base_mass
    delta_smb = current_smb .- previous.smb_ice
    delta_runoff = current_runoff .- previous.runoff

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

function _scatter_final_grids(final_state, domain, deltas, layout)
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

function _finalize_monthly_grids(monthly_sums, monthly_count::Vector{Int32}, layout)
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

"""
Internal runtime types and helpers used by execute_run! and the NetCDF writer.
These are not part of the public API.
"""

# ---------------------------------------------------------------------------
# Internal structs
# ---------------------------------------------------------------------------

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

struct NetCDFWriter
    dataset::NCDataset
    vars::Dict{Symbol, Any}
    max_steps::Int
    cycles::Int
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

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

@inline completed_cycle_count(history::Vector{NamedTuple}, ::Symbol, cycles::Int) =
    isempty(history) ? 0 : min(history[end].cycle, cycles)

@inline cycle_metrics_schedule_label(stride::Int) =
    stride == 0 ? "final cycle only" : stride == 1 ? "every cycle" : "every $(stride) cycles + final"

@inline _looks_like_directory_path(path::AbstractString) =
    !isempty(path) && (endswith(path, '/') || endswith(path, '\\'))

function _slug(name::AbstractString)
    slug = strip(replace(lowercase(strip(String(name))), r"[^a-z0-9]+" => "_"), '_')
    return isempty(slug) ? "run" : slug
end

_default_output_dir(name::AbstractString) = joinpath(pwd(), "run_output", _slug(name))

function resolve_netcdf_path(options::RunOptions)
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
)

const HISTORY_OUTPUT_SPECS = (
    (output=:history_mean_thickness, record=:mean_thickness),
    (output=:history_mean_wet_mass, record=:mean_wet_mass),
    (output=:history_mean_bulk_density, record=:mean_bulk_density),
    (output=:history_mean_base_mass, record=:mean_base_mass),
    (output=:history_mean_abs_delta_thickness, record=:mean_abs_delta_thickness),
    (output=:history_mean_abs_delta_wet_mass, record=:mean_abs_delta_wet_mass),
    (output=:history_mean_abs_delta_base_mass, record=:mean_abs_delta_base_mass),
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
    options::RunOptions,
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
    out = Dict(spec.output => fill(NaN, cycles) for spec in HISTORY_OUTPUT_SPECS)
    for rec in history, spec in HISTORY_OUTPUT_SPECS
        out[spec.output][rec.cycle] = getfield(rec, spec.record)
    end
    return out
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
end

using Dates
using Printf
using Base.Threads: nthreads
using TimerOutputs: TimerOutputs
using Statistics: mean
using CSV

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

@inline function _finite_mean(data)
    finite = filter(isfinite, data)
    return isempty(finite) ? NaN : mean(finite)
end

function _delta_stats(data)
    finite = filter(isfinite, data)
    isempty(finite) && return (mean_signed=NaN, mean_abs=NaN, max_abs=NaN)
    abs_vals = abs.(finite)
    return (
        mean_signed=mean(finite),
        mean_abs=mean(abs_vals),
        max_abs=maximum(abs_vals),
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

function write_run_history_csv(out_path::AbstractString, history::Vector{NamedTuple})
    mkpath(dirname(out_path))
    # Build a NamedTuple with the published column labels from the history records.
    cols = NamedTuple{Tuple(Symbol(spec.label) for spec in HISTORY_CSV_SPECS)}(
        Tuple(getfield.(history, spec.key) for spec in HISTORY_CSV_SPECS)
    )
    CSV.write(out_path, cols)
end

function write_run_summary(
    out_path::AbstractString,
    options::RunOptions,
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

function print_run_report(
    io::IO,
    options::RunOptions,
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
    println(io, "Backend         : ", String(options.backend))
    println(io, "Cycles          : ", completed_cycle_count(history, status, options.cycles))
    println(io, "Status          : ", string(status))
    println(io, "Cycle metrics   : ", cycle_metrics_schedule_label(options.history_stride))
    println(io, @sprintf("Simulation wall : %.3f s", simulation_wall_sec))
    println(io, @sprintf("Run wall total  : %.3f s", run_wall_sec))
    haskey(timings.to.inner_timers, "model_step_wall") && println(io, @sprintf("Model step wall : %.3f s", TimerOutputs.time(timings.to.inner_timers["model_step_wall"]) * 1e-9))
    println(io, "Output NetCDF   : ", options.write_netcdf ? abspath(nc_path) : "skipped (--no-nc)")
    if options.write_outputs
        println(io, "History CSV     : $(abspath(history_csv_path))")
        println(io, "Summary         : $(abspath(summary_path))")
    else
        println(io, "File outputs    : skipped (--no-output)")
    end
    print_timing_summary(io, timings; total_wall_sec=run_wall_sec)
end
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

function execute_run!(
    domain::SnowpackDomain,
    forcing::SnowpackForcing;
    grid::Union{Nothing, SnowpackGrid}=nothing,
    options::RunOptions=RunOptions(),
    io::IO=stdout,
    timings::StepTimingStats=StepTimingStats(),
    run_wall_t0::Integer=time_ns(),
)
    ncol = column_count(domain)
    size(forcing.air_temperature, 1) == ncol || error("Forcing column count must match the domain column count.")
    spatial_grid = has_spatial_coords(grid)
    options.write_netcdf && !spatial_grid && error("NetCDF output requires a grid with spatial coordinates.")
    spatial_grid && length(grid.js) != ncol && error("Grid point count must match the domain column count.")

    selected = Set(options.netcdf_variables)
    is_gpu = options.backend == :gpu
    write_final_fields = options.write_netcdf && !isnothing(grid)
    need_step_outputs = options.write_netcdf && any(var -> var in selected, OUTPUT_GROUPS.step)
    need_monthly_outputs = options.write_netcdf && any(var -> var in selected, OUTPUT_GROUPS.monthly)
    need_layer_outputs = options.write_netcdf && any(var -> var in selected, OUTPUT_GROUPS.layers)
    need_last_cycle_smb_delta = options.write_netcdf && (:last_cycle_delta_ice_sheet_smb in selected)
    need_step_diagnostics = need_step_outputs || need_monthly_outputs

    initial_thickness_vec = if write_final_fields
        summary = allocate_cycle_summary_buffers(ncol)
        time_block!(timings, :prepare_initial_output_fields) do
            summarize_cycle_state!(
                summary.thickness,
                summary.wet_mass,
                summary.bulk_density,
                summary.base_mass,
                domain
            )
        end
        copy(summary.thickness)
    else
        Float64[]
    end

    domain, step_fields, workspace = _prepare_backend!(timings, domain, forcing; is_gpu=is_gpu)
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
                domain.Ntot,
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

    prev = allocate_cycle_summary_buffers(ncol)
    final = allocate_cycle_summary_buffers(ncol)
    backend_cycle_summary = allocate_cycle_summary_buffers(domain, ncol)
    step_summary = need_step_diagnostics ? allocate_summary_buffers(ncol) : nothing
    backend_step_summary = need_step_diagnostics ? allocate_summary_buffers(domain, ncol) : nothing
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

    time_block!(timings, :summarize_columns_initial) do
        summarize_cycle_state!(
            backend_cycle_summary.thickness,
            backend_cycle_summary.wet_mass,
            backend_cycle_summary.bulk_density,
            backend_cycle_summary.base_mass,
            domain
            )
        _copy_summary_fields!(prev, backend_cycle_summary, CYCLE_BUFFER_NAMES)
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
    progress = Progress(options.cycles; desc="Running cycles: ", output=io, showspeed=true)
    for cycle in 1:options.cycles
        for t in eachindex(forcing.time_values)
            month_idx = need_monthly_outputs ? (cycle - 1) * schedule.nmonth_per_cycle + schedule.step_month[t] : 0
            time_counted_block!(timings, :model_step_wall, ncol) do
                step!(domain, step_fields, t, workspace)
            end
            if need_step_diagnostics
                time_counted_block!(timings, :step_diagnostics, ncol) do
                    backend_step_summary === nothing && error("Missing backend summary buffers.")
                    summarize_domain_state!(
                        backend_step_summary.thickness,
                        backend_step_summary.wet_mass,
                        backend_step_summary.bulk_density,
                        backend_step_summary.base_mass,
                        backend_step_summary.smb_ice,
                        backend_step_summary.liquid_water,
                        backend_step_summary.runoff,
                        domain
            )
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

        time_block!(timings, :summarize_columns_cycle) do
            summarize_cycle_state!(
                backend_cycle_summary.thickness,
                backend_cycle_summary.wet_mass,
                backend_cycle_summary.bulk_density,
                backend_cycle_summary.base_mass,
                domain
            )
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
                _update_cycle_smb_delta!(deltas.ice_sheet_smb, previous_cycle_smb_ice, domain)
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
            _scatter_final_grids(prev, domain, deltas, grid)
        end
        layer_grids = if need_layer_outputs
            final_domain = is_gpu ? time_block!(timings, :gpu_transfer) do
                cpu_domain(domain)
            end : domain
            time_block!(timings, :collect_final_layer_grids) do
                collect_final_layer_grids(final_domain, grid.js, grid.is, _grid_shape(grid), final_domain.Ntot)
            end
        else
            _empty_layer_grids()
        end
        monthly_grids = need_monthly_outputs ? time_block!(timings, :aggregate_monthly_outputs) do
            _finalize_monthly_grids(monthly_sums, monthly_count, grid)
        end : empty_monthly_grids()
        time_block!(timings, :write_netcdf) do
            finalize_netcdf!(writer, final_grids, layer_grids, history, monthly_grids, :cycles, length(history), steps_written)
        end
    end

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
        domain=domain,
    )
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
    sim::Simulation{<:BESSIModel};
    options::SimulationOptions=sim.options,
    output::OutputOptions=sim.output,
    io::IO=stdout,
)
    options = RunOptions(
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

    result = execute_run!(sim.model.domain, sim.forcing; grid=sim.model.grid, options=options, io=io)
    if options.backend == :gpu
        _copy_domain_state!(sim.model.domain, cpu_domain(result.domain))
    end
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

function run!(::Simulation{<:PDDModel}; io::IO=stdout)
    error("PDDModel is not yet implemented. Physics coming soon.")
end

function run!(::Simulation{<:ITMModel}; io::IO=stdout)
    error("ITMModel is not yet implemented. Physics coming soon.")
end
