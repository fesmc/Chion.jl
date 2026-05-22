# Output grids, NetCDF schema, and NetCDF writer helpers for Simulation runs.

const OUTPUT_VARIABLE_SPECS = (
    (key=:final_thickness, group=:final, source=:final_state, field=:thickness, name="final_thickness", dims=("x", "y"), long_name="Final snow thickness", units="m", integer=false),
    (key=:final_wet_mass, group=:final, source=:final_state, field=:wet_mass, name="final_wet_mass", dims=("x", "y"), long_name="Final snow wet mass", units="mmWE", integer=false),
    (key=:final_bulk_density, group=:final, source=:final_state, field=:bulk_density, name="final_bulk_density", dims=("x", "y"), long_name="Final bulk snow density", units="kg m-3", integer=false),
    (key=:final_base_mass, group=:final, source=:final_state, field=:base_mass, name="final_base_mass", dims=("x", "y"), long_name="Cumulative firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:final_ice_sheet_smb, group=:final, source=:domain, field=:smb_ice, name="final_ice_sheet_smb", dims=("x", "y"), long_name="Cumulative net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:final_runoff, group=:final, source=:domain, field=:runoff, name="final_runoff", dims=("x", "y"), long_name="Final cumulative runoff", units="mmWE", integer=false),
    (key=:last_year_delta_thickness, group=:final, source=:deltas, field=:thickness, name="last_year_delta_thickness", dims=("x", "y"), long_name="Last year snow-thickness change", units="m", integer=false),
    (key=:last_year_delta_wet_mass, group=:final, source=:deltas, field=:wet_mass, name="last_year_delta_wet_mass", dims=("x", "y"), long_name="Last year wet-mass change", units="mmWE", integer=false),
    (key=:last_year_delta_base_mass, group=:final, source=:deltas, field=:base_mass, name="last_year_delta_base_mass", dims=("x", "y"), long_name="Last year firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:last_year_delta_ice_sheet_smb, group=:final, source=:deltas, field=:ice_sheet_smb, name="last_year_delta_ice_sheet_smb", dims=("x", "y"), long_name="Last year net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:n_active, group=:layers, name="n_active", dims=("x", "y"), long_name="Number of active Chion layers", units="", integer=true),
    (key=:layer_density, group=:layers, name="layer_density", dims=("layer", "x", "y"), long_name="Final Chion layer density", units="kg m-3", integer=false),
    (key=:layer_thickness, group=:layers, name="layer_thickness", dims=("layer", "x", "y"), long_name="Final Chion layer thickness", units="m", integer=false),
    (key=:layer_snow_mass, group=:layers, name="layer_snow_mass", dims=("layer", "x", "y"), long_name="Final Chion layer snow mass", units="kg m-2", integer=false),
    (key=:layer_liquid_mass, group=:layers, name="layer_liquid_mass", dims=("layer", "x", "y"), long_name="Final Chion layer liquid-water mass", units="kg m-2", integer=false),
    (key=:layer_temperature_c, group=:layers, name="layer_temperature_c", dims=("layer", "x", "y"), long_name="Final Chion layer temperature", units="C", integer=false),
    (key=:history_mean_thickness, group=:history, record=:mean_thickness, name="history_mean_thickness", dims=("year",), long_name="Year-mean snow thickness", units="m", integer=false),
    (key=:history_mean_wet_mass, group=:history, record=:mean_wet_mass, name="history_mean_wet_mass", dims=("year",), long_name="Year-mean snow wet mass", units="mmWE", integer=false),
    (key=:history_mean_bulk_density, group=:history, record=:mean_bulk_density, name="history_mean_bulk_density", dims=("year",), long_name="Year-mean bulk snow density", units="kg m-3", integer=false),
    (key=:history_mean_base_mass, group=:history, record=:mean_base_mass, name="history_mean_base_mass", dims=("year",), long_name="Year-mean firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:history_mean_abs_delta_thickness, group=:history, record=:mean_abs_delta_thickness, name="history_mean_abs_delta_thickness", dims=("year",), long_name="Year mean absolute snow-thickness change", units="m", integer=false),
    (key=:history_mean_abs_delta_wet_mass, group=:history, record=:mean_abs_delta_wet_mass, name="history_mean_abs_delta_wet_mass", dims=("year",), long_name="Year mean absolute wet-mass change", units="mmWE", integer=false),
    (key=:history_mean_abs_delta_base_mass, group=:history, record=:mean_abs_delta_base_mass, name="history_mean_abs_delta_base_mass", dims=("year",), long_name="Year mean absolute firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:monthly_mean_thickness, group=:monthly, source=:summary, field=:thickness, aggregate=:mean, name="monthly_mean_thickness", dims=("month", "x", "y"), long_name="Monthly mean snow thickness", units="m", integer=false),
    (key=:monthly_mean_wet_mass, group=:monthly, source=:summary, field=:wet_mass, aggregate=:mean, name="monthly_mean_wet_mass", dims=("month", "x", "y"), long_name="Monthly mean snow wet mass", units="mmWE", integer=false),
    (key=:monthly_mean_bulk_density, group=:monthly, source=:summary, field=:bulk_density, aggregate=:mean, name="monthly_mean_bulk_density", dims=("month", "x", "y"), long_name="Monthly mean bulk snow density", units="kg m-3", integer=false),
    (key=:monthly_mean_base_mass, group=:monthly, source=:summary, field=:base_mass, aggregate=:sum, name="monthly_mean_base_mass", dims=("month", "x", "y"), long_name="Monthly mean cumulative firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:monthly_smb, group=:monthly, source=:computed, field=:climatic_smb, aggregate=:sum, name="monthly_smb", dims=("month", "x", "y"), long_name="Monthly climatic surface mass balance", units="mmWE", integer=false),
    (key=:monthly_export_to_ice, group=:monthly, source=:delta, field=:base_mass, aggregate=:sum, name="monthly_export_to_ice", dims=("month", "x", "y"), long_name="Monthly firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:monthly_net_ice_sheet_forcing, group=:monthly, source=:delta, field=:smb_ice, aggregate=:sum, name="monthly_net_ice_sheet_forcing", dims=("month", "x", "y"), long_name="Monthly net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:monthly_runoff, group=:monthly, source=:delta, field=:runoff, aggregate=:sum, name="monthly_runoff", dims=("month", "x", "y"), long_name="Monthly runoff production", units="mmWE", integer=false),
    (key=:monthly_melt, group=:monthly, source=:delta, field=:melt, aggregate=:sum, name="monthly_melt", dims=("month", "x", "y"), long_name="Monthly melt production", units="mmWE", integer=false),
    (key=:monthly_refreezing, group=:monthly, source=:delta, field=:refreezing, aggregate=:sum, name="monthly_refreezing", dims=("month", "x", "y"), long_name="Monthly refreezing", units="mmWE", integer=false),
    (key=:monthly_sublimation, group=:monthly, source=:delta, field=:sublimation, aggregate=:sum, name="monthly_sublimation", dims=("month", "x", "y"), long_name="Monthly sublimation mass loss from turbulent latent heat flux", units="mmWE", integer=false),
    (key=:monthly_mean_latent_heat_flux, group=:monthly, source=:delta, field=:latent_heat_flux_sum, aggregate=:mean, name="monthly_mean_latent_heat_flux", dims=("month", "x", "y"), long_name="Monthly mean turbulent latent heat flux", units="W m-2", integer=false),
    (key=:monthly_mean_albedo, group=:monthly, source=:summary, field=:albedo, aggregate=:mean, name="monthly_mean_albedo", dims=("month", "x", "y"), long_name="Monthly mean surface albedo", units="1", integer=false),
    (key=:step_export_to_ice, group=:step, source=:delta, field=:base_mass, name="step_export_to_ice", dims=("step", "x", "y"), long_name="Annual firn mass exported to the ice model for each written output interval", units="mmWE", integer=false),
    (key=:step_ice_sheet_smb, group=:step, source=:delta, field=:smb_ice, name="step_ice_sheet_smb", dims=("step", "x", "y"), long_name="Annual net mass forcing to the ice sheet for each written output interval", units="mmWE", integer=false),
    (key=:step_pdd, group=:step, source=:delta, field=:pdd, name="step_pdd", dims=("step", "x", "y"), long_name="Annual positive degree days for each written output interval", units="degC d", integer=false),
    (key=:daily_latent_heat_flux, group=:daily, source=:delta, field=:latent_heat_flux_sum, name="daily_latent_heat_flux", dims=("day", "x", "y"), long_name="Daily turbulent latent heat flux", units="W m-2", integer=false),
)

_specs_for(group::Symbol) = Tuple(spec for spec in OUTPUT_VARIABLE_SPECS if spec.group === group)
_keys_for(specs) = Tuple(spec.key for spec in specs)

const FINAL_OUTPUT_SPECS = _specs_for(:final)
const LAYER_OUTPUT_SPECS = _specs_for(:layers)
const HISTORY_OUTPUT_SPECS = _specs_for(:history)
const MONTHLY_OUTPUT_SPECS = _specs_for(:monthly)
const STEP_OUTPUT_SPECS = _specs_for(:step)
const DAILY_OUTPUT_SPECS = _specs_for(:daily)

const OUTPUT_GROUPS = (
    final=_keys_for(FINAL_OUTPUT_SPECS),
    layers=_keys_for(LAYER_OUTPUT_SPECS),
    history=_keys_for(HISTORY_OUTPUT_SPECS),
    monthly=_keys_for(MONTHLY_OUTPUT_SPECS),
    step=_keys_for(STEP_OUTPUT_SPECS),
    daily=_keys_for(DAILY_OUTPUT_SPECS),
)
const NETCDF_VARIABLES = unique(Symbol[var for group in values(OUTPUT_GROUPS) for var in group])
const FINAL_GRID_KEYS = OUTPUT_GROUPS.final
const LAYER_GRID_KEYS = OUTPUT_GROUPS.layers
const MONTHLY_GRID_KEYS = OUTPUT_GROUPS.monthly
const DAILY_GRID_KEYS = OUTPUT_GROUPS.daily
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
empty_daily_grids() = NamedTuple{DAILY_GRID_KEYS}(ntuple(_ -> Array{Float64}(undef, 0, 0, 0), length(DAILY_GRID_KEYS)))

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
_allocate_daily_vectors(active::Bool, ncol::Int) = NamedTuple{DAILY_GRID_KEYS}(ntuple(_ -> active ? zeros(Float64, ncol) : Float64[], length(DAILY_GRID_KEYS)))
_allocate_monthly_sums(active::Bool, nmonth_total::Int, ncol::Int) = NamedTuple{MONTHLY_GRID_KEYS}(ntuple(_ -> active ? zeros(Float64, nmonth_total, ncol) : Matrix{Float64}(undef, 0, 0), length(MONTHLY_GRID_KEYS)))

function _step_output_grids(step_vectors, layout)
    return NamedTuple{OUTPUT_GROUPS.step}(ntuple(i -> scatter_to_grid(getfield(step_vectors, OUTPUT_GROUPS.step[i]), layout.js, layout.is, _grid_shape(layout)), length(OUTPUT_GROUPS.step)))
end

function _daily_output_grids(daily_vectors, layout)
    return NamedTuple{DAILY_GRID_KEYS}(ntuple(i -> scatter_to_grid(getfield(daily_vectors, DAILY_GRID_KEYS[i]), layout.js, layout.is, _grid_shape(layout)), length(DAILY_GRID_KEYS)))
end

function _reset_step_vectors!(step_vectors)
    for key in OUTPUT_GROUPS.step
        isempty(getfield(step_vectors, key)) || fill!(getfield(step_vectors, key), 0.0)
    end
    return
end

function _summary_deltas(summary, previous)
    return (
        base_mass=summary.base_mass .- previous.base_mass,
        wet_mass=summary.wet_mass .- previous.wet_mass,
        smb_ice=summary.smb_ice .- previous.smb_ice,
        runoff=summary.runoff .- previous.runoff,
        pdd=summary.pdd .- previous.pdd,
        melt=summary.melt .- previous.melt,
        refreezing=summary.refreezing .- previous.refreezing,
        sublimation=summary.sublimation .- previous.sublimation,
        latent_heat_flux_sum=summary.latent_heat_flux_sum .- previous.latent_heat_flux_sum,
    )
end

function _diagnostic_vector(spec, summary, delta, computed)
    spec.source === :summary && return getfield(summary, spec.field)
    spec.source === :delta && return getfield(delta, spec.field)
    spec.source === :computed && return getfield(computed, spec.field)
    error("Unsupported diagnostic source `$(spec.source)` for `$(spec.key)`.")
end

function _copy_previous_summary!(previous, summary)
    for key in propertynames(previous)
        getfield(previous, key) .= getfield(summary, key)
    end
    return nothing
end

function _accumulate_step_diagnostics!(
    summary,
    previous,
    monthly_sums,
    step_vectors,
    daily_vectors,
    month_idx::Int,
    need_monthly_outputs::Bool,
    need_step_outputs::Bool,
    need_daily_outputs::Bool,
)
    delta = _summary_deltas(summary, previous)
    computed = (
        climatic_smb=delta.wet_mass .+ delta.smb_ice,
    )

    if need_monthly_outputs
        for spec in MONTHLY_OUTPUT_SPECS
            getfield(monthly_sums, spec.key)[month_idx, :] .+= _diagnostic_vector(spec, summary, delta, computed)
        end
    end
    if need_step_outputs
        for spec in STEP_OUTPUT_SPECS
            getfield(step_vectors, spec.key) .+= _diagnostic_vector(spec, summary, delta, computed)
        end
    end
    if need_daily_outputs
        for spec in DAILY_OUTPUT_SPECS
            getfield(daily_vectors, spec.key) .= _diagnostic_vector(spec, summary, delta, computed)
        end
    end
    _copy_previous_summary!(previous, summary)
    return
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
    max_days::Int
    years::Int
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

const NC_SPECS = OUTPUT_VARIABLE_SPECS

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
    month_of_year::Vector{Int32},
    source_month_code::Vector{Int32},
    annual_output_source_indices::Vector{Int32},
    annual_output_source_codes::Vector{Int32},
)
    mkpath(dirname(netcdf_path))
    isdir(netcdf_path) && error("NetCDF output path '$(abspath(netcdf_path))' is a directory; pass a file path ending in `.nc`.")
    ny, nx = _grid_shape(layout)
    max_steps = options.years * length(annual_output_source_indices)
    max_days = options.years * length(time_values)
    selected = Set(options.netcdf_variables)
    step_year = Int32[yr for yr in 1:options.years for _ in annual_output_source_indices]
    step_source_index = Int32[idx for _ in 1:options.years for idx in annual_output_source_indices]
    step_source_code = Int32[code for _ in 1:options.years for code in annual_output_source_codes]

    ds = NCDataset(netcdf_path, "c")
    for (name, len) in (("x", nx), ("y", ny), ("layer", max(nlayer, 1)), ("year", options.years), ("month", length(month_of_year)), ("day", max_days), ("point", length(layout.js)), ("step", max_steps))
        defDim(ds, name, len)
    end

    for spec in (
        (name="x", dims=("x",), meta=(key=:x, long_name="X coordinate", units="km", integer=false), value=layout.x),
        (name="y", dims=("y",), meta=(key=:y, long_name="Y coordinate", units="km", integer=false), value=layout.y),
        (name="layer", dims=("layer",), meta=(key=:layer, long_name="Chion internal layer index from surface downward", units="", integer=true), value=Int32.(collect(1:max(nlayer, 1)))),
        (name="year", dims=("year",), meta=(key=:year, long_name="Repeated forcing year index", units="", integer=true), value=Int32.(collect(1:options.years))),
        (name="month", dims=("month",), meta=(key=:month, long_name="Sequential monthly output index", units="", integer=true), value=Int32.(collect(1:length(month_of_year)))),
        (name="day", dims=("day",), meta=(key=:day, long_name="Sequential daily output index", units="", integer=true), value=Int32.(collect(1:max_days))),
        (name="month_of_year", dims=("month",), meta=(key=:month_of_year, long_name="Calendar month of the repeated forcing", units="", integer=true), value=month_of_year),
        (name="source_month_code", dims=("month",), meta=(key=:source_month_code, long_name="Source forcing month code YYYYMM", units="", integer=true), value=source_month_code),
        (name="point", dims=("point",), meta=(key=:point, long_name="Compact valid cell index", units="", integer=true), value=Int32.(collect(1:length(layout.js)))),
        (name="point_j", dims=("point",), meta=(key=:point_j, long_name="1-based y-index for each compact valid cell", units="", integer=true), value=Int32.(layout.js)),
        (name="point_i", dims=("point",), meta=(key=:point_i, long_name="1-based x-index for each compact valid cell", units="", integer=true), value=Int32.(layout.is)),
        (name="point_y_km", dims=("point",), meta=(key=:point_y_km, long_name="Y coordinate for each compact valid cell", units="km", integer=false), value=layout.y[layout.js]),
        (name="point_x_km", dims=("point",), meta=(key=:point_x_km, long_name="X coordinate for each compact valid cell", units="km", integer=false), value=layout.x[layout.is]),
        (name="step", dims=("step",), meta=(key=:step, long_name="Sequential yearly output index across repeated years", units="", integer=true), value=Int32.(collect(1:max_steps))),
        (name="step_year", dims=("step",), meta=(key=:step_year, long_name="Repeated forcing year index for each yearly output", units="", integer=true), value=step_year),
        (name="step_source_index", dims=("step",), meta=(key=:step_source_index, long_name="1-based index of the last forcing step included in each yearly output", units="", integer=true), value=step_source_index),
        (name="step_source_code", dims=("step",), meta=(key=:step_source_code, long_name="Source forcing timestamp code YYYYMMDDHH for the final step included in each yearly output", units="", integer=true), value=step_source_code),
        (name="domain_mask", dims=("x", "y"), meta=(key=:domain_mask, long_name="Domain mask", units="1", integer=false), value=Float32.(permutedims(layout.mask, (2, 1)))),
        (name="initial_thickness", dims=("x", "y"), meta=(key=:initial_thickness, long_name="Initial snow thickness", units="m", integer=false), value=Float32.(permutedims(initial_thickness, (2, 1)))),
    )
        _write_nc_var!(ds, spec.name, spec.dims, spec.meta, spec.value)
    end

    vars = _define_selected_nc_variables!(ds, selected)
    any(key -> key in selected, OUTPUT_GROUPS.step) && (vars[:step_valid] = _def_nc_var(ds, "step_valid", ("step",), (key=:step_valid, long_name="1 where a yearly output record was completed and written, 0 for unused trailing slots", units="", integer=true)))
    any(key -> key in selected, OUTPUT_GROUPS.daily) && (vars[:day_valid] = _def_nc_var(ds, "day_valid", ("day",), (key=:day_valid, long_name="1 where a daily output record was completed and written, 0 for unused trailing slots", units="", integer=true)))

    ds.attrib["title"] = options.name
    ds.attrib["source_model"] = "Chion"
    ds.attrib["input_label"] = isempty(options.input_label) ? "not provided" : options.input_label
    ds.attrib["forcing_start"] = string(first(time_values))
    ds.attrib["forcing_end"] = string(last(time_values))
    ds.attrib["years_completed"] = "pending"
    ds.attrib["status"] = "pending"
    ds.attrib["created"] = string(now())
    return NetCDFWriter(ds, vars, max_steps, max_days, options.years)
end

@inline _write_dataset_var!(var, data::AbstractVector) = (var[:] = eltype(var) <: Integer ? data : Float32.(data))
@inline _write_dataset_var!(var, data::AbstractMatrix) = (var[:, :] = eltype(var) <: Integer ? permutedims(data, (2, 1)) : Float32.(permutedims(data, (2, 1))))
@inline _write_dataset_var!(var, data::AbstractArray{<:Real, 3}) = (var[:, :, :] = Float32.(permutedims(data, (1, 3, 2))))

maybe_write_step_output!(writer::NetCDFWriter, step_index::Int, key::Symbol, data::AbstractMatrix{<:Real}) =
    (haskey(writer.vars, key) && (writer.vars[key][step_index, :, :] = Float32.(permutedims(data, (2, 1))); nothing))

maybe_write_daily_output!(writer::NetCDFWriter, day_index::Int, key::Symbol, data::AbstractMatrix{<:Real}) =
    (haskey(writer.vars, key) && (writer.vars[key][day_index, :, :] = Float32.(permutedims(data, (2, 1))); nothing))

maybe_write_output!(writer::NetCDFWriter, key::Symbol, data) = (haskey(writer.vars, key) && (_write_dataset_var!(writer.vars[key], data); nothing))

function _history_vectors(history::Vector{NamedTuple}, years::Int)
    out = Dict(spec.key => fill(NaN, years) for spec in HISTORY_OUTPUT_SPECS)
    for rec in history, spec in HISTORY_OUTPUT_SPECS
        out[spec.key][rec.year] = getfield(rec, spec.record)
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
    years_completed::Int,
    days_written::Int,
    steps_written::Int,
)
    _write_output_group!(writer, FINAL_GRID_KEYS, final_grids)
    _write_output_group!(writer, LAYER_GRID_KEYS, layer_grids)
    for (key, values) in _history_vectors(history, writer.years)
        maybe_write_output!(writer, key, values)
    end
    _write_output_group!(writer, MONTHLY_GRID_KEYS, monthly_grids)
    if haskey(writer.vars, :step_valid)
        step_valid = zeros(Int32, writer.max_steps)
        step_valid[1:steps_written] .= 1
        writer.vars[:step_valid][:] = step_valid
    end
    if haskey(writer.vars, :day_valid)
        day_valid = zeros(Int32, writer.max_days)
        day_valid[1:days_written] .= 1
        writer.vars[:day_valid][:] = day_valid
    end
    writer.dataset.attrib["years_completed"] = string(years_completed)
    writer.dataset.attrib["status"] = string(status)
    writer.dataset.attrib["days_written"] = string(days_written)
    writer.dataset.attrib["steps_written"] = string(steps_written)
    close(writer.dataset)
end

function _prepare_output_schedule(time_values::Vector{DateTime}, years::Int)
    write_output = falses(length(time_values))
    output_slot = zeros(Int, length(time_values))
    source_indices = Int32[]
    source_codes = Int32[]
    source_years = unique(year.(time_values))
    for (slot, yr) in enumerate(source_years)
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
        annual_output=(write_output=write_output, output_slot=output_slot, source_indices=source_indices, source_codes=source_codes, source_years=source_years),
        nstep_per_year=length(time_values),
        step_month=[month_lookup[(year(t), month(t))] for t in time_values],
        nmonth_per_year=length(month_keys),
        nmonth_total=years * length(month_keys),
        month_of_year=Int32[key[2] for _ in 1:years for key in month_keys],
        source_month_code=Int32[key[1] * 100 + key[2] for _ in 1:years for key in month_keys],
    )
end
