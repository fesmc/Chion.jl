# Output grids, NetCDF schema, and output helpers for Simulation runs.

const io_dict = Dict{Symbol, NamedTuple}(
    :thickness => (name="thickness", long_name="Snow thickness", units="m"),
    :wet_mass => (name="wet_mass", long_name="Snow wet mass", units="mmWE"),
    :bulk_density => (name="bulk_density", long_name="Bulk snow density", units="kg m-3"),
    :liquid_water => (name="liquid_water", long_name="Liquid water mass", units="kg m-2"),
    :mass_base => (name="mass_base", long_name="Firn mass exported to the ice model", units="mmWE"),
    :smb_ice => (name="smb_ice", long_name="Net mass forcing to the ice sheet", units="mmWE"),
    :runoff => (name="runoff", long_name="Cumulative runoff", units="mmWE"),
    :melt => (name="melt", long_name="Cumulative melt", units="mmWE"),
    :refreezing => (name="refreezing", long_name="Cumulative refreezing", units="mmWE"),
    :sublimation => (name="sublimation", long_name="Cumulative sublimation", units="mmWE"),
    :latent_heat_flux => (name="latent_heat_flux", long_name="Monthly mean turbulent latent heat flux", units="W m-2"),
    :latent_heat_flux_sum => (name="latent_heat_flux_sum", long_name="Integrated turbulent latent heat flux", units="W m-2"),
    :Tsrf => (name="Tsrf", long_name="Surface temperature", units="K"),
    :albedo => (name="albedo", long_name="Surface albedo", units="1"),
    :N => (name="N", long_name="Number of active snow layers", units="1"),
    :mass => (name="mass", long_name="Layer snow mass", units="kg m-2"),
    :mass_w => (name="mass_w", long_name="Layer liquid-water mass", units="kg m-2"),
    :density => (name="density", long_name="Layer density", units="kg m-3"),
    :temperature => (name="temperature", long_name="Layer temperature", units="K"),
)

const DEFAULT_STATE_OUTPUT_VARS = [:thickness, :wet_mass, :bulk_density, :mass_base, :smb_ice, :runoff, :melt, :refreezing, :sublimation, :albedo]
const STATE_FIELD_OUTPUT_VARS = [:mass, :mass_w, :density, :temperature, :N, :liquid_water, :latent_heat_flux_sum, :Tsrf]
const MONTHLY_OUTPUT_VARS = [:smb_ice, :runoff, :melt, :refreezing, :sublimation, :latent_heat_flux, :albedo]
const NETCDF_VARIABLES = unique(vcat(DEFAULT_STATE_OUTPUT_VARS, STATE_FIELD_OUTPUT_VARS))
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

function collect_final_layer_grids(
    domain::AbstractSnowpackDomain,
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
struct NetcdfOutput
    dataset::NCDataset
    vars::Dict{Symbol, Any}
    buffer_point_float::Matrix{Float32}
    buffer_point_int::Matrix{Int32}
    buffer_layer_float::Array{Float32, 3}
    buffer_time_point_float::Array{Float32, 3}
    max_steps::Int
    max_days::Int
    years::Int
end

_meta_value(meta, key::Symbol, default) = hasproperty(meta, key) ? getproperty(meta, key) : default

function state_output_vars(selected)
    isempty(selected) && return Symbol[]
    vars = Symbol[]
    for key in selected
        if key == :monthly
            append!(vars, MONTHLY_OUTPUT_VARS)
        elseif key == :all
            append!(vars, DEFAULT_STATE_OUTPUT_VARS)
        elseif key in (:mass, :mass_w, :density, :temperature, :N, :thickness, :wet_mass, :bulk_density, :liquid_water, :mass_base, :smb_ice, :runoff, :melt, :refreezing, :sublimation, :latent_heat_flux_sum, :Tsrf, :albedo)
            push!(vars, key)
        end
    end
    return unique(vars)
end

function init_state_netcdf(path::AbstractString, options, time_values::Vector{DateTime}, domain, state, vars::Vector{Symbol}; ntime::Integer=options.years * length(time_values), nlayer::Integer=getproperty(state, :Ntot))
    mkpath(dirname(path))
    ds = NCDataset(path, "c")
    ny, nx = _grid_shape(domain)
    defDim(ds, "t", Int(ntime))
    defDim(ds, "x", nx)
    defDim(ds, "y", ny)
    defDim(ds, "layer", Int(nlayer))
    _write_nc_var!(ds, "x", ("x",), (key=:x, long_name="X coordinate", units="km", integer=false), domain.x)
    _write_nc_var!(ds, "y", ("y",), (key=:y, long_name="Y coordinate", units="km", integer=false), domain.y)
    _write_nc_var!(ds, "domain_mask", ("x", "y"), (key=:domain_mask, long_name="Domain mask", units="1", integer=false), Float32.(permutedims(domain.mask, (2, 1))))
    handles = Dict{Symbol, Any}()
    for key in vars
        values = key == :latent_heat_flux && hasfield(typeof(state), :latent_heat_flux_sum) ?
            getfield(state, :latent_heat_flux_sum) :
            getfield(state, key)
        meta = get(io_dict, key, (name=String(key), long_name=String(key), units="", integer=eltype(values) <: Integer))
        if values isa AbstractVector
            handles[key] = _def_nc_var(ds, String(_meta_value(meta, :name, String(key))), ("t", "x", "y"), (key=key, long_name=_meta_value(meta, :long_name, String(key)), units=_meta_value(meta, :units, ""), integer=eltype(values) <: Integer))
        elseif values isa AbstractMatrix
            handles[key] = _def_nc_var(ds, String(_meta_value(meta, :name, String(key))), ("t", "layer", "x", "y"), (key=key, long_name=_meta_value(meta, :long_name, String(key)), units=_meta_value(meta, :units, ""), integer=eltype(values) <: Integer))
        end
    end
    ds.attrib["title"] = options.name
    ds.attrib["source_model"] = "Chion"
    ds.attrib["created"] = string(now())
    return NetcdfOutput(
        ds,
        handles,
        fill(NaN32, nx, ny),
        zeros(Int32, nx, ny),
        fill(NaN32, Int(nlayer), nx, ny),
        fill(NaN32, min(Int(ntime), 240), nx, ny),
        ntime,
        ntime,
        options.years,
    )
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
    for token in tokens
        stripped = strip(token)
        isempty(stripped) && continue
        key = Symbol(lowercase(stripped))
        key == :n && (key = :N)
        key == :tsrf && (key = :Tsrf)
        if key == :all
            append!(selected, NETCDF_VARIABLES)
        elseif key == :none
            continue
        elseif key == :monthly
            push!(selected, :monthly)
        elseif key in NETCDF_VARIABLES || key in MONTHLY_OUTPUT_VARS
            push!(selected, key)
        else
            error("Unsupported NetCDF variable selector '$token'. Use `all`, `none`, or a CurrentState field name.")
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

function _scatter_vector_to_grid!(dest::Matrix{Float32}, values, layout)
    fill!(dest, NaN32)
    host = _host_vector(values; copy_array=false)
    @inbounds for col in eachindex(layout.js)
        dest[layout.is[col], layout.js[col]] = Float32(host[col])
    end
    return dest
end

function _scatter_vector_to_grid!(dest::Matrix{Int32}, values, layout)
    fill!(dest, Int32(0))
    host = Array(values)
    @inbounds for col in eachindex(layout.js)
        dest[layout.is[col], layout.js[col]] = Int32(host[col])
    end
    return dest
end

function _scatter_matrix_to_grid!(dest::Array{Float32, 3}, values, layout)
    fill!(dest, NaN32)
    host = Array(values)
    @inbounds for col in eachindex(layout.js), layer in axes(host, 1)
        dest[layer, layout.is[col], layout.js[col]] = Float32(host[layer, col])
    end
    return dest
end

function _scatter_monthly_matrix_to_grid!(dest::Array{Float32, 3}, values::AbstractMatrix, first_row::Int, nrecord::Int, layout)
    fill!(dest, NaN32)
    host = Array(values)
    @inbounds for row in 1:nrecord, col in eachindex(layout.js)
        dest[row, layout.is[col], layout.js[col]] = Float32(host[first_row + row - 1, col])
    end
    return view(dest, 1:nrecord, :, :)
end

function _state_output_grid(state, key::Symbol, layout)
    values = getfield(state, key)
    if values isa AbstractVector
        return scatter_to_grid(_host_vector(values; copy_array=true), layout.js, layout.is, _grid_shape(layout))
    elseif values isa AbstractMatrix
        host = Array(values)
        grids = fill(NaN, size(host, 1), _grid_shape(layout)...)
        @inbounds for col in eachindex(layout.js), layer in axes(host, 1)
            grids[layer, layout.js[col], layout.is[col]] = host[layer, col]
        end
        return grids
    end
    error("Cannot write state field `$key` with type $(typeof(values)).")
end

function write_monthly_year_nc!(out::NetcdfOutput, state, first_record_index::Int, layout)
    nrecord = state.count
    nrecord == 0 && return out
    chunk_len = size(out.buffer_time_point_float, 1)
    for key in keys(out.vars)
        hasfield(typeof(state), key) || continue
        values = getfield(state, key)
        var = out.vars[key]
        offset = 0
        while offset < nrecord
            n = min(chunk_len, nrecord - offset)
            record_start = first_record_index + offset
            record_stop = record_start + n - 1
            var[record_start:record_stop, :, :] =
                _scatter_monthly_matrix_to_grid!(out.buffer_time_point_float, values, offset + 1, n, layout)
            offset += n
        end
    end
    return out
end

function write_nc!(out::NetcdfOutput, state, record_index::Int, layout)
    for key in keys(out.vars)
        hasfield(typeof(state), key) || continue
        values = getfield(state, key)
        var = out.vars[key]
        if values isa AbstractVector
            if eltype(values) <: Integer
                var[record_index, :, :] = _scatter_vector_to_grid!(out.buffer_point_int, values, layout)
            else
                var[record_index, :, :] = _scatter_vector_to_grid!(out.buffer_point_float, values, layout)
            end
        elseif values isa AbstractMatrix
            var[record_index, :, :, :] = _scatter_matrix_to_grid!(out.buffer_layer_float, values, layout)
        end
    end
    return out
end

function close_output!(out::NetcdfOutput, status::Symbol, records_written::Int)
    out.dataset.attrib["status"] = string(status)
    out.dataset.attrib["records_written"] = string(records_written)
    close(out.dataset)
    return nothing
end
