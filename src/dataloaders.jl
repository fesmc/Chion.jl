"""
Generic NetCDF forcing loaders for Chion's simulation-first API.
"""

function infer_dt_days(time_values::Vector{DateTime})
    ntime = length(time_values)
    ntime == 0 && error("forcing file contains no timesteps.")
    ntime == 1 && return [1.0]
    dt_days = Vector{Float64}(undef, ntime)
    for idx in 1:ntime-1
        dt_days[idx] = Dates.value(time_values[idx + 1] - time_values[idx]) / 86_400_000
    end
    dt_days[end] = dt_days[end - 1]
    return dt_days
end

function _read_time_values(ds::NCDataset, time_name::AbstractString, ntime::Int)
    if haskey(ds, time_name)
        decoded = ds[time_name][:]
        if eltype(decoded) <: DateTime
            return DateTime.(decoded)
        end
    end
    if haskey(ds, "TIME")
        decoded = ds["TIME"][:]
        if eltype(decoded) <: DateTime
            return DateTime.(decoded)
        end
    end
    if all(name -> haskey(ds, name), ("YYYY", "MM", "DD", "HH"))
        yyyy = round.(Int, vec(ds["YYYY"].var[:]))
        mm = round.(Int, vec(ds["MM"].var[:]))
        dd = round.(Int, vec(ds["DD"].var[:]))
        hh = round.(Int, vec(ds["HH"].var[:]))
        length(yyyy) == ntime || error("Time metadata length does not match forcing fields.")
        return [DateTime(yyyy[i], mm[i], dd[i], hh[i]) for i in 1:ntime]
    end
    base = DateTime(2001, 1, 1, 12)
    return [base + Dates.Day(idx - 1) for idx in 1:ntime]
end

function _infer_time_length(ds::NCDataset, time_name::AbstractString, data, ny::Int, nx::Int, name::AbstractString)
    haskey(ds, time_name) && return length(ds[time_name].var[:])
    haskey(ds, "TIME") && return length(ds["TIME"].var[:])
    haskey(ds, "YYYY") && return length(ds["YYYY"].var[:])

    A = _drop_singleton_layer(data, name)
    dims = size(A)
    ndims(A) == 3 || error("Cannot infer time length for `$name` from shape $(dims).")
    candidates = Int[d for d in dims if d != nx && d != ny]
    length(candidates) == 1 && return only(candidates)
    error("Cannot infer time dimension for `$name`; pass `time_name` or provide a TIME/YYYY variable.")
end

function _read_variable_data(ds::NCDataset, name::AbstractString)
    haskey(ds, name) || error("Variable '$name' not found in forcing file.")
    var = ds[name].var
    data = Float64.(var[ntuple(_ -> (:), ndims(var))...])
    return data, dimnames(var)
end

function _layer_dim_from_names(dim_names::Tuple)
    lower_names = lowercase.(String.(dim_names))
    return findfirst(name ->
        name == "z" ||
        name == "level" ||
        name == "lev" ||
        name == "lvl" ||
        occursin("height", name) ||
        occursin("level", name),
        lower_names,
    )
end

function _drop_singleton_layer(data, name::AbstractString)
    ndims(data) == 4 || return data
    singleton_dim = findfirst(==(1), size(data))
    if !isnothing(singleton_dim)
        return dropdims(data; dims=singleton_dim)
    end
    level_dim = findfirst(==(2), size(data))
    isnothing(level_dim) && error("`$name` is 4-D, but no singleton or 2-level layer dimension could be identified.")
    return Array(selectdim(data, level_dim, 1))
end

function _drop_singleton_layer(data, dim_names::Tuple, name::AbstractString)
    ndims(data) == 4 || return data, dim_names
    singleton_dim = findfirst(==(1), size(data))
    if !isnothing(singleton_dim)
        kept = ntuple(i -> i < singleton_dim ? dim_names[i] : dim_names[i + 1], ndims(data) - 1)
        return dropdims(data; dims=singleton_dim), kept
    end
    level_dim = _layer_dim_from_names(dim_names)
    if isnothing(level_dim) || size(data, level_dim) != 2
        level_dim = findfirst(==(2), size(data))
    end
    isnothing(level_dim) && error("`$name` is 4-D, but no singleton or 2-level layer dimension could be identified.")
    kept = ntuple(i -> i < level_dim ? dim_names[i] : dim_names[i + 1], ndims(data) - 1)
    return Array(selectdim(data, level_dim, 1)), kept
end

function _as_time_y_x(data, ntime::Int, ny::Int, nx::Int, name::AbstractString)
    A = _drop_singleton_layer(Float64.(data), name)
    nd = ndims(A)
    nd == 3 || error("`$name` must be a 3-D field, or a 4-D field with one singleton layer dimension.")
    dims = size(A)
    if dims == (ntime, ny, nx)
        return A
    elseif dims == (nx, ny, ntime)
        return permutedims(A, (3, 2, 1))
    elseif dims == (ny, nx, ntime)
        return permutedims(A, (3, 1, 2))
    elseif dims == (ntime, nx, ny)
        return permutedims(A, (1, 3, 2))
    end
    error("`$name` must have shape `(time, y, x)` or a simple permutation of it; got $(dims).")
end

function _as_time_y_x(data, dim_names::Tuple, ntime::Int, ny::Int, nx::Int, name::AbstractString)
    A, names = _drop_singleton_layer(Float64.(data), dim_names, name)
    if ndims(A) == 3
        lower_names = lowercase.(String.(names))
        tdim = findfirst(n -> n == "time" || n == "t" || occursin("time", n), lower_names)
        ydim = findfirst(==("y"), lower_names)
        xdim = findfirst(==("x"), lower_names)
        if !isnothing(tdim) && !isnothing(ydim) && !isnothing(xdim)
            dims = size(A)
            dims[tdim] == ntime || error("`$name` time dimension must have length $ntime, got $(dims[tdim]).")
            dims[ydim] == ny || error("`$name` y dimension must have length $ny, got $(dims[ydim]).")
            dims[xdim] == nx || error("`$name` x dimension must have length $nx, got $(dims[xdim]).")
            return permutedims(A, (tdim, ydim, xdim))
        end
    end
    return _as_time_y_x(A, ntime, ny, nx, name)
end

function _read_time_y_x(ds::NCDataset, name::AbstractString, ntime::Int, ny::Int, nx::Int)
    data, names = _read_variable_data(ds, name)
    return _as_time_y_x(data, names, ntime, ny, nx, name)
end

function _column_matrix(field::Array{Float64, 3})
    ntime, ny, nx = size(field)
    return reshape(permutedims(field, (3, 2, 1)), nx * ny, ntime)
end

function _as_y_x(data, dim_names::Tuple, ny::Int, nx::Int, name::AbstractString)
    A = Float64.(data)
    ndims(A) == 2 || error("`$name` must be a 2-D latitude/coordinate field, got shape $(size(A)).")
    lower_names = lowercase.(String.(dim_names))
    ydim = findfirst(==("y"), lower_names)
    xdim = findfirst(==("x"), lower_names)
    if !isnothing(ydim) && !isnothing(xdim)
        size(A, ydim) == ny || error("`$name` y dimension must have length $ny, got $(size(A, ydim)).")
        size(A, xdim) == nx || error("`$name` x dimension must have length $nx, got $(size(A, xdim)).")
        return permutedims(A, (ydim, xdim))
    elseif size(A) == (ny, nx)
        return A
    elseif size(A) == (nx, ny)
        return permutedims(A, (2, 1))
    end
    error("`$name` must have shape `(y, x)` or `(x, y)`; got $(size(A)).")
end

_column_vector_y_x(field::AbstractMatrix{<:Real}) = vec(permutedims(Float64.(field), (2, 1)))

"""
    read_forcing_columns(ds, name, rows, ntime, ny, nx; x_name="x", y_name="y")

Read the gridded variable `name` from the open dataset `ds` as a
`(length(rows), ntime)` `Float64` matrix in one pass. `rows` are flat column
indices `i + (j - 1) * nx`. The x and y dimensions are found by the names
`x_name`/`y_name` (or `x`/`y`) and the time dimension by name; any other
dimension (vertical level, sector) is read at its first index only. Variables
without recognisable x/y/time dimensions fall back to shape matching.
"""
function read_forcing_columns(ds::NCDataset, name::AbstractString, rows::AbstractVector{<:Integer},
        ntime::Int, ny::Int, nx::Int; x_name::AbstractString="x", y_name::AbstractString="y")
    haskey(ds, name) || error("Variable '$name' not found in forcing file.")
    var = ds[name].var
    names = lowercase.(String.(dimnames(var)))
    tdim = findfirst(n -> n == "time" || n == "t" || occursin("time", n), names)
    ydim = findfirst(n -> n == "y" || n == lowercase(y_name), names)
    xdim = findfirst(n -> n == "x" || n == lowercase(x_name), names)
    if isnothing(tdim) || isnothing(ydim) || isnothing(xdim)
        return _column_matrix(_read_time_y_x(ds, name, ntime, ny, nx))[rows, :]
    end
    sizes = size(var)
    (sizes[tdim], sizes[ydim], sizes[xdim]) == (ntime, ny, nx) ||
        error("`$name` must have time, y and x lengths ($ntime, $ny, $nx); got $(sizes).")
    raw = var[ntuple(d -> d in (tdim, ydim, xdim) ? Colon() : 1, ndims(var))...]
    kept = [d for d in 1:ndims(var) if d in (tdim, ydim, xdim)]
    order = (findfirst(==(xdim), kept), findfirst(==(ydim), kept), findfirst(==(tdim), kept))
    columns = reshape(order == (1, 2, 3) ? raw : permutedims(raw, order), nx * ny, ntime)
    return _gather_rows(columns, rows)
end

# Function barrier: the raw element type is only known at run time.
function _gather_rows(columns::AbstractMatrix, rows::AbstractVector{<:Integer})
    out = Matrix{Float64}(undef, length(rows), size(columns, 2))
    @inbounds for t in axes(columns, 2), (k, row) in enumerate(rows)
        out[k, t] = columns[row, t]
    end
    return out
end

"""Optional field as `(values, available)`: non-finite values are zero and unavailable.
An absent field is `(nothing, nothing)`, which `SnowpackForcing` stores as a constant."""
function _read_optional_columns(ds::NCDataset, name, rows, ntime::Int, ny::Int, nx::Int; kwargs...)
    (isnothing(name) || !haskey(ds, name)) && return nothing, nothing
    values = read_forcing_columns(ds, name, rows, ntime, ny, nx; kwargs...)
    available = isfinite.(values)
    values[.!available] .= 0.0
    return values, available
end

"""Static `(y, x)` field as a vector over the flat column indices `rows`."""
function _read_static_columns(ds::NCDataset, name::AbstractString, rows, ny::Int, nx::Int)
    data, dim_names = _read_variable_data(ds, name)
    return _column_vector_y_x(_as_y_x(data, dim_names, ny, nx, name))[rows]
end

"""
    load_forcing_file(path; kwargs...) -> (grid, forcing)

Load a generic NetCDF forcing file into a `SnowpackGrid` and
`SnowpackForcing`. Required variables default to `x`, `y`, `time`, `TT`, `SF`,
`RF`, and `SWD`; optional turbulent/radiative flux fields and prescribed
surface albedo are loaded when present and requested. If no explicit pressure
field is requested, the static `surface_height_name` field is used to derive
air pressure with the barometric formula. Set `mask_name` and `mask_threshold`
to load only columns selected by a spatial mask variable.
"""
function load_forcing_file(
    path::AbstractString;
    x_name::AbstractString="x",
    y_name::AbstractString="y",
    time_name::AbstractString="time",
    air_temperature_name::AbstractString="TT",
    snowfall_name::AbstractString="SF",
    rainfall_name::AbstractString="RF",
    shortwave_name::AbstractString="SWD",
    wind_speed_name::Union{Nothing, AbstractString}=nothing,
    q_lw_down_name::Union{Nothing, AbstractString}="LWD",
    q_sh_name::Union{Nothing, AbstractString}="SHF",
    q_lh_name::Union{Nothing, AbstractString}="LHF",
    relative_humidity_name::Union{Nothing, AbstractString}="RHZ",
    air_pressure_name::Union{Nothing, AbstractString}=nothing,
    surface_height_name::Union{Nothing, AbstractString}="SH",
    air_pressure_temperature_mode=:annual_mean,
    prescribed_albedo_name::Union{Nothing, AbstractString}=nothing,
    coszm_name::Union{Nothing, AbstractString}=nothing,
    cloud_name::Union{Nothing, AbstractString}=nothing,
    dust_deposition_name::Union{Nothing, AbstractString}=nothing,
    z_sur_std_name::Union{Nothing, AbstractString}=nothing,
    prescribed_ice_albedo_name::Union{Nothing, AbstractString}=nothing,
    latitude_name::Union{Nothing, AbstractString}="LAT",
    mask_name::Union{Nothing, AbstractString}=nothing,
    mask_threshold::Real=0.0,
    air_temperature_in_celsius::Bool=true,
    precipitation_in_mmwe_day::Bool=true,
    air_pressure_default::Float64=DEFAULT_SEA_LEVEL_AIR_PRESSURE,
    wind_default::Float64=5.0,
)
    ds = NCDataset(path)
    try
        x = Float64.(vec(ds[x_name][:]))
        y = Float64.(vec(ds[y_name][:]))
        nx, ny = length(x), length(y)

        # Select the columns first so that every field is read straight into them.
        mask = ones(Float64, ny, nx)
        rows = collect(1:(nx * ny))
        if !isnothing(mask_name)
            mask_data, mask_dims = _read_variable_data(ds, mask_name)
            mask = _as_y_x(mask_data, mask_dims, ny, nx, mask_name)
            mask_columns = _column_vector_y_x(mask)
            rows = findall(isfinite.(mask_columns) .& (mask_columns .>= mask_threshold))
            isempty(rows) && error("`$mask_name` selected no columns.")
        end
        ncol = length(rows)

        ntime = if any(name -> haskey(ds, name), (time_name, "TIME", "YYYY"))
            _infer_time_length(ds, time_name, nothing, ny, nx, air_temperature_name)
        else
            _infer_time_length(ds, time_name, first(_read_variable_data(ds, air_temperature_name)), ny, nx, air_temperature_name)
        end
        time_values = _read_time_values(ds, time_name, ntime)
        dt_days = infer_dt_days(time_values)
        columns(name) = read_forcing_columns(ds, name, rows, ntime, ny, nx; x_name, y_name)
        optional(name) = _read_optional_columns(ds, name, rows, ntime, ny, nx; x_name, y_name)

        air_temperature = columns(air_temperature_name)
        air_temperature_in_celsius && (air_temperature .+= 273.15)
        snowfall_rate = columns(snowfall_name)
        rainfall_rate = columns(rainfall_name)
        if precipitation_in_mmwe_day
            snowfall_rate ./= 86_400.0
            rainfall_rate ./= 86_400.0
        end
        shortwave_down = columns(shortwave_name)
        wind_speed = isnothing(wind_speed_name) || !haskey(ds, wind_speed_name) ?
            fill(wind_default, ncol, ntime) : columns(wind_speed_name)
        latitude_deg = !isnothing(latitude_name) && haskey(ds, latitude_name) ?
            _read_static_columns(ds, latitude_name, rows, ny, nx) : nothing

        q_lw_m, has_q_lw_m = optional(q_lw_down_name)
        q_sh_m, has_q_sh_m = optional(q_sh_name)
        q_lh_m, has_q_lh_m = optional(q_lh_name)
        relative_humidity_m, has_relative_humidity_m = optional(relative_humidity_name)

        has_surface_height = !isnothing(surface_height_name) && haskey(ds, surface_height_name)
        surface_height_m = fill(NaN, ncol, ntime)
        has_surface_height && (surface_height_m .= _read_static_columns(ds, surface_height_name, rows, ny, nx))
        air_pressure_m = fill(air_pressure_default, ncol, ntime)
        if !isnothing(air_pressure_name) && haskey(ds, air_pressure_name)
            pressure_raw = columns(air_pressure_name)
            finite_pressure = isfinite.(pressure_raw)
            air_pressure_m[finite_pressure] .= pressure_raw[finite_pressure]
        elseif has_surface_height
            air_pressure_m .= air_pressure_from_surface_height(
                surface_height_m,
                air_temperature;
                dt_days=dt_days,
                time_values=time_values,
                temperature_mode=air_pressure_temperature_mode,
                sea_level_pressure=air_pressure_default,
            )
        end

        prescribed_albedo_m, has_prescribed_albedo_m = optional(prescribed_albedo_name)
        coszm_m, has_coszm_m = optional(coszm_name)
        cloud_m, has_cloud_m = optional(cloud_name)
        dust_deposition_m, has_dust_deposition_m = optional(dust_deposition_name)
        z_sur_std_m, has_z_sur_std_m = optional(z_sur_std_name)
        prescribed_ice_albedo_m, has_prescribed_ice_albedo_m = optional(prescribed_ice_albedo_name)

        js = repeat(collect(1:ny), inner=nx)[rows]
        is = repeat(collect(1:nx), outer=ny)[rows]
        grid = SnowpackGrid(ncol; x=x, y=y, js=js, is=is, mask=mask)
        forcing = SnowpackForcing(
            time_values=time_values,
            dt_days=dt_days,
            air_temperature=air_temperature,
            snowfall_rate=snowfall_rate,
            rainfall_rate=rainfall_rate,
            shortwave_down=shortwave_down,
            wind_speed=wind_speed,
            q_lw_down=q_lw_m,
            has_q_lw_down=has_q_lw_m,
            q_sh=q_sh_m,
            has_q_sh=has_q_sh_m,
            q_lh=q_lh_m,
            has_q_lh=has_q_lh_m,
            relative_humidity=relative_humidity_m,
            has_relative_humidity=has_relative_humidity_m,
            surface_height=surface_height_m,
            air_pressure=air_pressure_m,
            prescribed_albedo=prescribed_albedo_m,
            has_prescribed_albedo=has_prescribed_albedo_m,
            coszm=coszm_m, has_coszm=has_coszm_m,
            cloud=cloud_m, has_cloud=has_cloud_m,
            dust_deposition=dust_deposition_m, has_dust_deposition=has_dust_deposition_m,
            z_sur_std=z_sur_std_m, has_z_sur_std=has_z_sur_std_m,
            prescribed_ice_albedo=prescribed_ice_albedo_m, has_prescribed_ice_albedo=has_prescribed_ice_albedo_m,
            latitude_deg=latitude_deg,
        )
        return (grid=grid, forcing=forcing)
    finally
        close(ds)
    end
end
