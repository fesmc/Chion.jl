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
        raw = ds[time_name].var[:]
        if eltype(raw) <: DateTime
            return DateTime.(raw)
        end
    end
    if haskey(ds, "TIME")
        raw = ds["TIME"].var[:]
        if eltype(raw) <: DateTime
            return DateTime.(raw)
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
    base = DateTime(2000, 1, 1, 12)
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

function _drop_singleton_layer(data, name::AbstractString)
    ndims(data) == 4 || return data
    singleton_dim = findfirst(==(1), size(data))
    isnothing(singleton_dim) && error("`$name` is 4-D, but no singleton layer dimension could be identified.")
    return dropdims(data; dims=singleton_dim)
end

function _drop_singleton_layer(data, dim_names::Tuple, name::AbstractString)
    ndims(data) == 4 || return data, dim_names
    singleton_dim = findfirst(==(1), size(data))
    isnothing(singleton_dim) && error("`$name` is 4-D, but no singleton layer dimension could be identified.")
    kept = ntuple(i -> i < singleton_dim ? dim_names[i] : dim_names[i + 1], ndims(data) - 1)
    return dropdims(data; dims=singleton_dim), kept
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

function _wind_matrix(ds::NCDataset, wind_speed_name, ntime::Int, ny::Int, nx::Int, wind_default::Float64)
    isnothing(wind_speed_name) && return fill(wind_default, nx * ny, ntime)
    haskey(ds, wind_speed_name) || return fill(wind_default, nx * ny, ntime)
    raw = _read_time_y_x(ds, wind_speed_name, ntime, ny, nx)
    return _column_matrix(raw)
end

"""
    load_forcing_file(path; kwargs...) -> (grid, forcing)

Load a generic NetCDF forcing file into a `SnowpackGrid` and
`SnowpackForcing`. Required variables default to `x`, `y`, `time`, `TT`, `SF`,
`RF`, and `SWD`; optional turbulent/radiative flux fields are loaded when
present.
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
    air_temperature_in_celsius::Bool=true,
    precipitation_in_mmwe_day::Bool=true,
    wind_default::Float64=5.0,
    device=CPU(),
)
    ds = NCDataset(path)
    try
        x = Float64.(vec(ds[x_name][:]))
        y = Float64.(vec(ds[y_name][:]))
        nx, ny = length(x), length(y)

        tair_raw, tair_names = _read_variable_data(ds, air_temperature_name)
        ntime = _infer_time_length(ds, time_name, tair_raw, ny, nx, air_temperature_name)
        tair = _as_time_y_x(tair_raw, tair_names, ntime, ny, nx, air_temperature_name)
        snow = _read_time_y_x(ds, snowfall_name, ntime, ny, nx)
        rain = _read_time_y_x(ds, rainfall_name, ntime, ny, nx)
        sw = _read_time_y_x(ds, shortwave_name, ntime, ny, nx)

        time_values = _read_time_values(ds, time_name, ntime)
        dt_days = infer_dt_days(time_values)

        tair_m = _column_matrix(tair)
        snow_m = _column_matrix(snow)
        rain_m = _column_matrix(rain)
        sw_m = _column_matrix(sw)
        wind_m = _wind_matrix(ds, wind_speed_name, ntime, ny, nx, wind_default)

        q_lw_m = zeros(Float64, nx * ny, ntime)
        has_q_lw_m = fill(false, nx * ny, ntime)
        if !isnothing(q_lw_down_name) && haskey(ds, q_lw_down_name)
            q_lw_raw = _column_matrix(_read_time_y_x(ds, q_lw_down_name, ntime, ny, nx))
            q_lw_m .= q_lw_raw
            has_q_lw_m .= isfinite.(q_lw_raw)
            q_lw_m[.!has_q_lw_m] .= 0.0
        end

        q_sh_m = zeros(Float64, nx * ny, ntime)
        has_q_sh_m = fill(false, nx * ny, ntime)
        if !isnothing(q_sh_name) && haskey(ds, q_sh_name)
            q_sh_raw = _column_matrix(_read_time_y_x(ds, q_sh_name, ntime, ny, nx))
            q_sh_m .= q_sh_raw
            has_q_sh_m .= isfinite.(q_sh_raw)
            q_sh_m[.!has_q_sh_m] .= 0.0
        end

        q_lh_m = zeros(Float64, nx * ny, ntime)
        has_q_lh_m = fill(false, nx * ny, ntime)
        if !isnothing(q_lh_name) && haskey(ds, q_lh_name)
            q_lh_raw = _column_matrix(_read_time_y_x(ds, q_lh_name, ntime, ny, nx))
            q_lh_m .= q_lh_raw
            has_q_lh_m .= isfinite.(q_lh_raw)
            q_lh_m[.!has_q_lh_m] .= 0.0
        end

        air_temperature = air_temperature_in_celsius ? tair_m .+ 273.15 : tair_m
        snowfall_rate = precipitation_in_mmwe_day ? snow_m ./ 86_400.0 : snow_m
        rainfall_rate = precipitation_in_mmwe_day ? rain_m ./ 86_400.0 : rain_m

        js = repeat(collect(1:ny), inner=nx)
        is = repeat(collect(1:nx), outer=ny)
        grid = SnowpackGrid(device, nx * ny; x=x, y=y, js=js, is=is, mask=ones(ny, nx))
        forcing = SnowpackForcing(
            time_values=time_values,
            dt_days=dt_days,
            air_temperature=air_temperature,
            snowfall_rate=snowfall_rate,
            rainfall_rate=rainfall_rate,
            shortwave_down=sw_m,
            wind_speed=wind_m,
            q_lw_down=q_lw_m,
            has_q_lw_down=has_q_lw_m,
            q_sh=q_sh_m,
            has_q_sh=has_q_sh_m,
            q_lh=q_lh_m,
            has_q_lh=has_q_lh_m,
        )
        return (grid=grid, forcing=forcing)
    finally
        close(ds)
    end
end
