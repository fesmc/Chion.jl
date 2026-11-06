"""
`Forcing`: where an `Experiment` takes its forcing from, one calendar year at a time.

    Forcing("/data/MAR/MARv3.14.3-10km-daily-ERA5-{year}.nc"; format=:mar, surface_fluxes=:parameterized)
    Forcing("/data/site.nc"; variables=(air_temperature="Tair", ...))     # one file, any length
    Forcing(snowpack_forcing; grid)                                         # in memory
    Forcing(year -> my_forcing(year), 1979:2025; grid, area)                # user loader

Only the year that is being run is held in memory, so long multi-file records fit on a GPU node.
"""

abstract type AbstractForcingSource end

"""Files with one year each (`{year}` in the path pattern)."""
struct YearFileSource <: AbstractForcingSource
    files::Dict{Int, String}
    loader::Any                  # path -> (grid, forcing)
end

"""One file (or in-memory forcing) covering any period; years are sliced from it."""
mutable struct WholeSource <: AbstractForcingSource
    loader::Any                  # () -> (grid, forcing)
    cache::Any
end

"""A user function `year -> SnowpackForcing`."""
struct FunctionSource <: AbstractForcingSource
    f::Any
    years::Vector{Int}
end

struct Forcing{S <: AbstractForcingSource}
    source::S
    surface_fluxes::Symbol       # :prescribed, :parameterized or :as_loaded
    transform::Any               # nothing or (forcing, year) -> forcing
    description::String
    grid_ref::Base.RefValue{Any}
    area_ref::Base.RefValue{Any}
    area_loader::Any             # nothing or grid -> Vector{Float64}
end

"""
    Forcing(path; format=:auto, surface_fluxes=:prescribed, years=nothing, files=Dict(),
            mask=nothing, wind=nothing, area=:auto, transform=nothing, variables=nothing, kwargs...)

- `path`: a file, or a pattern with `{year}` for one file per year (missing years are skipped).
- `format`: `:mar` (MAR daily output: variable names, the ice mask `MSK` ≥ 50 %, a constant wind
  of 5 m s⁻¹ and the cell area `AREA`), `:alma` (single-site forcing in the ALMA convention of
  land-surface model intercomparisons such as ESM-SnowMIP: `Tair`, `Qair`, `Psurf`, `Wind`,
  `Rainf`, `Snowf`, `SWdown`, `LWdown` in SI units) or `:generic` (names from `variables` and
  `kwargs`, as in `load_forcing_file`). `:auto` picks `:mar` when the file has `MSK` and
  `TTZ`/`TT`, and `:alma` when it has `Tair` and `Qair`.
- `surface_fluxes`: `:prescribed` uses the surface energy fluxes of the file (for MAR the net
  shortwave from SWD − SWU, LWD, SHF, LHF and the albedo AL2), `:parameterized` drops them so
  that Chion computes them, `:albedo` keeps only the albedo (Chion computes the energy fluxes
  with it; use with `albedo=:prescribed` in the model), `:as_loaded` keeps whatever the loader
  provides.
- `files`: `Dict(year => path)` overrides single years (e.g. a repaired last year).
- `mask`: `(name, threshold)` selects columns, e.g. `("MSK", 50.0)`.
- `wind`: constant wind speed (m s⁻¹) or the name of a wind variable.
- `area`: column areas (km²): `:auto` (MAR `AREA`), a variable name, a vector, or `nothing`.
- `transform`: `(forcing, year) -> forcing`, applied to every year (e.g. a temperature change).
- `variables`: NamedTuple of variable names for `:generic` files, e.g.
  `(air_temperature="t2m", snowfall="sf", rainfall="rf", shortwave_down="ssrd")`.
"""
function Forcing(
    path::AbstractString;
    format::Symbol=:auto,
    surface_fluxes::Symbol=:prescribed,
    years=nothing,
    files=Dict{Int, String}(),
    mask=nothing,
    wind=nothing,
    area=:auto,
    area_file=nothing,
    transform=nothing,
    variables=nothing,
    kwargs...,
)
    surface_fluxes in (:prescribed, :parameterized, :albedo, :as_loaded) ||
        error("`surface_fluxes` must be :prescribed, :parameterized, :albedo or :as_loaded.")
    yearly = occursin("{year}", path)
    file_years = Dict{Int, String}()
    if yearly
        dir, pat = dirname(path), basename(path)
        rx = Regex("^" * replace(replace(pat, "." => "\\."), "{year}" => "(\\d{4})") * "\$")
        isdir(dir) || error("Forcing directory not found: $dir")
        for f in readdir(dir)
            m = match(rx, f)
            m === nothing || (file_years[parse(Int, m.captures[1])] = joinpath(dir, f))
        end
        for (y, p) in files
            file_years[Int(y)] = String(p)
        end
        years === nothing || filter!(kv -> first(kv) in years, file_years)
        isempty(file_years) && error("No forcing files match $path")
    else
        isfile(path) || error("Forcing file not found: $path")
    end
    probe = yearly ? file_years[minimum(keys(file_years))] : path
    fmt = format == :auto ? _detect_format(probe) : format
    fmt in (:mar, :alma, :generic) || error("`format` must be :auto, :mar, :alma or :generic.")
    loader = p -> _load_file(p, fmt, surface_fluxes; mask, wind, variables, kwargs...)
    area_loader = _area_loader(area, fmt, probe, area_file; kwargs...)
    desc = "$(yearly ? "yearly files " : "file ")$path (format $fmt, surface fluxes $surface_fluxes)"
    source = yearly ? YearFileSource(file_years, loader) : WholeSource(() -> loader(path), nothing)
    return Forcing(source, surface_fluxes, transform, desc, Ref{Any}(nothing), Ref{Any}(nothing), area_loader)
end

"""
    Forcing(forcing::SnowpackForcing; grid=SnowpackGrid(ncol), area=nothing, transform=nothing)

Wrap forcing that is already in memory. Years are sliced from its time axis.
"""
function Forcing(forcing::SnowpackForcing; grid=nothing, area=nothing, transform=nothing, description="in-memory forcing")
    ncol = size(forcing.air_temperature, 1)
    g = grid === nothing ? SnowpackGrid(ncol) : grid
    a = area === nothing ? nothing : Float64.(collect(area))
    return Forcing(WholeSource(() -> (grid=g, forcing=forcing), nothing), :as_loaded, transform, description,
        Ref{Any}(g), Ref{Any}(a), nothing)
end

"""
    Forcing(f, years; grid, area=nothing, transform=nothing)

Forcing from a user function `f(year) -> SnowpackForcing` (for preprocessing that the file
loaders do not cover). `grid` is the `SnowpackGrid` of the columns.
"""
function Forcing(f::Function, years; grid, area=nothing, transform=nothing, description="user forcing function")
    a = area === nothing ? nothing : Float64.(collect(area))
    return Forcing(FunctionSource(f, collect(Int, years)), :as_loaded, transform, description, Ref{Any}(grid), Ref{Any}(a), nothing)
end

# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

function _detect_format(path)
    NCDataset(path) do ds
        haskey(ds, "MSK") && (haskey(ds, "TTZ") || haskey(ds, "TT")) && return :mar
        haskey(ds, "Tair") && haskey(ds, "Qair") && return :alma
        return :generic
    end
end

"""Copy of `values` with every non-finite (or negative, if `nonnegative`) value replaced by the last
valid one, or `fallback` before the first."""
function _forward_fill(values, fallback; nonnegative=false)
    out = copy(values)
    previous = fallback
    for i in eachindex(out)
        valid = isfinite(out[i]) && (!nonnegative || out[i] >= 0)
        valid ? (previous = out[i]) : (out[i] = previous)
    end
    return out
end

"""
Single-site ALMA forcing (`Tair` K, `Qair` kg kg⁻¹, `Psurf` Pa, `Wind` m s⁻¹, `Rainf` and `Snowf`
kg m⁻² s⁻¹, `SWdown` and `LWdown` W m⁻²). Gaps in temperature, wind, pressure and humidity are
filled with the last valid value; missing precipitation and shortwave are zero; longwave is
marked unavailable where missing (Chion then computes it). Specific humidity becomes relative
humidity with respect to ice.
"""
function _load_alma(path)
    met = NCDataset(path) do ds
        v(name) = (haskey(ds, name) || error("`$name` is missing from $(basename(path)) (ALMA forcing).");
                   Float64.(coalesce.(vec(ds[name][:]), NaN)))
        (dates=DateTime.(ds["time"][:]), tair=v("Tair"), snow=max.(v("Snowf"), 0.0), rain=max.(v("Rainf"), 0.0),
         sw=max.(v("SWdown"), 0.0), lw=v("LWdown"), wind=max.(v("Wind"), 0.0), qair=v("Qair"), psurf=v("Psurf"))
    end
    n = length(met.dates)
    n > 1 || error("ALMA forcing needs at least two time records.")
    dt_days = [Dates.value(met.dates[i + 1] - met.dates[i]) / 86_400_000 for i in 1:(n - 1)]
    push!(dt_days, dt_days[end])
    air_temperature = _forward_fill(met.tair, 268.0)
    air_pressure = _forward_fill(met.psurf, _median(filter(isfinite, met.psurf)); nonnegative=true)
    # q_air = 0.622 RH e_s,ice(T_air) / p
    e_ice = _bessi_ice_saturation_vapor_pressure.(air_temperature, 273.15)
    relative_humidity = clamp.(_forward_fill(met.qair, 0.0; nonnegative=true) .* air_pressure ./ (0.622 .* e_ice), 0.0, 1.0)
    forcing = SnowpackForcing(dt_days=dt_days, ncol=1, air_temperature=air_temperature,
        snowfall_rate=replace(met.snow, NaN => 0.0), rainfall_rate=replace(met.rain, NaN => 0.0),
        shortwave_down=replace(met.sw, NaN => 0.0), q_lw_down=replace(met.lw, NaN => 0.0),
        has_q_lw_down=isfinite.(met.lw), wind_speed=_forward_fill(met.wind, 5.0; nonnegative=true),
        relative_humidity=relative_humidity, has_relative_humidity=trues(n), air_pressure=air_pressure,
        time_values=met.dates)
    return (grid=SnowpackGrid(1), forcing=forcing)
end

function _median(v)
    isempty(v) && return DEFAULT_SEA_LEVEL_AIR_PRESSURE
    s = sort(v); n = length(s)
    return isodd(n) ? s[(n + 1) ÷ 2] : s[n ÷ 2] / 2 + s[n ÷ 2 + 1] / 2      # as Statistics.median
end

"""x and y dimension names of a MAR file (from the mask), e.g. `x`/`y` or `X18_215`/`Y15_176`."""
function _mar_xy_names(ds)
    names = dimnames(ds["MSK"])
    x = findfirst(n -> lowercase(n) == "x" || startswith(uppercase(n), "X"), names)
    y = findfirst(n -> lowercase(n) == "y" || startswith(uppercase(n), "Y"), names)
    (x === nothing || y === nothing) && error("Cannot identify the x/y dimensions of MSK: $(names)")
    return String(names[x]), String(names[y])
end

const _GENERIC_NAME_KEYS = (air_temperature=:air_temperature_name, snowfall=:snowfall_name, rainfall=:rainfall_name,
    shortwave_down=:shortwave_name, wind_speed=:wind_speed_name, longwave_down=:q_lw_down_name,
    sensible_heat_flux=:q_sh_name, latent_heat_flux=:q_lh_name, relative_humidity=:relative_humidity_name,
    air_pressure=:air_pressure_name, surface_height=:surface_height_name, albedo=:prescribed_albedo_name,
    latitude=:latitude_name, time=:time_name, x=:x_name, y=:y_name)

function _load_file(path, fmt, surface_fluxes; mask, wind, variables, kwargs...)
    if fmt == :alma
        loaded = _load_alma(path)
        surface_fluxes == :parameterized && return (grid=loaded.grid, forcing=_without_surface_fluxes(loaded.forcing))
        return loaded
    end
    if fmt == :mar
        xn, yn, tname, has_swu = NCDataset(path) do ds
            (_mar_xy_names(ds)..., haskey(ds, "TTZ") ? "TTZ" : "TT", haskey(ds, "SWU"))
        end
        mname, mthr = mask === nothing ? ("MSK", 50.0) : mask
        wind_name = wind isa AbstractString ? String(wind) : nothing
        wind_default = wind isa Real ? Float64(wind) : 5.0
        loaded = load_forcing_file(path; x_name=xn, y_name=yn, time_name="TIME", air_temperature_name=tname,
            wind_speed_name=wind_name, wind_default, mask_name=mname, mask_threshold=mthr, kwargs...)
        surface_fluxes == :prescribed && return (grid=loaded.grid, forcing=_mar_prescribed(path, loaded, has_swu, xn, yn))
        surface_fluxes == :albedo && return (grid=loaded.grid,
            forcing=_without_surface_fluxes(_mar_prescribed(path, loaded, has_swu, xn, yn); keep_albedo=true))
        surface_fluxes == :parameterized && return (grid=loaded.grid, forcing=_without_surface_fluxes(loaded.forcing))
        surface_fluxes == :albedo && return (grid=loaded.grid, forcing=_without_surface_fluxes(loaded.forcing; keep_albedo=true))
        return loaded
    end
    names = Dict{Symbol, Any}()
    if variables !== nothing
        for (k, v) in pairs(variables)
            haskey(_GENERIC_NAME_KEYS, k) || error("Unknown forcing variable `$k`; use $(join(string.(keys(_GENERIC_NAME_KEYS)), ", ")).")
            names[getfield(_GENERIC_NAME_KEYS, k)] = v
        end
    end
    if mask !== nothing
        names[:mask_name], names[:mask_threshold] = mask
    end
    wind isa AbstractString && (names[:wind_speed_name] = wind)
    wind isa Real && (names[:wind_default] = Float64(wind))
    loaded = load_forcing_file(path; names..., kwargs...)
    surface_fluxes == :parameterized && return (grid=loaded.grid, forcing=_without_surface_fluxes(loaded.forcing))
    surface_fluxes == :albedo && return (grid=loaded.grid, forcing=_without_surface_fluxes(loaded.forcing; keep_albedo=true))
    return loaded
end

"""MAR with every surface-energy flux and the albedo prescribed."""
function _mar_prescribed(path, loaded, has_swu, xn, yn)
    grid, forcing = loaded.grid, loaded.forcing
    ntime = length(forcing.time_values)
    rows = grid.is .+ (grid.js .- 1) .* length(grid.x)
    read(name) = NCDataset(ds -> read_forcing_columns(ds, name, rows, ntime, length(grid.y), length(grid.x);
        x_name=xn, y_name=yn), path)
    albedo = clamp.(read("AL2"), 0.0, 1.0)
    # Recent MAR files (e.g. Antarctica from 2025) omit SWU; derive it from AL2.
    shortwave_up = has_swu ? read("SWU") : forcing.shortwave_down .* albedo
    return SnowpackForcing(time_values=forcing.time_values, dt_days=forcing.dt_days,
        air_temperature=forcing.air_temperature, snowfall_rate=forcing.snowfall_rate, rainfall_rate=forcing.rainfall_rate,
        shortwave_down=forcing.shortwave_down, wind_speed=forcing.wind_speed,
        q_sw_net=forcing.shortwave_down .- shortwave_up, has_q_sw_net=true,
        q_lw_down=forcing.q_lw_down, has_q_lw_down=true, q_sh=forcing.q_sh, has_q_sh=true,
        q_lh=forcing.q_lh, has_q_lh=true,
        relative_humidity=forcing.relative_humidity, has_relative_humidity=forcing.has_relative_humidity,
        air_pressure=forcing.air_pressure, surface_height=forcing.surface_height,
        prescribed_albedo=albedo, has_prescribed_albedo=true, latitude_deg=forcing.latitude_deg)
end

"""Meteorology and precipitation only: every surface-energy flux (and, unless `keep_albedo`, the
albedo) is computed by Chion."""
function _without_surface_fluxes(forcing; keep_albedo::Bool=false)
    return SnowpackForcing(time_values=forcing.time_values, dt_days=forcing.dt_days,
        air_temperature=forcing.air_temperature, snowfall_rate=forcing.snowfall_rate, rainfall_rate=forcing.rainfall_rate,
        shortwave_down=forcing.shortwave_down, wind_speed=forcing.wind_speed,
        q_sw_net=forcing.q_sw_net, has_q_sw_net=false, q_lw_down=forcing.q_lw_down, has_q_lw_down=false,
        q_sh=forcing.q_sh, has_q_sh=false, q_lh=forcing.q_lh, has_q_lh=false,
        relative_humidity=forcing.relative_humidity, has_relative_humidity=forcing.has_relative_humidity,
        air_pressure=forcing.air_pressure, surface_height=forcing.surface_height,
        prescribed_albedo=forcing.prescribed_albedo, has_prescribed_albedo=keep_albedo && forcing.has_prescribed_albedo,
        latitude_deg=forcing.latitude_deg)
end

function _area_loader(area, fmt, probe, area_file; kwargs...)
    area === nothing && return nothing
    area isa AbstractVector && (v = Float64.(area); return _ -> v)
    name = area === :auto ? (fmt == :mar ? "AREA" : nothing) : String(area)
    name === nothing && return nothing
    return function (grid)
        for file in (probe, area_file)
            file === nothing && continue
            has = NCDataset(ds -> haskey(ds, name), file)
            has || continue
            field = NCDataset(file) do ds
                raw, dims = _read_variable_data(ds, name)
                _as_y_x(raw, dims, length(grid.y), length(grid.x), name)
            end
            return [Float64(field[grid.js[c], grid.is[c]]) for c in eachindex(grid.js)]
        end
        @warn "No `$name` in the forcing; area-aggregated output needs `area` in the Output."
        return nothing
    end
end

# ---------------------------------------------------------------------------
# Access by year
# ---------------------------------------------------------------------------

"""Years the forcing can provide."""
forcing_years(f::Forcing{YearFileSource}) = sort(collect(keys(f.source.files)))
forcing_years(f::Forcing{FunctionSource}) = f.source.years
forcing_years(f::Forcing{WholeSource}) = sort(unique(year.(_whole(f).forcing.time_values)))

function _whole(f::Forcing{WholeSource})
    f.source.cache === nothing && (f.source.cache = f.source.loader())
    return f.source.cache
end

_slice_field(x::ConstantForcingMatrix, r, nt) = ConstantForcingMatrix(x.value, (size(x, 1), length(r)))
_slice_field(x::ColumnForcingMatrix, r, nt) = ColumnForcingMatrix(copy(x.values), length(r))
_slice_field(x::TimeForcingMatrix, r, nt) = TimeForcingMatrix(x.values[r], x.ncol)
_slice_field(x, r, nt) = x isa OptionalForcingField ? OptionalForcingField(_slice_field(x.values, r, nt), _slice_field(x.available, r, nt)) :
    x isa AbstractMatrix && size(x, 2) == nt ? x[:, r] : x isa AbstractVector && length(x) == nt ? x[r] : x

"""The forcing restricted to the time indices `r`."""
function slice_forcing(f::SnowpackForcing, r)
    nt = length(f.calendar.time_values)
    cal = f.calendar
    calendar = typeof(cal).name.wrapper((_slice_field(getfield(cal, k), r, nt) for k in fieldnames(typeof(cal)))...)
    return SnowpackForcing(calendar, map(x -> _slice_field(x, r, nt), f.fields))
end

"""
    relabel_year(forcing, year)

`forcing` with its time stamps moved to `year` (29 February becomes 28 February), for running the
weather of one year as another, e.g. base years under a projected warming. The day of year and the
solar geometry stay those of the original year.
"""
function relabel_year(f::SnowpackForcing, y::Integer)
    cal = f.calendar
    shift = Year(y - year(first(cal.time_values)))
    calendar = typeof(cal).name.wrapper(cal.time_values .+ shift, cal.dt_days, cal.day_of_year, cal.solar_longitude_deg)
    return SnowpackForcing(calendar, f.fields)
end

"""`(grid, forcing)` of one calendar year (before `transform`)."""
function _load_year(f::Forcing{YearFileSource}, y::Int)
    haskey(f.source.files, y) || error("No forcing file for $y.")
    return f.source.loader(f.source.files[y])
end
function _load_year(f::Forcing{WholeSource}, y::Int)
    w = _whole(f)
    r = findall(t -> year(t) == y, w.forcing.time_values)
    isempty(r) && error("The forcing has no time steps in $y.")
    r == 1:length(w.forcing.time_values) && return w
    return (grid=w.grid, forcing=slice_forcing(w.forcing, first(r):last(r)))
end
function _load_year(f::Forcing{FunctionSource}, y::Int)
    forcing = f.source.f(y)
    return (grid=f.grid_ref[], forcing=forcing)
end

"""Forcing of one year, checked against the grid of the first year and transformed."""
function load_forcing_year(f::Forcing, y::Integer)
    loaded = _load_year(f, Int(y))
    g = forcing_grid(f; loaded)
    _same_columns(g, loaded.grid) || error("The columns of $y differ from those of the first forcing year.")
    forcing = f.transform === nothing ? loaded.forcing : f.transform(loaded.forcing, Int(y))
    return forcing
end

_same_columns(a::SnowpackGrid, b::SnowpackGrid) = ncols(a) == ncols(b) &&
    (!has_spatial_coords(a) || !has_spatial_coords(b) || (a.x == b.x && a.y == b.y && a.js == b.js && a.is == b.is))

"""The `SnowpackGrid` of the forcing columns."""
function forcing_grid(f::Forcing; loaded=nothing)
    if f.grid_ref[] === nothing
        l = loaded === nothing ? _load_year(f, first(forcing_years(f))) : loaded
        f.grid_ref[] = l.grid
    end
    return f.grid_ref[]
end

"""Column areas in km², or `nothing`."""
function forcing_area(f::Forcing)
    if f.area_ref[] === nothing && f.area_loader !== nothing
        f.area_ref[] = f.area_loader(forcing_grid(f))
    end
    return f.area_ref[]
end

Base.show(io::IO, f::Forcing) = print(io, "Forcing(", f.description, ")")
