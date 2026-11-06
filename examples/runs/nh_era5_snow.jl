#!/usr/bin/env julia

# Northern Hemisphere seasonal snow forced by daily 1° ERA5 (preprocess_era5_nh_daily.sh).
# Land cells 25.5–89.5°N (land-sea mask ≥ 0.5) without the Greenland ice sheet, the SnowMIP
# land configuration (no glacier-ice substrate, HTESSEL densification, parameterised fresh-snow
# density), ERA5 longwave down, 10 m wind, 2 m dewpoint (as relative humidity) and surface
# pressure, daily-to-diurnal substeps of the defaults.
# One spin-up year (CHION_NH_SPINUP_YEAR, default 1979), then CHION_NH_FIRST..CHION_NH_LAST.
# Output (output/nh_era5_snow):
#   daily_snow.csv   NH / Eurasia / North America snow-covered area (km², SWE ≥ 5 kg m⁻²) and snow
#                    mass (Gt) per day
#   monthly_snow.nc  monthly mean SWE and snow-cover frequency per 1° cell (x = longitude, y = latitude)
#
#   julia --project=. --threads 16 examples/runs/nh_era5_snow.jl
using Chion, NCDatasets, Dates, Statistics, Printf

const ERA = get(ENV, "CHION_ERA5_NH_DIR", "/p/projects/ou/labs/ai/Nils/ERA5_nh_daily_1deg")
const NOUT = get(ENV, "CHION_NH_OUTPUT_DIR", joinpath(@__DIR__, "..", "..", "output", "nh_era5_snow"))
const SPINUP_YEAR = parse(Int, get(ENV, "CHION_NH_SPINUP_YEAR", "1979"))
const YEARS = parse(Int, get(ENV, "CHION_NH_FIRST", "1979")):parse(Int, get(ENV, "CHION_NH_LAST", "2025"))
const SNOW_THRESHOLD = 5.0                       # kg m⁻² SWE for a snow-covered cell
const EARTH_RADIUS_KM = 6371.0

"""Greenland, excluded as in the Rutgers North America series without Greenland (coarse box
mask). East of 335°E only north of 67°N, so that Iceland (63–67°N, 336–347°E) stays in."""
is_greenland(lon, lat) = (lat >= 59.5 && 300 <= lon <= 335) || (lat >= 67 && 300 <= lon <= 350) ||
                         (lat >= 76 && 288 <= lon <= 300)
is_north_america(lon) = 190 <= lon < 310

function domain()
    lon, lat, lsm, z = NCDataset(joinpath(ERA, "static.nc")) do ds
        # land-sea mask (GRIB parameter 172) and surface geopotential
        Float64.(ds["lon"][:]), Float64.(ds["lat"][:]), Float64.(coalesce.(ds["var172"][:, :, 1], 0)),
        Float64.(coalesce.(ds["z"][:, :, 1], 0))
    end
    cells = [(i, j) for j in eachindex(lat), i in eachindex(lon)
             if lsm[i, j] >= 0.5 && !is_greenland(lon[i], lat[j])]
    is, js = first.(cells), last.(cells)
    clat = lat[js]; clon = lon[is]
    area = [deg2rad(1)^2 * EARTH_RADIUS_KM^2 * cosd(clat[k]) * lsm[is[k], js[k]] for k in eachindex(is)]
    return (; lon, lat, is, js, clat, clon, area, height=[z[is[k], js[k]] / 9.80665 for k in eachindex(is)])
end

"""(column, day) matrix of one variable and year."""
function read_var(name, year, d)
    NCDataset(joinpath(ERA, "$(name)_$(year).nc")) do ds
        raw = ds[name][:, :, :]                        # (lon, lat, day)
        out = Matrix{Float64}(undef, length(d.is), size(raw, 3))
        @inbounds for t in axes(raw, 3), k in eachindex(d.is)
            out[k, t] = coalesce(raw[d.is[k], d.js[k], t], NaN)
        end
        return out, DateTime.(Date.(ds["valid_time"][:]))
    end
end

function year_forcing(year, d)
    T, dates = read_var("t2m", year, d)
    sf, _ = read_var("sf", year, d); tp, _ = read_var("tp", year, d)
    sw, _ = read_var("ssrd", year, d); lw, _ = read_var("strd", year, d); ws, _ = read_var("ws", year, d)
    td, _ = read_var("d2m", year, d); sp, _ = read_var("sp", year, d)
    to_rate(x) = max.(x, 0.0) .* (1000 / 86_400)        # m w.e. per day -> kg m⁻² s⁻¹
    snow = to_rate(sf); rain = max.(to_rate(tp) .- snow, 0.0)
    fill_nan(x, v) = replace(x, NaN => v)
    T = fill_nan(T, 263.15); td = fill_nan(td, 258.15)
    # Relative humidity with respect to ice, the form the model expects: the vapour pressure from
    # the dewpoint (saturation over water, as in ERA5) over the ice saturation pressure at T.
    e_air = @. 611.21 * exp(17.502 * (td - 273.16) / (td - 32.19))
    rh = clamp.(e_air ./ Chion._bessi_ice_saturation_vapor_pressure.(T, 273.15), 0.0, 1.0)
    return SnowpackForcing(
        dt_days=ones(length(dates)), ncol=length(d.is), time_values=dates,
        air_temperature=T, snowfall_rate=fill_nan(snow, 0.0), rainfall_rate=fill_nan(rain, 0.0),
        shortwave_down=fill_nan(max.(sw, 0.0) ./ 86_400, 0.0), q_lw_down=fill_nan(lw ./ 86_400, 250.0),
        wind_speed=fill_nan(ws, 3.0), relative_humidity=rh, has_relative_humidity=trues(size(rh)),
        air_pressure=fill_nan(sp, 101_325.0), latitude_deg=d.clat, surface_height=d.height,
    )
end

function nh_main()
    d = domain()
    ncol = length(d.is)
    @info "NH domain" columns=ncol area_million_km2=sum(d.area) / 1e6
    grid = SnowpackGrid(ncol; x=d.lon, y=d.lat, js=d.js, is=d.is)
    forcing = Forcing(y -> year_forcing(y, d), SPINUP_YEAR:last(YEARS); grid, area=d.area,
        description="ERA5 daily 1°, $(ERA)")
    model = BESSIModel(forcing; densification=:htessel, fresh_snow_density=:parameterized, ice_substrate_layers=0)
    na = is_north_america.(d.clon)
    regions = (nh=trues(ncol), eurasia=.!na, north_america=na)
    experiment = Experiment(model, forcing; years=YEARS, spinup=CycledSpinup(SPINUP_YEAR:SPINUP_YEAR; cycles=1),
        name="nh_era5_snow", outputs=[
            # snow-covered area (km²) and snow mass (Gt) per region, every day
            Output(joinpath(NOUT, "daily_snow.csv"); variables=[:snow_cover, :swe], frequency=:daily,
                aggregate=:area_sum, regions, snow_cover_threshold=SNOW_THRESHOLD),
            # monthly mean SWE and snow-cover frequency (fraction of days with SWE ≥ threshold) per cell
            Output(joinpath(NOUT, "monthly_snow.nc"); variables=[:swe, :snow_cover], frequency=:monthly,
                snow_cover_threshold=SNOW_THRESHOLD)])
    run!(experiment)
    @info "NH run written" NOUT
end

if abspath(PROGRAM_FILE) == @__FILE__
    nh_main()
end
