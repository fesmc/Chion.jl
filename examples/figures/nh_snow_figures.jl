#!/usr/bin/env julia

# Figures for the ERA5-forced Northern Hemisphere snow run (run_nh_era5_snow.jl).
#   N1  snow-covered area against the Rutgers/NOAA CDR monthly snow-cover extent
#       (Estilow et al., 2015; Robinson et al., 2012): seasonal cycles and NH April/June anomalies,
#       with the six-product blend of Mudryk et al. (2020; data/mudryk2020, ECCC Data Mart)
#   N2  polar maps of mean March SWE and the trend in snow-cover duration from Chion; the variant
#       N2_maps_reference adds GlobSnow v3 (Luojus et al., 2021; data/globsnow_v3) and the gridded
#       weekly NOAA CDR (data/noaa_cdr_sce) as panels c and d
# Rutgers series (data/rutgers_sce, https://climate.rutgers.edu/snowcover/files/moncov.*.txt):
# eurasia, and nam = North America without Greenland. Chion excludes Greenland, so its NH is
# compared with eurasia + nam (nhland and namgnld include Greenland, whose snow-covered area
# varies seasonally between about 1.87 and 2.15 million km²).
# Trends are Sen slopes, their significance a two-sided Mann–Kendall test (5 %).
#
#   julia --project=. examples/figures/nh_snow_figures.jl
using Plots, Printf, Statistics, Dates, NCDatasets
using Plots: RGB

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const NIN = get(ENV, "CHION_NH_OUTPUT_DIR", joinpath(ROOT, "output", "nh_era5_snow"))
const NOUT = joinpath(NIN, "figures")
const RUTGERS = joinpath(ROOT, "data", "rutgers_sce")
const ERA5_NH = get(ENV, "CHION_ERA5_NH_DIR", "/p/projects/ou/labs/ai/Nils/ERA5_nh_daily_1deg")
const BASE = 1981:2010
const YEARS = 1980:2024                   # complete years of both records

const C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
const INK, MUTED = "#3a3a38", "#8a8a86"
Plots.default(fontfamily="DejaVu Sans", titlefontsize=9, guidefontsize=8, tickfontsize=7,
    legendfontsize=7, framestyle=:box, grid=false, linewidth=1.5, foreground_color_legend=nothing,
    background_color_legend=nothing, titlelocation=:left, dpi=300)
panel_tag(i) = "($(Char('a' + i - 1)))"
function note!(p, label; loc=:topleft, size=8, row=1)
    (x0, x1), (y0, y1) = Plots.xlims(p), Plots.ylims(p)
    left, top, centre = loc in (:topleft, :bottomleft), loc in (:topleft, :topright, :topcenter), loc == :topcenter
    annotate!(p, x0 + (centre ? 0.5 : left ? 0.03 : 0.97) * (x1 - x0), y0 + (top ? 1.05 - 0.08 * row : 0.08 * row - 0.05) * (y1 - y0),
        text(label, size, INK, centre ? :center : left ? :left : :right, top ? :top : :bottom))
    return p
end
tag!(p, i, label=""; loc=:topleft) = (title!(p, panel_tag(i)); isempty(label) ? p : note!(p, label; loc))
"""Region name and the annual-mean snow-covered area of each dataset, centred at the top of a
seasonal-cycle panel (above the summer minimum, where the curves leave room)."""
function annual_means!(p, region, items)
    lines = split(region, '\n')
    for (k, l) in enumerate(lines); note!(p, l; loc=:topcenter, size=8, row=k); end
    note!(p, "annual mean:"; loc=:topcenter, size=7, row=length(lines) + 1.2)
    for (k, (name, v)) in enumerate(items)
        note!(p, @sprintf("%s %.1f", name, v); loc=:topcenter, size=7, row=length(lines) + 1.2 + 0.85k)
    end
    return p
end
function save_figure(fig, name)
    mkpath(NOUT)
    savefig(fig, joinpath(NOUT, "$name.pdf")); savefig(fig, joinpath(NOUT, "$name.png"))
end
const STATS = NamedTuple[]
stat!(figure, quantity, value, unit) = push!(STATS, (; figure, quantity, value, unit))
"""Sen slope (median of all pairwise slopes), per decade."""
trend_per_decade(x, y) = 10 * median([(y[j] - y[i]) / (x[j] - x[i]) for i in eachindex(x) for j in (i + 1):lastindex(x)])
"""Mann–Kendall z statistic with tie correction and continuity correction."""
function mann_kendall_z(y)
    n = length(y)
    s = sum(sign(y[j] - y[i]) for i in 1:n for j in (i + 1):n)
    ties = values(Dict(v => count(==(v), y) for v in unique(y)))
    var_s = (n * (n - 1) * (2n + 5) - sum(t * (t - 1) * (2t + 5) for t in ties)) / 18
    return s == 0 || var_s <= 0 ? 0.0 : (s - sign(s)) / sqrt(var_s)
end
const Z_CRIT = 1.959964                                  # two-sided 5 %

function rutgers(name)
    out = Dict{Tuple{Int,Int},Float64}()
    for line in eachline(joinpath(RUTGERS, "moncov.$(name).txt"))
        f = split(line); length(f) == 3 || continue
        out[(parse(Int, f[1]), parse(Int, f[2]))] = parse(Float64, f[3])
    end
    return out
end

"""Mudryk et al. (2020) blended NH snow-cover extent (six products, Greenland excluded), million km²:
(year, month) => (best estimate, lower, upper 95 % bound)."""
mudryk_blend() = NCDataset(joinpath(ROOT, "data", "mudryk2020", "SCE_timeseries.nc")) do ds
    Dict((year(d), month(d)) => (Float64(ds["sce"][i]), Float64(ds["sce_dn"][i]), Float64(ds["sce_up"][i])) for (i, d) in enumerate(ds["time"][:]))
end

# Columns of daily_snow.csv (written by the Experiment output) and the names used below.
const DAILY_COLUMNS = Dict("snow_cover_nh (km2)" => :sca_nh_km2, "snow_cover_eurasia (km2)" => :sca_eurasia_km2,
    "snow_cover_north_america (km2)" => :sca_north_america_km2, "swe_nh (Gt)" => :snow_mass_nh_gt,
    "swe_eurasia (Gt)" => :snow_mass_eurasia_gt, "swe_north_america (Gt)" => :snow_mass_north_america_gt)

"""Monthly means of the daily Chion series (km² and Gt)."""
function chion_monthly()
    lines = readlines(joinpath(NIN, "daily_snow.csv"))
    header = vcat(:date, [DAILY_COLUMNS[String(h)] for h in split(lines[1], ',')[3:end]])
    sums = Dict{Tuple{Int,Int},Vector{Float64}}(); counts = Dict{Tuple{Int,Int},Int}()
    for l in lines[2:end]
        f = split(l, ','); d = Date(DateTime(f[1])); key = (year(d), month(d))
        v = parse.(Float64, f[3:end])
        sums[key] = get(sums, key, zeros(length(v))) .+ v; counts[key] = get(counts, key, 0) + 1
    end
    return Dict(k => Dict(zip(header[2:end], s ./ counts[k])) for (k, s) in sums), header
end

"""Monthly SWE and snow-cover frequency maps, cell areas, and the seasonal-snow mask (cells
snow covered on fewer than half of the August days on average)."""
function monthly_maps()
    lon, lat, swe, freq, time = NCDataset(joinpath(NIN, "monthly_snow.nc")) do ds
        Float64.(ds["x"][:]), Float64.(ds["y"][:]), Float64.(coalesce.(ds["swe"][:, :, :], NaN)),
        Float64.(coalesce.(ds["snow_cover"][:, :, :], NaN)), ds["time"][:]
    end
    lsm = NCDataset(ds -> Float64.(coalesce.(ds["var172"][:, :, 1], 0)), joinpath(ERA5_NH, "static.nc"))
    area = [isfinite(swe[i, j, 1]) ? deg2rad(1)^2 * 6371.0^2 * cosd(lat[j]) * lsm[i, j] : 0.0 for i in eachindex(lon), j in eachindex(lat)]
    yrs, mons = year.(time), month.(time)
    august = dropdims(mean(freq[:, :, mons .== 8]; dims=3); dims=3)
    seasonal = isfinite.(august) .& (august .< 0.5)
    na = [190 <= lon[i] < 310 for i in eachindex(lon), j in eachindex(lat)]
    index(y, m) = findfirst(k -> yrs[k] == y && mons[k] == m, eachindex(time))
    return (; lon, lat, swe, freq, time, yrs, mons, area, seasonal, na, index)
end

function figure_extent()
    chion, _ = chion_monthly()
    eur, nam = rutgers("eurasia"), rutgers("nam")
    blend = mudryk_blend()
    nh = Dict(k => eur[k] + nam[k] for k in keys(eur) if haskey(nam, k))   # NH land without Greenland
    keys_ok(y, m) = haskey(chion, (y, m)) && haskey(eur, (y, m)) && haskey(nam, (y, m))
    panels = Any[]
    # (a) NH seasonal cycle with the Mudryk et al. (2020) blend, which is published for the NH only
    nh_ok(y, m) = haskey(chion, (y, m)) && haskey(nh, (y, m)) && haskey(blend, (y, m))
    nhc(f) = [f([chion[(y, m)][:sca_nh_km2] for y in BASE if nh_ok(y, m)]) / 1e6 for m in 1:12]
    nho(f) = [f([nh[(y, m)] for y in BASE if nh_ok(y, m)]) / 1e6 for m in 1:12]
    nhb(i) = [mean(blend[(y, m)][i] for y in BASE if nh_ok(y, m)) for m in 1:12]
    cm, cs, om, os = nhc(mean), nhc(std), nho(mean), nho(std)
    bm, bdn, bup = nhb(1), nhb(2), nhb(3)
    p = plot(; xlabel="Month", ylabel="Snow-covered area (million km²)", xticks=(1:12, first.(monthname.(1:12), 1)),
        xlims=(0.5, 12.5), legend=false)
    plot!(p, 1:12, om .- os; fillrange=om .+ os, c=MUTED, fillalpha=0.25, lw=0)
    plot!(p, 1:12, bdn; fillrange=bup, c=C2, fillalpha=0.25, lw=0)
    plot!(p, 1:12, bm; c=C2, lw=1.2)
    plot!(p, 1:12, om; c=INK, lw=1.6, ls=:dash)
    plot!(p, 1:12, cm .- cs; fillrange=cm .+ cs, c=C1, fillalpha=0.18, lw=0)
    plot!(p, 1:12, cm; c=C1, lw=1.6, marker=:circle, ms=2.5, msw=0)
    ylims!(p, (0, 1.4 * Plots.ylims(p)[2]))
    title!(p, panel_tag(1))
    annual_means!(p, "Northern Hemisphere\n(without Greenland)", [("Chion", mean(cm)), ("NOAA CDR", mean(om)), ("Mudryk et al.", mean(bm))])
    for m in (1, 4, 5, 6, 10)
        stat!("N1", "NH, $(first(BASE))–$(last(BASE)): $(monthname(m)) Chion / Mudryk blend (95 % range) / Rutgers", NaN,
            @sprintf("%.2f / %.2f (%.2f–%.2f) / %.2f million km²", cm[m], bm[m], bdn[m], bup[m], om[m]))
    end
    push!(panels, p)
    for (i, (label, ckey, obs)) in enumerate((("Eurasia", :sca_eurasia_km2, eur),
                                              ("North America\n(without Greenland)", :sca_north_america_km2, nam)))
        p = plot(; xlabel="Month", xticks=(1:12, first.(monthname.(1:12), 1)), xlims=(0.5, 12.5), legend=false)
        # climatology over the reference period BASE (1981–2010), the anomaly baseline of (c, d)
        cm = [mean(chion[(y, m)][ckey] for y in BASE if keys_ok(y, m)) / 1e6 for m in 1:12]
        om = [mean(obs[(y, m)] for y in BASE if keys_ok(y, m)) / 1e6 for m in 1:12]
        cs = [std([chion[(y, m)][ckey] for y in BASE if keys_ok(y, m)]) / 1e6 for m in 1:12]
        os = [std([obs[(y, m)] for y in BASE if keys_ok(y, m)]) / 1e6 for m in 1:12]
        plot!(p, 1:12, om .- os; fillrange=om .+ os, c=MUTED, fillalpha=0.25, lw=0)
        plot!(p, 1:12, om; c=INK, lw=1.6, ls=:dash)
        plot!(p, 1:12, cm .- cs; fillrange=cm .+ cs, c=C1, fillalpha=0.18, lw=0)
        plot!(p, 1:12, cm; c=C1, lw=1.6, marker=:circle, ms=2.5, msw=0)
        ylims!(p, (0, 1.4 * Plots.ylims(p)[2]))
        title!(p, panel_tag(i + 1))
        annual_means!(p, label, [("Chion", mean(cm)), ("NOAA CDR", mean(om))])
        label = first(split(label, '\n'))
        per = "$(first(BASE))–$(last(BASE))"
        stat!("N1", "$label, $per: mean annual snow-covered area Chion / Rutgers", NaN, @sprintf("%.1f / %.1f million km²", mean(cm), mean(om)))
        for m in (1, 4, 5, 6, 10)
            stat!("N1", "$label, $per: $(monthname(m)) Chion / Rutgers / bias", NaN,
                @sprintf("%.2f / %.2f / %+.2f million km²", cm[m], om[m], cm[m] - om[m]))
        end
        push!(panels, p)
    end

    # NH April and June anomalies relative to 1981–2010, one panel each
    for (k, mon) in enumerate((4, 6))
        name = monthname(mon)
        ys = [y for y in YEARS if haskey(chion, (y, mon)) && haskey(nh, (y, mon))]
        cv = [chion[(y, mon)][:sca_nh_km2] for y in ys] ./ 1e6; ov = [nh[(y, mon)] for y in ys] ./ 1e6
        base = in(BASE).(ys)
        ca = cv .- mean(cv[base]); oa = ov .- mean(ov[base])
        r, tc, to = cor(ca, oa), trend_per_decade(ys, cv), trend_per_decade(ys, ov)
        zc, zo = mann_kendall_z(cv), mann_kendall_z(ov)
        # Mudryk et al. (2020) six-product blend (NH without Greenland, 1967–2018), anomalies to 1981–2010
        ym = [y for y in ys if haskey(blend, (y, mon))]
        bm = [blend[(y, mon)] for y in ym]; bbase = mean(first.(bm)[in(BASE).(ym)])
        ba, bdn, bup = first.(bm) .- bbase, getindex.(bm, 2) .- bbase, getindex.(bm, 3) .- bbase
        p = plot(; xlabel="Year", ylabel=k == 1 ? "Anomaly (million km²)" : "", legend=false, xlims=(first(YEARS), last(YEARS)))
        hline!(p, [0]; c=MUTED, lw=0.6)
        plot!(p, ym, bdn; fillrange=bup, c=C2, fillalpha=0.18, lw=0)
        plot!(p, ym, ba; c=C2, lw=1.2)
        plot!(p, ys, oa; c=INK, lw=1.3, ls=:dash)
        plot!(p, ys, ca; c=C1, lw=1.5)
        lim = 1.1 * maximum(abs, vcat(ca, oa, bdn, bup)); ylims!(p, (-lim, 1.95lim))
        # common period of all three records for the in-panel numbers
        sel = in(ym).(ys); xc = ys[sel]
        tb, tc2, to2 = trend_per_decade(ym, first.(bm)), trend_per_decade(xc, cv[sel]), trend_per_decade(xc, ov[sel])
        tag!(p, 3 + k, "NH, $name"; loc=:topleft)
        note!(p, @sprintf("r of Chion with Mudryk et al. %.2f, with NOAA CDR %.2f", cor(ca[sel], ba), cor(ca[sel], oa[sel])); loc=:topleft, size=7, row=2)
        note!(p, "trend $(first(ym))–$(last(ym)) (million km² per decade):"; loc=:topleft, size=7, row=3)
        note!(p, replace(@sprintf("Chion %+.2f, Mudryk et al. %+.2f, NOAA CDR %+.2f", tc2, tb, to2), "-" => "−"); loc=:topleft, size=7, row=4)
        stat!("N1", "NH $name, $(first(ym))–$(last(ym)): r of Chion with Mudryk et al. / with NOAA CDR (NOAA CDR with Mudryk et al.)", NaN,
            @sprintf("%.2f / %.2f (%.2f)", cor(ca[sel], ba), cor(ca[sel], oa[sel]), cor(oa[sel], ba)))
        stat!("N1", "NH $name, $(first(ym))–$(last(ym)): Sen slope Chion / Mudryk blend / Rutgers", NaN,
            @sprintf("%+.2f / %+.2f / %+.2f million km² per decade (MK z blend %.1f)", tc2, tb, to2, mann_kendall_z(first.(bm))))
        stat!("N1", "NH $name, 1981–2010 mean: Chion / Mudryk blend / Rutgers", NaN,
            @sprintf("%.2f / %.2f / %.2f million km²", mean(cv[base]), bbase, mean(ov[base])))
        stat!("N1", "NH $name anomaly: r Chion vs Rutgers $(first(ys))–$(last(ys))", r, "1")
        stat!("N1", "NH $name, $(first(ys))–$(last(ys)): Sen slope Chion / Rutgers", NaN, @sprintf("%+.2f / %+.2f million km² per decade", tc, to))
        stat!("N1", "NH $name, $(first(ys))–$(last(ys)): Mann–Kendall z Chion / Rutgers", NaN, @sprintf("%.2f / %.2f (|z| > 1.96: p < 0.05)", zc, zo))
        push!(panels, p)
    end

    # March snow mass (statistics only; no observational reference at this resolution), without
    # perennial-snow (glacier) cells: their firn column is still building up after the one-year
    # spin-up and is not seasonal snow.
    m = monthly_maps()
    for (label, sel) in (("NH", trues(size(m.area))), ("Eurasia", .!m.na), ("North America", m.na))
        ys = [y for y in YEARS if !isnothing(m.index(y, 3))]
        v = [sum(filter(isfinite, (m.swe[:, :, m.index(y, 3)] .* m.area)[sel .& m.seasonal])) * 1e-6 for y in ys]
        stat!("N1", "$label March snow mass, seasonal-snow cells, $(first(ys))–$(last(ys)): mean / Sen slope", NaN,
            @sprintf("%.0f Gt / %+.0f Gt per decade", mean(v), trend_per_decade(ys, v)))
    end
    stat!("N1", "perennial-snow cells excluded from snow mass (area)", sum(m.area[.!m.seasonal .& (m.area .> 0)]) / 1e6, "million km²")

    legend = plot(; framestyle=:none, legend=:top, legend_column=3, xlims=(0, 1), ylims=(0, 1))
    plot!(legend, [NaN], [NaN]; c=C1, lw=1.6, marker=:circle, ms=2.5, msw=0, label="Chion")
    plot!(legend, [NaN], [NaN]; c=INK, lw=1.6, ls=:dash, label="NOAA CDR")
    plot!(legend, [NaN], [NaN]; c=C2, lw=1.2, label="Mudryk et al. (2020)")
    layout = @layout [grid(1, 3); grid(1, 2); l{0.05h}]
    save_figure(plot(panels..., legend; layout, size=(820, 600), left_margin=4Plots.mm), "N1_snow_cover")
end

# ---------- polar maps ----------
# Lambert azimuthal equal-area projection centred on the North Pole (unit sphere), Greenwich
# meridian pointing down; this is also the projection of the NSIDC EASE-Grid North (version 1, spherical; GlobSnow).
# Maps are rendered as images (one pixel samples the grid cell it falls in); coastlines come
# from the 0.25° ERA5 land-sea mask, and the graticule is drawn into the image below the data.
const MAP_LAT0 = 25.0
const MAP_PX = 1100
const MAP_R = 2sind((90 - MAP_LAT0) / 2)
const ERA5_LSM = "/p/projects/climate_data_central/reanalysis/ERA5/land_sea_mask_ECMWF-ERA5_observation.nc"
const LAND = RGB(0.86, 0.86, 0.84)
const MASKED = RGB(0.70, 0.70, 0.68)                       # excluded by the reference product
const GRID = RGB(0.74, 0.74, 0.72)
laea(lon, lat) = (ρ = 2sind((90 - lat) / 2); (ρ * sind(lon), -ρ * cosd(lon)))
function laea_inverse(x, y)
    ρ = hypot(x, y)
    return (mod(atand(x, -y), 360), 90 - 2asind(min(ρ / 2, 1)))
end
to_pixel(x, y) = ((x + MAP_R) / 2MAP_R * (MAP_PX - 1) + 1, (MAP_R - y) / 2MAP_R * (MAP_PX - 1) + 1)
const MAP_XS = range(-MAP_R, MAP_R; length=MAP_PX)

"""Longitude/latitude of every map pixel (NaN outside the 25°N circle)."""
function map_pixels()
    lon = fill(NaN, MAP_PX, MAP_PX); lat = fill(NaN, MAP_PX, MAP_PX)
    for c in 1:MAP_PX, r in 1:MAP_PX
        x, y = MAP_XS[c], -MAP_XS[r]
        hypot(x, y) <= MAP_R || continue
        lon[r, c], lat[r, c] = laea_inverse(x, y)
    end
    return lon, lat
end

"""0.25° ERA5 land fraction sampled onto the map pixels."""
function map_land(plon, plat)
    lsm = NCDataset(ds -> Float64.(coalesce.(ds["var172"][:, :, 1], 0)), ERA5_LSM)    # (1440, 721), lat 90 → −90
    return [isfinite(plat[k]) ? lsm[mod(round(Int, plon[k] / 0.25), 1440) + 1, round(Int, (90 - plat[k]) / 0.25) + 1] : 0.0
            for k in CartesianIndices(plon)]
end

"""Dotted graticule (40, 60 and 80°N; meridians every 30° up to 80°N) as a pixel mask."""
function graticule(plon, plat)
    px = (MAP_PX - 1) / 2MAP_R                                # pixels per projected unit
    mask = falses(MAP_PX, MAP_PX)
    for k in CartesianIndices(plat)
        isfinite(plat[k]) || continue
        r = 2sind((90 - plat[k]) / 2) * px                    # distance from the pole in pixels
        for la in (40, 60, 80)
            rl = 2sind((90 - la) / 2) * px
            abs(r - rl) < 0.7 && iseven(floor(Int, deg2rad(plon[k]) * rl / 3)) && (mask[k] = true)
        end
        plat[k] <= 80 || continue
        for lo in 0:30:330
            d = deg2rad(plon[k] - lo)
            abs(r * sin(d)) < 0.7 && cos(d) > 0 && iseven(floor(Int, r / 3)) && (mask[k] = true)
        end
    end
    return mask
end

"""Value of the 1° field (lon 0.5…359.5, lat 25.5…89.5) under every map pixel."""
sample_cells(field, plon, plat) = [isfinite(plat[k]) && plat[k] >= 25 ?
    field[clamp(floor(Int, plon[k]) + 1, 1, 360), clamp(floor(Int, plat[k] - 25) + 1, 1, 65)] : NaN
    for k in CartesianIndices(plon)]

"""Value of an EASE-Grid North field (GlobSnow, 25 km) under every map pixel."""
function sample_ease(field, g, plat)
    dx, dy = g.x[2] - g.x[1], g.y[2] - g.y[1]
    return [isfinite(plat[k]) ? begin
                i = round(Int, (EASE_R * MAP_XS[k[2]] - g.x[1]) / dx) + 1
                j = round(Int, (-EASE_R * MAP_XS[k[1]] - g.y[1]) / dy) + 1
                (1 <= i <= length(g.x) && 1 <= j <= length(g.y)) ? Float64(field[i, j]) : NaN
            end : NaN for k in CartesianIndices(plat)]
end

"""Value of a NOAA CDR field (88 × 88 polar stereographic, true at 60°N, 80°W at the bottom)
under every map pixel; a spherical approximation places the cell centres within 8 km."""
function sample_cdr(field, g, plon, plat)
    return [isfinite(plat[k]) ? begin
                ρ = 6378.137 * (1 + sind(60)) * tand(45 - plat[k] / 2)
                i = round(Int, (ρ * sind(plon[k] + 80) - g.x[1]) / (g.x[2] - g.x[1])) + 1
                j = round(Int, (-ρ * cosd(plon[k] + 80) - g.y[1]) / (g.y[2] - g.y[1])) + 1
                (1 <= i <= length(g.x) && 1 <= j <= length(g.y)) ? field[i, j] : NaN
            end : NaN for k in CartesianIndices(plat)]
end

bin_index(v, edges) = searchsortedlast(edges, v) + 1          # 1 … length(edges) + 1

function polar_map(values, base, edges, colours; label, tag, stipple=nothing, masked=nothing)
    land, grid = base
    img = fill(RGB(1.0, 1.0, 1.0), MAP_PX, MAP_PX)
    for k in CartesianIndices(img)
        land[k] >= 0.5 && (img[k] = LAND)
        grid[k] && (img[k] = GRID)
        !isnothing(masked) && masked[k] && (img[k] = MASKED)
        isfinite(values[k]) && (img[k] = colours[bin_index(values[k], edges)])
    end
    p = plot(img; framestyle=:none, ticks=false, aspect_ratio=:equal, xlims=(0.5, MAP_PX + 0.5),
        ylims=(0.5, MAP_PX + 0.5), legend=false, title=panel_tag(tag))
    contour!(p, 1:MAP_PX, 1:MAP_PX, land; levels=[0.5], c=INK, lw=0.35, colorbar=false)
    for la in (40, 60, 80)                                     # labels over the North Atlantic
        x, y = to_pixel(laea(330, la)...)
        annotate!(p, x, y, text("$(la)°N", 6, MUTED, :center))
    end
    plot!(p, project(0:1:360, fill(MAP_LAT0, 361))...; c=INK, lw=0.6)
    if !isnothing(stipple)
        scatter!(p, project(first.(stipple), last.(stipple))...; ms=0.8, msw=0, c=INK, alpha=0.6)
    end
    annotate!(p, 0.0, 0.0, text(label, 7, INK, :left, :top))
    return p
end
project(lons, lats) = (xy = [to_pixel(laea(lo, la)...) for (lo, la) in zip(lons, lats)]; (first.(xy), last.(xy)))

"""Horizontal colour bar for binned colours; ticks at the inner bin edges."""
function colour_bar(edges, colours, title; low=nothing)
    fmt(e) = e == round(e) ? replace(string(Int(e)), "-" => "−") : string(e)
    ticks = isnothing(low) ? (1:length(edges), fmt.(edges)) : (0:length(edges), fmt.(vcat(low, edges)))
    p = plot(; framestyle=:box, legend=false, xlims=(0, length(colours)), ylims=(0, 1), yticks=false,
        xticks=ticks, xlabel=title, tickfontsize=6, guidefontsize=7, grid=false,
        left_margin=14Plots.mm, right_margin=14Plots.mm, top_margin=-1Plots.mm)
    for (k, c) in enumerate(colours)
        plot!(p, Shape([k - 1, k, k, k - 1], [0, 0, 1, 1]); c, lw=0, linecolor=c)
    end
    return p
end

"""Snow-cover duration trend (Sen slope, days per decade) and its Mann–Kendall significance for
every cell of a series of snow-cover days per hydrological year; cells with fewer than 5 snow days
a year on average or with perennial snow are left out (NaN)."""
function duration_trend(ndays, keep)
    x = collect((first(YEARS) + 1):last(YEARS))
    trend = fill(NaN, size(keep)); sig = falses(size(keep))
    Threads.@threads for k in findall(keep)
        v = ndays[k, :]
        (all(isfinite, v) && mean(v) > 5) || continue
        trend[k] = trend_per_decade(x, v)
        sig[k] = abs(mann_kendall_z(v)) > Z_CRIT
    end
    return trend, sig
end

# ---------- reference data for the maps ----------
const GLOBSNOW = joinpath(ROOT, "data", "globsnow_v3")
const CDR = joinpath(ROOT, "data", "noaa_cdr_sce", "nhsce_v01r01_19661004_20261005.nc")
const EASE_R = 6371228.0                                     # EASE-Grid (version 1) sphere radius (m)

"""GlobSnow v3 March SWE 1981–2010 (mean of the monthly files, pixels with valid SWE in at least
two thirds of the years) and its mountain mask (flag −2; −1 is water, ice sheets or no data)."""
function globsnow_march()
    files = [joinpath(GLOBSNOW, @sprintf("%d03_northern_hemisphere_monthly_swe_0.25grid.nc", y)) for y in BASE]
    x, y = NCDataset(ds -> (Float64.(ds["x"][:]), Float64.(ds["y"][:])), files[1])
    total = zeros(length(x), length(y)); n = zeros(Int, length(x), length(y)); mountain = falses(length(x), length(y))
    for f in files
        s = NCDataset(ds -> Float64.(coalesce.(ds["swe"][:, :], -1)), f)
        ok = s .>= 0
        total[ok] .+= s[ok]; n[ok] .+= 1; mountain .|= s .== -2
    end
    swe = [n[k] >= 2length(files) ÷ 3 ? total[k] / n[k] : NaN for k in CartesianIndices(total)]
    return (; x, y, swe, mountain)
end

"""NOAA CDR weekly snow cover: snow-cover days per hydrological year (fraction of weekly charts
with snow × days in the year), Sen-slope trend and Mann–Kendall significance per land cell."""
function cdr_duration()
    x, y, lat, lon, land, area, t, sce = NCDataset(CDR) do ds
        Float64.(ds["x"][:]), Float64.(ds["y"][:]), Float64.(ds["latitude"][:, :]), mod.(Float64.(ds["longitude"][:, :]), 360),
        ds["land"][:, :] .== 1, Float64.(ds["area"][:, :]), Date.(ds["time"][:]), Float64.(coalesce.(ds["snow_cover_extent"][:, :, :], NaN))
    end
    hyears = (first(YEARS) + 1):last(YEARS)
    ndays = fill(NaN, size(lat)..., length(hyears))
    for (h, hy) in enumerate(hyears)
        weeks = findall(d -> Date(hy - 1, 9, 1) <= d < Date(hy, 9, 1), t)
        ndays[:, :, h] = dropdims(mean(sce[:, :, weeks]; dims=3); dims=3) .* Dates.value(Date(hy, 9, 1) - Date(hy - 1, 9, 1))
    end
    aug = dropdims(mean(sce[:, :, [month(d) == 8 && year(d) in YEARS for d in t]]; dims=3); dims=3)
    greenland(lo, la) = (la >= 59.5 && 300 <= lo <= 335) || (la >= 67 && 300 <= lo <= 350) || (la >= 76 && 288 <= lo <= 300)
    keep = land .& (lat .>= MAP_LAT0) .& (aug .< 0.5) .& .!greenland.(lon, lat)
    trend, sig = duration_trend(ndays, keep)
    # area-mean snow days in autumn (Sep–Dec) and in Jan–Aug of each hydrological year
    function season(h, autumn)
        lo, hi = autumn ? (Date(h - 1, 9, 1), Date(h, 1, 1)) : (Date(h, 1, 1), Date(h, 9, 1))
        frac = dropdims(mean(sce[:, :, findall(d -> lo <= d < hi, t)]; dims=3); dims=3)
        return sum((frac .* Dates.value(hi - lo) .* area)[keep]) / sum(area[keep])
    end
    split = (autumn=[season(h, true) for h in hyears], spring=[season(h, false) for h in hyears])
    return (; x, y, lat, lon, area, trend, sig, split)
end

function figure_maps()
    lon, lat, swe, freq, time = NCDataset(joinpath(NIN, "monthly_snow.nc")) do ds
        Float64.(ds["x"][:]), Float64.(ds["y"][:]), Float64.(coalesce.(ds["swe"][:, :, :], NaN)),
        Float64.(coalesce.(ds["snow_cover"][:, :, :], NaN)), ds["time"][:]
    end
    yrs = year.(time); mons = month.(time)
    mm = monthly_maps(); perennial = .!mm.seasonal; m_area = mm.area
    march_all = dropdims(mean(swe[:, :, (mons .== 3) .& in(BASE).(yrs)]; dims=3); dims=3)
    march_all[perennial] .= NaN
    march = copy(march_all); march[march .< 5] .= NaN                  # snow-free (< 5 kg m⁻²) and glacier cells: land colour
    # Chion: snow-covered days per hydrological year (Sep–Aug), Sen slope and Mann–Kendall test
    hyears = (first(YEARS) + 1):last(YEARS)
    ndays = zeros(length(lon), length(lat), length(hyears))
    for (k, hy) in enumerate(hyears)
        for (y, m) in [(hy - 1, mm) for mm in 9:12] ∪ [(hy, mm) for mm in 1:8]
            t = findfirst(i -> yrs[i] == y && mons[i] == m, eachindex(time))
            isnothing(t) || (ndays[:, :, k] .+= freq[:, :, t] .* Dates.daysinmonth(y, m))
        end
    end
    trend, sig = duration_trend(ndays, .!perennial .& isfinite.(freq[:, :, 1]))
    w = ones(length(lon)) .* cosd.(lat'); ok = isfinite.(trend)
    stat!("N2", "Chion: area-weighted mean trend of snow-covered days (cells with > 5 days)", sum(trend[ok] .* w[ok]) / sum(w[ok]), "days per decade")
    stat!("N2", "Chion: area fraction of trend cells with p < 0.05", sum(w[sig]) / sum(w[ok]), "1")
    stat!("N2", "Chion: area fraction of trend cells with fewer snow days", sum(w[ok][trend[ok] .< 0]) / sum(w[ok]), "1")

    # reference data
    gs = globsnow_march()
    cdr = cdr_duration()
    okc = isfinite.(cdr.trend); wc = cdr.area
    stat!("N2", "NOAA CDR: area-weighted mean trend of snow-covered days", sum(cdr.trend[okc] .* wc[okc]) / sum(wc[okc]), "days per decade")
    stat!("N2", "NOAA CDR: area fraction of trend cells with p < 0.05", sum(wc[cdr.sig]) / sum(wc[okc]), "1")
    stat!("N2", "NOAA CDR: area fraction of trend cells with fewer snow days", sum(wc[okc][cdr.trend[okc] .< 0]) / sum(wc[okc]), "1")
    # autumn (Sep–Dec) and Jan–Aug contributions, area means over the kept cells of each dataset
    keepc = .!perennial .& (m_area .> 0)
    chion_season(h, ms, off) = sum(sum(freq[:, :, findfirst(i -> yrs[i] == h - off && mons[i] == mm, eachindex(time))] .*
        Dates.daysinmonth(h - off, mm) for mm in ms)[keepc] .* m_area[keepc]) / sum(m_area[keepc])
    x = collect(hyears)
    for (name, a, sp) in (("Chion", [chion_season(h, 9:12, 1) for h in hyears], [chion_season(h, 1:8, 0) for h in hyears]),
                          ("NOAA CDR", cdr.split.autumn, cdr.split.spring))
        stat!("N2", "$name: area-mean snow days autumn (Sep–Dec) / Jan–Aug", NaN, @sprintf("%.1f / %.1f days", mean(a), mean(sp)))
        stat!("N2", "$name: Sen slope autumn / Jan–Aug (Mann–Kendall z)", NaN, @sprintf("%+.2f (%.1f) / %+.2f (%.1f) days per decade",
            trend_per_decade(x, a), mann_kendall_z(a), trend_per_decade(x, sp), mann_kendall_z(sp)))
    end
    stat!("N2", "NOAA CDR: autumn snow days, mean of hydrological years ≤ 1999 / ≥ 2000 (IMS since June 1999)", NaN,
        @sprintf("%.1f / %.1f days", mean(cdr.split.autumn[x .<= 1999]), mean(cdr.split.autumn[x .>= 2000])))
    # March SWE on the GlobSnow pixels that are valid in GlobSnow and inside the Chion domain
    pix = 25067.525^2                                        # m², equal-area pixels
    chion_on_gs = fill(NaN, size(gs.swe))
    for j in eachindex(gs.y), i in eachindex(gs.x)
        isfinite(gs.swe[i, j]) || continue
        lo, la = laea_inverse(gs.x[i] / EASE_R, gs.y[j] / EASE_R)
        la >= 25 || continue
        chion_on_gs[i, j] = march_all[clamp(floor(Int, lo) + 1, 1, 360), clamp(floor(Int, la - 25) + 1, 1, 65)]
    end
    common = isfinite.(chion_on_gs) .& isfinite.(gs.swe)
    stat!("N2", "March SWE 1981–2010 on common non-mountain pixels: Chion / GlobSnow", NaN,
        @sprintf("%.0f / %.0f kg m⁻² (mean); %.0f / %.0f Gt (mass); r = %.2f; %.1f million km²",
            mean(chion_on_gs[common]), mean(gs.swe[common]), sum(chion_on_gs[common]) * pix * 1e-12, sum(gs.swe[common]) * pix * 1e-12,
            cor(chion_on_gs[common], gs.swe[common]), count(common) * pix * 1e-12))

    plon, plat = map_pixels(); base = (map_land(plon, plat), graticule(plon, plat))
    swe_edges = [25.0, 50, 100, 150, 200, 300, 500]
    swe_colours = [get(cgrad(:oslo; rev=true), f) for f in range(0.25, 0.97; length=length(swe_edges) + 1)]
    trend_edges = [-10.0, -6, -4, -2, -1, 1, 2, 4, 6, 10]
    trend_colours = [get(cgrad(:BrBG), f) for f in range(0.05, 0.95; length=length(trend_edges) + 1)]
    stipple_b = [(lon[k[1]], lat[k[2]]) for k in findall(sig) if iseven(k[1]) && iseven(k[2])]   # every other cell
    chion_swe, chion_trend = sample_cells(march, plon, plat), sample_cells(trend, plon, plat)
    gs_swe = copy(gs.swe); gs_swe[gs_swe .< 5] .= NaN
    ref_swe, ref_mask = sample_ease(gs_swe, gs, plat), sample_ease(gs.mountain, gs, plat) .== 1
    ref_trend = sample_cdr(cdr.trend, cdr, plon, plat)
    bars() = (colour_bar(swe_edges, swe_colours, "March SWE $(first(BASE))–$(last(BASE)) (kg m⁻²)"; low=5.0),
              colour_bar(trend_edges, trend_colours, "Snow-cover duration trend $(first(hyears))–$(last(hyears)) (days per decade)"))
    chion_panels() = (polar_map(chion_swe, base, swe_edges, swe_colours; label="Chion", tag=1),
                      polar_map(chion_trend, base, trend_edges, trend_colours; label="Chion", tag=2, stipple=stipple_b))
    # N2_maps: Chion only; N2_maps_reference: with GlobSnow v3 and the gridded NOAA CDR (c, d)
    layout = @layout [grid(1, 2){0.935h}; grid(1, 2)]
    save_figure(plot(chion_panels()..., bars()...; layout, size=(760, 430), bottom_margin=4Plots.mm), "N2_maps")
    c = polar_map(ref_swe, base, swe_edges, swe_colours; label="GlobSnow v3", tag=3, masked=ref_mask)
    d = polar_map(ref_trend, base, trend_edges, trend_colours; label="NOAA CDR", tag=4,
        stipple=[(cdr.lon[k], cdr.lat[k]) for k in findall(cdr.sig)])
    layout = @layout [grid(2, 2){0.955h}; grid(1, 2)]
    save_figure(plot(chion_panels()..., c, d, bars()...; layout, size=(760, 820), bottom_margin=4Plots.mm), "N2_maps_reference")
end

function nh_figures_main()
    mkpath(NOUT)
    figure_extent()
    figure_maps()
    open(joinpath(NOUT, "stats.csv"), "w") do io
        println(io, "figure,quantity,value,unit")
        for s in STATS; println(io, join((s.figure, '"' * s.quantity * '"', s.value, '"' * s.unit * '"'), ',')); end
    end
    @info "NH figures written" NOUT
end

if abspath(PROGRAM_FILE) == @__FILE__
    nh_figures_main()
end
