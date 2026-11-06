#!/usr/bin/env julia

# ESM-SnowMIP reference sites (Ménard et al., 2019; Krinner et al., 2018) with in situ
# forcing: Chion with the default model, paper-ready figures and statistics.
#   S1  SWE time series at the ten sites
#   S2  mean water-year cycle of SWE and snow depth
#   S3  SWE normalised RMSE against the ESM-SnowMIP models (Krinner et al., 2018, Fig. 4b)
#   S4  peak SWE and melt-out date, simulated against observed, per site and year
# Data: data/ESM-SnowMIP_all (PANGAEA, doi:10.1594/PANGAEA.897575). Model results: output/snowmip/runs
# (examples/runs/snowmip_sites.jl). Output: output/snowmip.
#
#   julia --project=. examples/runs/snowmip_sites.jl
#   julia --project=. examples/figures/snowmip_figures.jl
using Dates, NCDatasets, Plots, Printf, Statistics

const DATA_DIR = joinpath(@__DIR__, "..", "..", "data", "ESM-SnowMIP_all")

const SOUT = get(ENV, "CHION_SNOWMIP_DIR", joinpath(@__DIR__, "..", "..", "output", "snowmip"))
const SITES = [  # code, name, elevation (m), forested
    ("cdp", "Col de Porte", 1325, false), ("oas", "BERMS Old Aspen", 600, true),
    ("obs", "BERMS Old Black Spruce", 629, true), ("ojp", "BERMS Old Jack Pine", 579, true),
    ("rme", "Reynolds Mountain East", 2060, false), ("sap", "Sapporo", 15, false),
    ("snb", "Senator Beck", 3714, false), ("sod", "Sodankylä", 179, false),
    ("swa", "Swamp Angel", 3371, false), ("wfj", "Weissfluhjoch", 2540, false),
]

# ---------- style (as in evaluation_figures.jl) ----------
const C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
const INK, MUTED = "#3a3a38", "#8a8a86"
Plots.default(fontfamily="DejaVu Sans", titlefontsize=9, guidefontsize=8, tickfontsize=7,
    legendfontsize=7, framestyle=:box, grid=false, linewidth=1.5, foreground_color_legend=nothing,
    background_color_legend=nothing, titlelocation=:left, dpi=300)
panel_tag(i) = "($(Char('a' + i - 1)))"
function note!(p, label; loc=:topleft, size=8, row=1)
    (x0, x1), (y0, y1) = Plots.xlims(p), Plots.ylims(p)
    xv = x0 isa Real ? x0 : Dates.value(x0)
    left, top = loc in (:topleft, :bottomleft), loc in (:topleft, :topright)
    xpos = x0 + (left ? 0.03 : 0.97) * (x1 - x0)
    annotate!(p, xpos, y0 + (top ? 1.05 - 0.08 * row : 0.08 * row - 0.05) * (y1 - y0),
        text(label, size, INK, left ? :left : :right, top ? :top : :bottom))
    return p
end
tag!(p, i, label=""; loc=:topleft) = (title!(p, panel_tag(i)); isempty(label) ? p : note!(p, label; loc))
function save_figure(fig, name)
    mkpath(SOUT)
    savefig(fig, joinpath(SOUT, "$name.pdf")); savefig(fig, joinpath(SOUT, "$name.png"))
end
const STATS = NamedTuple[]
stat!(site, quantity, value, unit) = push!(STATS, (; site, quantity, value, unit))

# ---------- per site-year metrics ----------
"""Peak SWE and melt-out day (first day after the peak with SWE < 5 kg m⁻²) per water year,
from daily means; years with observations on less than 60 % of the days from November to
July are skipped (this keeps the sites with automatic, near-daily SWE)."""
function season_metrics(dates, obs, sim)
    od, sd = daily_means(dates, obs), daily_means(dates, sim)
    rows = NamedTuple[]
    for wy in sort(unique(water_year.(Date.(dates))))
        days = Date(wy - 1, 10, 1):Day(1):Date(wy, 9, 30)
        haskey(sd, first(days)) && haskey(sd, last(days)) || continue
        season = Date(wy - 1, 11, 1):Day(1):Date(wy, 7, 31)
        cover = count(d -> haskey(od, d), season) / length(season)
        cover >= 0.6 || continue
        o = [get(od, d, NaN) for d in days]; s = [get(sd, d, NaN) for d in days]
        peak(v) = (finite = filter(isfinite, v); isempty(finite) ? NaN : maximum(finite))
        function melt_out(v)
            ip = argmax(replace(v, NaN => -Inf))
            k = findfirst(i -> i > ip && isfinite(v[i]) && v[i] < 5, eachindex(v))
            isnothing(k) ? NaN : Float64(k)                     # days since 1 October
        end
        push!(rows, (; wy, obs_peak=peak(o), sim_peak=peak(s), obs_off=melt_out(o), sim_off=melt_out(s)))
    end
    return rows
end

# ---------- data ----------
function matching_file(data_dir, prefix, site)
    matches = filter(name -> startswith(name, "$(prefix)_$(site)_") && endswith(name, ".nc"), readdir(data_dir))
    length(matches) == 1 || error("Expected one $(prefix)_$(site)_*.nc file in $(data_dir), found $(length(matches)).")
    joinpath(data_dir, only(matches))
end
read_values(ds, name) = Float64.(coalesce.(ds[name][:], NaN))
first_available(ds, names) = something((n for n in names if haskey(ds, n))..., nothing)

"""Raw meteorology of a site."""
read_met(path) = NCDataset(path) do ds
    (dates=DateTime.(ds["time"][:]), tair=read_values(ds, "Tair"), snow=max.(read_values(ds, "Snowf"), 0.0),
     rain=max.(read_values(ds, "Rainf"), 0.0), sw=max.(read_values(ds, "SWdown"), 0.0), lw=read_values(ds, "LWdown"))
end

function read_obs(path)
    NCDataset(path) do ds
        depth_name = first_available(ds, ["snd_auto", "snd_can_auto", "snd_gap_auto", "snd_gap1_auto", "snd_gap2_auto", "snd_man"])
        swe_name = first_available(ds, ["snw_auto", "snw_man"])
        albedo_name = haskey(ds, "albs") ? "albs" : nothing
        albedo_long_name = isnothing(albedo_name) ? "" : lowercase(String(ds[albedo_name].attrib["long_name"]))
        isnothing(depth_name) && error("No supported snow-depth variable in $(path).")
        dates = DateTime.(ds["time"][:])
        return (dates=dates, depth=read_values(ds, depth_name),
            swe=isnothing(swe_name) ? fill(NaN, length(dates)) : read_values(ds, swe_name),
            albedo=isnothing(albedo_name) ? fill(NaN, length(dates)) : read_values(ds, albedo_name),
            has_surface_albedo=!isnothing(albedo_name) && !occursin("above-canopy", albedo_long_name))
    end
end

"""Observations on the forcing time steps (NaN where there is none)."""
function align_observations(target_dates, obs)
    index = Dict(date => i for (i, date) in pairs(obs.dates))
    pick(v) = [(j = get(index, d, 0); j == 0 ? NaN : v[j]) for d in target_dates]
    return pick(obs.depth), pick(obs.swe), pick(obs.albedo)
end

"""Chion results of a site (output/snowmip/runs/<site>.csv, the state after every forcing step),
shifted so that index i holds the state at the start of forcing step i, the time of the
observations (index 1: the snow-free initial state)."""
function read_simulation(code)
    rows = [parse.(Float64, split(line, ',')[3:5]) for line in Iterators.drop(eachline(joinpath(SOUT, "runs", "$code.csv")), 1)]
    n = length(rows)
    swe = vcat(0.0, [r[1] for r in rows[1:n-1]])
    depth = max.(vcat(0.0, [r[2] for r in rows[1:n-1]]), 0.0)
    albedo = vcat(NaN, [r[3] for r in rows[1:n-1]])
    return (; depth, swe, density=[d > 0 ? w / d : NaN for (w, d) in zip(swe, depth)], albedo)
end

# ---------- statistics ----------
function rmse(observed, simulated)
    valid = isfinite.(observed) .& isfinite.(simulated)
    any(valid) ? sqrt(mean((simulated[valid] .- observed[valid]) .^ 2)) : NaN
end

function daily_means(dates, values)
    sums = Dict{Date,Float64}(); counts = Dict{Date,Int}()
    for (datetime, value) in zip(dates, values)
        isfinite(value) || continue
        date = Date(datetime)
        sums[date] = get(sums, date, 0.0) + value
        counts[date] = get(counts, date, 0) + 1
    end
    Dict(date => total / counts[date] for (date, total) in sums)
end

"""RMSE of the daily means divided by the standard deviation of the observed daily means."""
function normalized_rmse(dates, observed, simulated)
    od, sd = daily_means(dates, observed), daily_means(dates, simulated)
    common = sort!(collect(intersect(keys(od), keys(sd))))
    length(common) > 1 || return NaN
    o = [od[d] for d in common]; s = [sd[d] for d in common]
    σ = std(o; corrected=false)
    σ > 0 ? sqrt(mean((s .- o) .^ 2)) / σ : NaN
end

"""Day of the water year (1 October = 1; 29 February = 0, skipped)."""
function climatology_day(date)
    month(date) == 2 && day(date) == 29 && return 0
    reference = month(date) >= 10 ? Date(2000, month(date), day(date)) : Date(2001, month(date), day(date))
    Dates.value(reference - Date(2000, 10, 1)) + 1
end
water_year(date) = month(date) >= 10 ? year(date) + 1 : year(date)

"""Daily mean curve (365 days) of every water year."""
function water_year_curves(dates, values)
    sums = Dict{Int,Vector{Float64}}(); counts = Dict{Int,Vector{Int}}()
    for (datetime, value) in zip(dates, values)
        isfinite(value) || continue
        date = Date(datetime)
        k = climatology_day(date)
        k == 0 && continue
        wy = water_year(date)
        get!(sums, wy, zeros(365))[k] += value
        get!(counts, wy, zeros(Int, 365))[k] += 1
    end
    Dict(wy => [counts[wy][d] > 0 ? sums[wy][d] / counts[wy][d] : NaN for d in 1:365] for wy in sort(collect(keys(sums))))
end

function curve_mean(curves)
    [(v = [c[d] for c in Base.values(curves) if isfinite(c[d])]; isempty(v) ? NaN : mean(v)) for d in 1:365]
end

function run_all()
    data_dir = abspath(get(ENV, "CHION_SNOWMIP_DATA", DATA_DIR))
    results = Vector{Any}(undef, length(SITES))
    for k in eachindex(SITES)
        code = SITES[k][1]
        met = read_met(matching_file(data_dir, "met_insitu", code))
        obs = read_obs(matching_file(data_dir, "obs_insitu", code))
        depth, swe, albedo = align_observations(met.dates, obs)
        sim = read_simulation(code)
        length(sim.swe) == length(met.dates) || error("The results of $code do not match its forcing; rerun snowmip_sites.jl.")
        results[k] = (; code, met, obs_depth=depth, obs_swe=swe, obs_albedo=albedo,
            has_albedo=obs.has_surface_albedo, sim)
    end
    return results
end

# ---------- figures ----------
daily_series(dates, values) = (d = daily_means(dates, values); k = sort(collect(keys(d))); (k, [d[x] for x in k]))
site_label(k) = "$(uppercase(SITES[k][1])), $(SITES[k][3]) m"

function figure_swe_series(results)
    panels = Any[]
    for (k, r) in enumerate(results)
        od, ov = daily_series(r.met.dates, r.obs_swe); sd, sv = daily_series(r.met.dates, r.sim.swe)
        years = year(first(sd)):year(last(sd)); step = length(years) > 12 ? 4 : 2
        ticks = [Date(y, 1, 1) for y in years if (y - first(years)) % step == 0 && Date(y, 1, 1) >= first(sd)]
        p = plot(; ylabel=isodd(k) ? "SWE (kg m⁻²)" : "", legend=false, xticks=(ticks, string.(year.(ticks))))
        scatter!(p, od, ov; ms=1.2, msw=0, c=MUTED, label="Observed")
        plot!(p, sd, sv; c=C1, lw=1.0, label="Chion")
        ylims!(p, (0, 1.3 * max(maximum(filter(isfinite, ov); init=1), maximum(filter(isfinite, sv); init=1))))
        nr = normalized_rmse(r.met.dates, r.obs_swe, r.sim.swe)
        tag!(p, k, site_label(k)); note!(p, @sprintf("nRMSE %.2f", nr); loc=:topright, size=7)
        push!(panels, p)
    end
    legend = plot(; framestyle=:none, legend=:top, legend_column=2, xlims=(0, 1), ylims=(0, 1))
    scatter!(legend, [NaN], [NaN]; ms=3, msw=0, c=MUTED, label="Observed")
    plot!(legend, [NaN], [NaN]; c=C1, label="Chion")
    layout = @layout [grid(5, 2); l{0.03h}]
    save_figure(plot(panels..., legend; layout, size=(760, 900), left_margin=3Plots.mm), "S1_swe_timeseries")
end

const WY_TICKS = ([1, 93, 183, 274, 365], ["Oct", "Jan", "Apr", "Jul", "Sep"])
function quantile_curve(curves, q)
    [(v = [c[d] for c in values(curves) if isfinite(c[d])]; isempty(v) ? NaN : quantile(v, q)) for d in 1:365]
end
function figure_climatology(results)
    panels = Any[]
    for (k, r) in enumerate(results)
        o = water_year_curves(r.met.dates, r.obs_swe); s = water_year_curves(r.met.dates, r.sim.swe)
        p = plot(; xticks=WY_TICKS, xlims=(1, 365), ylabel=isodd(k) ? "SWE (kg m⁻²)" : "", legend=false)
        # Chion: mean and 10–90 % band over years. Observations: every daily value as a dot,
        # and the mean only on days observed in at least half of the years.
        lo, hi = quantile_curve(s, 0.1), quantile_curve(s, 0.9); ok = isfinite.(lo) .& isfinite.(hi)
        plot!(p, (1:365)[ok], lo[ok]; fillrange=hi[ok], c=C1, fillalpha=0.18, lw=0, label="")
        for curve in values(o)
            scatter!(p, 1:365, curve; ms=0.9, msw=0, c=MUTED, alpha=0.5)
        end
        nyears = length(o)
        coverage = [count(c -> isfinite(c[d]), values(o)) for d in 1:365]
        obs_mean = [coverage[d] >= nyears / 2 ? v : NaN for (d, v) in enumerate(curve_mean(o))]
        plot!(p, 1:365, obs_mean; c=INK, lw=1.6)
        plot!(p, 1:365, curve_mean(s); c=C1, lw=1.6)
        ylims!(p, (0, 1.25 * Plots.ylims(p)[2])); tag!(p, k, site_label(k))
        push!(panels, p)
    end
    legend = plot(; framestyle=:none, legend=:top, legend_column=2, xlims=(0, 1), ylims=(0, 1))
    scatter!(legend, [NaN], [NaN]; ms=3, msw=0, c=MUTED, label="Observed (daily)")
    plot!(legend, [NaN], [NaN]; c=INK, label="Observed mean")
    plot!(legend, [NaN], [NaN]; c=C1, label="Chion mean, 10–90 % band")
    legend = plot!(legend; legend_column=3)
    layout = @layout [grid(5, 2); l{0.03h}]
    save_figure(plot(panels..., legend; layout, size=(760, 900), left_margin=3Plots.mm), "S2_swe_climatology")
end

# Approximate values digitised from the white markers in Fig. 4b of Krinner et al. (2018),
# doi:10.5194/gmd-11-5027-2018 (markers that overlap in the source cannot be separated).
const SNOWMIP_MODELS = Dict(
    "cdp" => [1.65, 1.60, 1.56, 1.41, 1.06, 0.94, 0.91, 0.88, 0.85, 0.80, 0.73, 0.66, 0.60, 0.54, 0.49, 0.44, 0.40],
    "oas" => [1.93, 1.78, 1.74, 1.59, 1.56, 1.52, 1.47, 1.36, 1.25, 1.09, 1.04, 0.96, 0.92, 0.88, 0.83, 0.78],
    "obs" => [2.03, 1.94, 1.90, 1.69, 1.65, 1.60, 1.55, 1.49, 1.39, 1.34, 1.28, 1.21, 1.16, 1.10, 1.04, 0.98, 0.92, 0.87, 0.82, 0.68],
    "ojp" => [1.72, 1.68, 1.47, 1.42, 1.38, 1.29, 1.23, 1.17, 1.11, 1.05, 0.99, 0.94, 0.90, 0.84, 0.68, 0.64, 0.60, 0.57],
    "rme" => [1.23, 1.18, 0.89, 0.73, 0.68, 0.63, 0.57, 0.52, 0.47, 0.42, 0.38, 0.34, 0.30],
    "sap" => [2.04, 2.01, 1.85, 1.75, 1.40, 1.11, 1.01, 0.81, 0.78, 0.75, 0.71, 0.67, 0.64, 0.60, 0.57, 0.54],
    "snb" => [2.06, 2.02, 1.61, 0.92, 0.88, 0.84, 0.80, 0.76, 0.72, 0.68, 0.64, 0.61],
    "sod" => [1.31, 1.10, 0.69, 0.66, 0.63, 0.54, 0.49, 0.44, 0.40, 0.36],
    "swa" => [2.36, 1.99, 1.21, 1.18, 0.85, 0.82, 0.78, 0.75, 0.71, 0.67, 0.63, 0.59, 0.56, 0.52, 0.49, 0.46, 0.42],
    "wfj" => [1.58, 1.54, 1.49, 0.95, 0.91, 0.80, 0.69, 0.65, 0.60, 0.56, 0.53, 0.50, 0.47, 0.44, 0.40, 0.36],
)
function figure_nrmse(results)
    p = plot(; ylabel="SWE normalised RMSE", xticks=(1:length(SITES), uppercase.(first.(SITES))),
        xlims=(0.4, length(SITES) + 0.6), ylims=(0, 2.5), legend=false)
    legend = plot(; framestyle=:none, legend=:top, legend_column=3, xlims=(0, 1), ylims=(0, 1))
    scatter!(legend, [NaN], [NaN]; marker=:circle, ms=3.2, mc=:white, msc=MUTED, msw=0.9, label="ESM-SnowMIP models, open sites")
    scatter!(legend, [NaN], [NaN]; marker=:utriangle, ms=3.2, mc=:white, msc=MUTED, msw=0.9, label="forested sites")
    scatter!(legend, [NaN], [NaN]; marker=:diamond, ms=5.5, c=C1, msw=0, label="Chion")
    hline!(p, [1.0]; c=MUTED, lw=0.8, ls=:dash, label="")
    for (k, (code, _, _, forest)) in enumerate(SITES)
        m = SNOWMIP_MODELS[code]
        scatter!(p, fill(k - 0.12, length(m)), m; marker=forest ? :utriangle : :circle, ms=3.2,
            mc=:white, msc=MUTED, msw=0.9, label="")
        nr = normalized_rmse(results[k].met.dates, results[k].obs_swe, results[k].sim.swe)
        rank = count(<(nr), m) + 1
        scatter!(p, [k + 0.12], [nr]; marker=:diamond, ms=5.5, c=C1, msw=0, label="")
        stat!(code, "SWE nRMSE", nr, "1"); stat!(code, "rank among $(length(m) + 1) (1 = best)", rank, "1")
        stat!(code, "ESM-SnowMIP median nRMSE", median(m), "1")
    end
    layout = @layout [a; l{0.08h}]
    save_figure(plot(p, legend; layout, size=(620, 350), left_margin=3Plots.mm, bottom_margin=3Plots.mm), "S3_swe_nrmse")
end

function figure_season_metrics(results)
    a = plot(; xlabel="Observed peak SWE (kg m⁻²)", ylabel="Chion peak SWE (kg m⁻²)", legend=false, aspect_ratio=:equal)
    b = plot(; xlabel="Observed melt-out date", ylabel="Chion melt-out date", legend=false, aspect_ratio=:equal)
    all_peaks = Float64[]; all_off = Float64[]
    palette = [C1, C2, C3, C4, "#e87ba4"]; used = Tuple{String,String}[]
    for (k, r) in enumerate(results)
        rows = season_metrics(r.met.dates, r.obs_swe, r.sim.swe)
        isempty(rows) && continue
        colour = palette[length(used) + 1]; push!(used, (uppercase(r.code), colour))
        op, sp = getproperty.(rows, :obs_peak), getproperty.(rows, :sim_peak)
        oo, so = getproperty.(rows, :obs_off), getproperty.(rows, :sim_off)
        scatter!(a, op, sp; ms=3.5, msw=0, c=colour)
        scatter!(b, oo, so; ms=3.5, msw=0, c=colour)
        append!(all_peaks, filter(isfinite, vcat(op, sp))); append!(all_off, filter(isfinite, vcat(oo, so)))
        okp = isfinite.(op) .& isfinite.(sp); oko = isfinite.(oo) .& isfinite.(so)
        stat!(r.code, "peak SWE bias (mean over years)", any(okp) ? mean(sp[okp] .- op[okp]) : NaN, "kg m⁻²")
        stat!(r.code, "melt-out date bias (mean over years)", any(oko) ? mean(so[oko] .- oo[oko]) : NaN, "days")
        stat!(r.code, "years evaluated", length(rows), "1")
    end
    lp = (0, 1.05 * maximum(all_peaks)); lo = (extrema(all_off) .+ (-10, 10))
    plot!(a, collect(lp), collect(lp); c=MUTED, ls=:dash, lw=0.8); xlims!(a, lp); ylims!(a, lp)
    plot!(b, collect(lo), collect(lo); c=MUTED, ls=:dash, lw=0.8); xlims!(b, lo); ylims!(b, lo)
    # day-of-water-year ticks as month starts
    starts = [(1, "Oct"), (32, "Nov"), (62, "Dec"), (93, "Jan"), (124, "Feb"), (152, "Mar"), (183, "Apr"), (213, "May"), (244, "Jun"), (274, "Jul"), (305, "Aug")]
    ticks = filter(t -> lo[1] <= t[1] <= lo[2], starts)
    plot!(b; xticks=(first.(ticks), last.(ticks)), yticks=(first.(ticks), last.(ticks)))
    tag!(a, 1, "Peak SWE"); tag!(b, 2, "Melt-out")
    legend = plot(; framestyle=:none, legend=:top, legend_column=length(used), xlims=(0, 1), ylims=(0, 1))
    for (name, colour) in used; scatter!(legend, [NaN], [NaN]; ms=4, msw=0, c=colour, label=name); end
    layout = @layout [grid(1, 2); l{0.07h}]
    save_figure(plot(a, b, legend; layout, size=(720, 380), left_margin=3Plots.mm, bottom_margin=3Plots.mm), "S4_season_metrics")
end

function write_stats(results)
    for (k, r) in enumerate(results)
        stat!(r.code, "SWE RMSE", rmse(r.obs_swe, r.sim.swe), "kg m⁻²")
        stat!(r.code, "snow depth RMSE", rmse(r.obs_depth, r.sim.depth), "m")
        ok = isfinite.(r.obs_swe) .& isfinite.(r.sim.swe)
        stat!(r.code, "SWE bias", any(ok) ? mean(r.sim.swe[ok] .- r.obs_swe[ok]) : NaN, "kg m⁻²")
        okd = isfinite.(r.obs_depth) .& isfinite.(r.sim.depth)
        stat!(r.code, "snow depth bias", any(okd) ? mean(r.sim.depth[okd] .- r.obs_depth[okd]) : NaN, "m")
    end
    open(joinpath(SOUT, "stats.csv"), "w") do io
        println(io, "site,quantity,value,unit")
        for s in STATS; println(io, join((s.site, '"' * s.quantity * '"', s.value, '"' * s.unit * '"'), ',')); end
    end
end

function snowmip_main()
    mkpath(SOUT)
    results = run_all()
    figure_swe_series(results)
    figure_climatology(results)
    figure_nrmse(results)
    figure_season_metrics(results)
    write_stats(results)
    @info "SnowMIP figures written" SOUT
end

if abspath(PROGRAM_FILE) == @__FILE__
    snowmip_main()
end
