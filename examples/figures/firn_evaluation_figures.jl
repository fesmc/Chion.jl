#!/usr/bin/env julia

# Firn evaluation against SUMup (2025 release; Vandecrux et al., doi:10.18739/A2M61BR5M),
# independent of MAR: density profiles and 10 m firn temperatures on the GrIS or, with
# CHION_DOMAIN=antarctica, on the AIS.
# Model: examples/configs/greenland_firn.toml and antarctica_firn.toml, run by
# run_paper_experiments.sh (output/greenland_firn_<mode>, output/antarctica_firn_<mode>), profiles at
# the end of every month (GrIS 1950–2023, AIS 1979–2025) for the parameterized and the prescribed-flux setup.
#   D1  mean density 0–1, 1–5 and 5–10 m, model against SUMup, and the 0–10 m bias against elevation
#   D2  example density profiles at well-known sites
#   D_firn_evaluation  D1, D2 and D3 in one figure
#   D3  10 m firn temperature, model against SUMup, and the bias against elevation
#   D4  (Antarctica) top-10 m firn air content, model against SUMup, and the bias against elevation
# Matching: nearest Chion column within 15 km (Antarctica: 20 km, the 27.5 km grid); nearest monthly
# snapshot in time (density; Antarctic profiles before 1979 are matched to 1979), with the profiles of
# one column and month averaged to one point (site-month); annual mean of the 12 snapshots of the
# year (10 m temperature, compared with the annual mean of the
# observations at 9.5–10.5 m; Antarctic site-years before 1979 use the 1979–2025 model mean).
#
#   julia --project=. --threads 8 examples/figures/firn_evaluation_figures.jl
include(joinpath(@__DIR__, "figure_style.jl"))
include(joinpath(@__DIR__, "run_outputs.jl"))
using NCDatasets

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DOMAIN = Symbol(get(ENV, "CHION_DOMAIN", "greenland"))
const AIS = DOMAIN == :antarctica
const REGION = AIS ? "antarctica" : "greenland"
run_dir(mode) = joinpath(ROOT, "output", "$(REGION)_firn_$mode")
const DOUT = joinpath(ROOT, "output", AIS ? "antarctica_firn" : "gris_firn", "figures")
const SUMUP = joinpath(ROOT, "data", "SUMup")
# parameterized_ntot30: the same parameterized run with 30 instead of the default 15 layers
# (included when present), to show the effect of the layer count on the deep firn.
const ALL_MODES = (("prescribed", "Chion, prescribed fluxes", C1), ("parameterized", "Chion, parameterized", C3),
                   ("parameterized_ntot30", "Chion, parameterized, 30 layers", C2))
const MODES = Tuple(m for m in ALL_MODES if isfile(joinpath(run_dir(m[1]), "firn_profiles.nc")))
const WINDOWS = ((0.0, 1.0), (1.0, 5.0), (5.0, 10.0))
const EXAMPLES = AIS ? ["WDC06A", "DFS2010_20100115_bulk", "BER11C95_25", "DML94C07_38"] :
    ["Summit", "Camp Century", "Dye 2 (Raven) Core A", "NEEM07S3"]
const EXAMPLE_LABELS = Dict("WDC06A" => "WAIS Divide", "DFS2010_20100115_bulk" => "Dome Fuji traverse",
    "BER11C95_25" => "Berkner Island", "DML94C07_38" => "coastal DML",
    "Dye 2 (Raven) Core A" => "Dye-2", "NEEM07S3" => "NEEM")
const OBS_YEARS = AIS ? (1950:2025) : (1950:2023)
const MAX_DISTANCE_KM = AIS ? 20.0 : 15.0
const ELEVATION_BANDS = AIS ? ((0, 1000), (1000, 2500), (2500, 4500)) : ((0, 2000), (2000, 4000))
const T_LIMS = AIS ? (-60, 0) : (-35, 0)

load_model(mode) = read_firn_profiles(run_dir(mode))

haversine(la1, lo1, la2, lo2) = 2 * 6371.0 * asin(sqrt(sind((la2 - la1) / 2)^2 + cosd(la1) * cosd(la2) * sind((lo2 - lo1) / 2)^2))
function nearest_column(m, lat, lon)
    best, d = 0, Inf
    @inbounds for k in eachindex(m.lat)
        abs(m.lat[k] - lat) > 0.3 && continue
        dk = haversine(lat, lon, m.lat[k], m.lon[k]); dk < d && ((best, d) = (k, dk))
    end
    return d <= MAX_DISTANCE_KM ? best : 0
end

function names_of(ds)
    n = ds.group["METADATA"]["name"][:, :]; nkeys = length(ds.group["METADATA"]["name_key"][:])
    rows = size(n, 2) == nkeys ? [n[:, i] for i in axes(n, 2)] : [n[i, :] for i in axes(n, 1)]   # (chars, keys)
    return [strip(String(filter(!=('\0'), collect(r)))) for r in rows]
end

"""SUMup density profiles (one per name, date and position) in OBS_YEARS."""
function sumup_profiles()
    NCDataset(joinpath(SUMUP, "SUMup_2025_density_$(REGION).nc")) do ds
        g = ds.group["DATA"]
        key = Int.(g["name_key"][:]); ts = DateTime.(g["timestamp"][:]); lat = Float64.(g["latitude"][:]); lon = Float64.(g["longitude"][:])
        rho = Float64.(coalesce.(g["density"][:], NaN)); z0 = Float64.(coalesce.(g["start_depth"][:], NaN)); z1 = Float64.(coalesce.(g["stop_depth"][:], NaN))
        mid = Float64.(coalesce.(g["midpoint"][:], NaN))
        names = Dict(zip(Int.(ds.group["METADATA"]["name_key"][:]), names_of(ds)))
        groups = Dict{Tuple{Int,DateTime,Float32,Float32},Vector{Int}}()
        for i in eachindex(key)
            (isfinite(rho[i]) && 100 < rho[i] < 950) || continue
            push!(get!(groups, (key[i], ts[i], Float32(lat[i]), Float32(lon[i])), Int[]), i)
        end
        profiles = []
        for ((k, t, la, lo), idx) in groups
            date = t
            year(date) in OBS_YEARS || continue
            top = [isfinite(z0[i]) ? z0[i] : mid[i] - 0.05 for i in idx]; bot = [isfinite(z1[i]) ? z1[i] : mid[i] + 0.05 for i in idx]
            push!(profiles, (; name=get(names, k, string(k)), date, lat=Float64(la), lon=Float64(lo), top, bot, rho=rho[idx]))
        end
        return profiles
    end
end

"""Mean observed density in [a, b] (overlap-weighted), if at least half of the window is covered."""
function obs_window(p, a, b)
    num = 0.0; cov = 0.0
    for i in eachindex(p.rho)
        ov = max(0.0, min(b, p.bot[i]) - max(a, p.top[i]))
        num += ov * p.rho[i]; cov += ov
    end
    return cov >= 0.5 * (b - a) ? num / cov : NaN
end
model_window(zr, prof, a, b) = (sel = (zr .> a) .& (zr .< b); mean(prof[sel]))
function depth_550(z, rho)
    k = findfirst(>=(550), rho)
    (isnothing(k) || k == 1) && return NaN
    return z[k - 1] + (550 - rho[k - 1]) / (rho[k] - rho[k - 1]) * (z[k] - z[k - 1])
end
function obs_depth_550(p)
    zc = (p.top .+ p.bot) ./ 2; ok = isfinite.(zc)
    count(ok) < 5 && return NaN
    order = sortperm(zc[ok])
    z = zc[ok][order]; r = p.rho[ok][order]
    maximum(r) < 550 && return NaN
    # smooth with a 1 m running median to avoid single ice lenses
    rs = [median(r[abs.(z .- zz) .<= 0.5]) for zz in z]
    return depth_550(z, rs)
end

"""A text-only panel with r and RMSE of each setup as coloured lines, placed under a scatter plot so
that no text covers data. `scores` holds (colour, r, rmse)."""
function score_panel(scores, unit; fmt="%.0f", size=7)
    p = plot(; framestyle=:none, legend=false, xlims=(0, 1), ylims=(0, 1))
    for (row, (colour, r, rmse)) in enumerate(scores)
        label = "r " * @sprintf("%.2f", r) * ", RMSE " * Printf.format(Printf.Format(fmt), rmse) * " " * unit
        annotate!(p, 0.02, 1 - (row - 0.5) / length(scores), text(label, size, colour, :left, :vcenter))
    end
    return p
end

"""
One point per site and month: the observations matched to the same Chion column and monthly
snapshot are averaged (`fo(x)` per matched profile), so repeated sampling at one site, e.g. a pit
every two days through a melt season, counts once. Returns the observed means, the model values
(`fx(x)`, the same for every profile of a group) and the column of each point.
"""
function by_site_month(matched, fo, fx)
    groups = Dict{Tuple{Int,Int},Vector{Float64}}()
    model = Dict{Tuple{Int,Int},Float64}()
    for x in matched
        o = fo(x); isfinite(o) || continue
        key = (x.c, x.s)
        push!(get!(groups, key, Float64[]), o)
        model[key] = fx(x)
    end
    keys_ = sort(collect(keys(groups)))
    return [mean(groups[k]) for k in keys_], [model[k] for k in keys_], first.(keys_)
end

function figure_density(models, profiles)
    panels = [plot(; xlabel="SUMup (kg m⁻³)", ylabel=i == 1 ? "Chion (kg m⁻³)" : "", legend=false, aspect_ratio=:equal) for i in 1:3]
    d550 = plot(; xlabel="Surface elevation (m)", ylabel="Chion − SUMup, 0–10 m (kg m⁻³)", legend=false)
    hline!(d550, [0]; c=MUTED, lw=0.8)
    matched = Dict(mode => [] for (mode, _, _) in MODES)
    scores = [Tuple[] for _ in WINDOWS]
    for (mode, label, colour) in MODES
        m = models[mode]
        for p in profiles
            c = nearest_column(m, p.lat, p.lon); c == 0 && continue
            s = argmin(abs.(Dates.value.(m.time .- p.date)))
            prof = Float64.(m.rho[c, :, s])
            push!(matched[mode], (; p, c, s, prof))
        end
        for (w, (a, b)) in enumerate(WINDOWS)
            o, x, _ = by_site_month(matched[mode], x -> obs_window(x.p, a, b), x -> model_window(m.zr, x.prof, a, b))
            ok = isfinite.(o) .& isfinite.(x)
            scatter!(panels[w], o[ok], x[ok]; c=colour, ms=2.2, msw=0, alpha=0.55)
            push!(scores[w], (colour, cor(x[ok], o[ok]), sqrt(mean((x[ok] .- o[ok]) .^ 2))))
            stat!("D1", "$label: density $(a)–$(b) m bias / RMSE / r (site-months)", NaN,
                @sprintf("%+.0f / %.0f kg m⁻³ / %.2f (%d)", mean(x[ok] .- o[ok]), sqrt(mean((x[ok] .- o[ok]) .^ 2)), cor(x[ok], o[ok]), count(ok)))
        end
        o, x, cols = by_site_month(matched[mode], x -> obs_window(x.p, 0.0, 10.0), x -> model_window(m.zr, x.prof, 0.0, 10.0))
        z = m.elev[cols]
        ok = isfinite.(o) .& isfinite.(x)
        scatter!(d550, z[ok], x[ok] .- o[ok]; c=colour, ms=2.2, msw=0, alpha=0.55)
        stat!("D1", "$label: density 0–10 m bias / RMSE (site-months)", NaN,
            @sprintf("%+.0f / %.0f kg m⁻³ (%d)", mean(x[ok] .- o[ok]), sqrt(mean((x[ok] .- o[ok]) .^ 2)), count(ok)))
        stat!("D1", "$label: SUMup profiles matched to a Chion column", length(matched[mode]), "1")
    end
    lims = ((200, 800), (250, 800), (350, 900))
    for (w, p) in enumerate(panels)
        plot!(p, collect(lims[w]), collect(lims[w]); c=MUTED, ls=:dash, lw=0.8); xlims!(p, lims[w]); ylims!(p, lims[w])
        title!(p, panel_tag(w) * @sprintf("  %g–%g m", WINDOWS[w]...))
    end
    yl = maximum(abs, Plots.ylims(d550)); ylims!(d550, (-yl, 1.35yl))
    tag!(d550, 4, "Bias, 0–10 m"; loc=:topright)
    entries = [(label, colour, (; seriestype=:scatter, marker=:circle, ms=4, msw=0)) for (_, label, colour) in MODES]
    # r and RMSE in a text row under each scatter panel (none under the bias panel)
    notes = [[score_panel(scores[w], "kg m⁻³") for w in eachindex(WINDOWS)]..., plot(; framestyle=:none)]
    layout = @layout [grid(1, 4){0.74h}; grid(1, 4){0.16h}; l{0.1h}]
    save_figure(plot(panels..., d550, notes..., legend_strip(entries); layout, size=(1000, 390), left_margin=3Plots.mm, bottom_margin=1Plots.mm),
        DOUT, "D1_density")
    return matched, [panels..., d550], notes
end

function figure_examples(models, matched)
    panels = Any[]
    for (i, site) in enumerate(EXAMPLES)
        p = plot(; xlabel="Density (kg m⁻³)", ylabel=i == 1 ? "Depth (m)" : "", yflip=true, legend=false, xlims=(150, 950), ylims=(0, 20))
        candidates = filter(x -> x.p.name == site, matched["prescribed"])
        isempty(candidates) && continue
        x = first(sort(candidates; by=x -> -length(x.p.rho)))
        for k in eachindex(x.p.rho)
            (x.p.bot[k] - x.p.top[k] <= 2 && x.p.rho[k] < 900) || continue   # skip coarse bulk intervals and ice values
            plot!(p, [x.p.rho[k], x.p.rho[k]], [x.p.top[k], x.p.bot[k]]; c=MUTED, lw=1.2)
        end
        for (mode, _, colour) in MODES
            m = models[mode]; s = argmin(abs.(Dates.value.(m.time .- x.p.date)))
            plot!(p, m.rho[x.c, :, s], m.zr; c=colour, lw=1.6, linetype=:steppre)
        end
        title!(p, panel_tag(i))
        # site label in the deep, low-density corner (visually bottom left), which no profile reaches
        annotate!(p, 175, 19.4, text("$(get(EXAMPLE_LABELS, site, site)),\n$(year(x.p.date))", 8, INK, :left, :bottom))
        push!(panels, p)
    end
    isempty(panels) && return panels
    entries = vcat([("SUMup", MUTED, (; lw=1.2))], [(label, colour, (; lw=1.6)) for (_, label, colour) in MODES])
    layout = @layout [grid(1, length(panels)); l{0.1h}]
    save_figure(plot(panels..., legend_strip(entries); layout, size=(250 * length(panels) + 60, 380), left_margin=3Plots.mm, bottom_margin=3Plots.mm),
        DOUT, "D2_example_profiles")
    return panels
end

function figure_temperature(models)
    obs = NCDataset(joinpath(SUMUP, "SUMup_2025_temperature_$(REGION).nc")) do ds
        g = ds.group["DATA"]
        depth = Float64.(coalesce.(g["depth"][:], NaN)); T = Float64.(coalesce.(g["temperature"][:], NaN))
        sel = findall(i -> 9.5 <= depth[i] <= 10.5 && isfinite(T[i]), eachindex(depth))
        key = Int.(g["name_key"][sel]); ts = DateTime.(g["timestamp"][sel])
        (; key, year=year.(ts),
           lat=Float64.(g["latitude"][sel]), lon=Float64.(g["longitude"][sel]), T=T[sel])
    end
    groups = Dict{Tuple{Int,Int},Vector{Int}}()
    for i in eachindex(obs.key)
        obs.year[i] in OBS_YEARS && push!(get!(groups, (obs.key[i], obs.year[i]), Int[]), i)
    end
    a = plot(; xlabel="SUMup 10 m temperature (°C)", ylabel="Chion (°C)", legend=false, aspect_ratio=:equal)
    b = plot(; xlabel="Surface elevation (m)", ylabel="Chion − SUMup (K)", legend=false)
    hline!(b, [0]; c=MUTED, lw=0.8)
    scores = Tuple[]
    for (mode, label, colour) in MODES
        m = models[mode]; k10 = findfirst(==(10.0), m.zt)
        o = Float64[]; x = Float64[]; z = Float64[]
        for ((key, y), idx) in groups
            c = nearest_column(m, mean(obs.lat[idx]), mean(obs.lon[idx])); c == 0 && continue
            s = findall(t -> year(t) == y, m.time)
            AIS && isempty(s) && y < year(first(m.time)) && (s = eachindex(m.time))   # pre-1979 sites: model mean
            isempty(s) && continue
            model_T = mean(m.temp[c, k10, s]); isfinite(model_T) || continue
            push!(o, mean(obs.T[idx])); push!(x, model_T - 273.15); push!(z, m.elev[c])
        end
        # observations in K or °C: SUMup temperatures are in °C
        scatter!(a, o, x; c=colour, ms=2.5, msw=0, alpha=0.6)
        scatter!(b, z, x .- o; c=colour, ms=2.5, msw=0, alpha=0.6)
        push!(scores, (colour, cor(x, o), sqrt(mean((x .- o) .^ 2))))
        stat!("D3", "$label: 10 m temperature bias / RMSE / r (site-years)", NaN,
            @sprintf("%+.2f / %.2f K / %.2f (%d)", mean(x .- o), sqrt(mean((x .- o) .^ 2)), cor(x, o), length(o)))
        for (lo, hi) in ELEVATION_BANDS
            sel = lo .<= z .< hi
            any(sel) && stat!("D3", "$label: 10 m temperature bias, $(lo)–$(hi) m", mean((x .- o)[sel]), "K")
        end
    end
    plot!(a, collect(T_LIMS), collect(T_LIMS); c=MUTED, ls=:dash, lw=0.8); xlims!(a, T_LIMS); ylims!(a, T_LIMS)
    tag!(a, 1, "10 m firn temperature"); tag!(b, 2, "Bias"; loc=:topright)
    entries = [(label, colour, (; seriestype=:scatter, marker=:circle, ms=4, msw=0)) for (_, label, colour) in MODES]
    layout = @layout [grid(1, 2){0.74h}; grid(1, 2){0.15h}; l{0.11h}]
    notes = score_panel(scores, "K"; fmt="%.2f")
    save_figure(plot(a, b, notes, plot(; framestyle=:none), legend_strip(entries);
            layout, size=(720, 420), left_margin=3Plots.mm, bottom_margin=1Plots.mm),
        DOUT, "D3_temperature_10m")
    return [a, b], notes
end

"""D1 (top), D2 (middle) and D3 (bottom) in one figure, panels relabelled in reading order."""
function figure_combined(density, density_notes, examples, temperature, temperature_notes)
    # D1 titles carry the depth window after the tag, the others only the tag; notes stay as drawn
    suffix = [[@sprintf("  %g–%g m", w...) for w in WINDOWS]; ""; fill("", length(examples) + 2)]
    panels = [density; examples; temperature]
    for (i, p) in enumerate(panels); title!(p, panel_tag(i) * suffix[i]); end
    plot!(density[end]; left_margin=10Plots.mm)       # room for the bias axis label next to (c)
    blank() = plot(; framestyle=:none)
    examples_row = isempty(examples) ? [blank()] : examples
    layout = grid(6, 1; heights=[0.27, 0.06, 0.31, 0.27, 0.06, 0.03])
    layout[1, 1] = grid(1, 4); layout[2, 1] = grid(1, 4)
    layout[3, 1] = grid(1, length(examples_row))
    layout[4, 1] = grid(1, 2; widths=[0.25, 0.75]); layout[5, 1] = grid(1, 2; widths=[0.25, 0.75])
    entries = vcat([("SUMup profile", MUTED, (; lw=1.2))],
        [(label, colour, (; lw=1.6, marker=:circle, ms=4, msw=0)) for (_, label, colour) in MODES])
    legend_row = legend_strip(entries)
    save_figure(plot(density..., density_notes..., examples_row..., temperature..., temperature_notes, blank(), legend_row;
            layout, size=(1060, 1080), left_margin=3Plots.mm, bottom_margin=1Plots.mm),
        DOUT, "D_firn_evaluation")
end

"""Top-10 m firn air content (m) of an observed profile on a 0.1 m grid, if it covers ≥ 90 % of 0–10 m."""
function obs_fac10(p)
    z = 0.05:0.1:9.95; ρ = fill(NaN, length(z))
    for i in eachindex(p.rho), (k, zz) in enumerate(z)
        p.top[i] <= zz < p.bot[i] && (ρ[k] = p.rho[i])
    end
    count(isfinite, ρ) >= 0.9 * length(z) || return NaN
    return 10.0 * mean(1 .- min.(filter(isfinite, ρ), 917.0) ./ 917.0)
end
model_fac10(zr, prof) = (sel = zr .< 10.0; sum(0.5 .* (1 .- min.(prof[sel], 917.0) ./ 917.0)))

function figure_fac10(models, matched)
    a = plot(; xlabel="SUMup FAC₁₀ (m)", ylabel="Chion FAC₁₀ (m)", legend=false, aspect_ratio=:equal)
    b = plot(; xlabel="Surface elevation (m)", ylabel="Chion − SUMup (m)", legend=false)
    hline!(b, [0]; c=MUTED, lw=0.8)
    for (mode, label, colour) in MODES
        m = models[mode]
        o, x, cols = by_site_month(matched[mode], x -> obs_fac10(x.p), x -> model_fac10(m.zr, x.prof))
        z = m.elev[cols]
        ok = isfinite.(o) .& isfinite.(x)
        scatter!(a, o[ok], x[ok]; c=colour, ms=2.5, msw=0, alpha=0.6)
        scatter!(b, z[ok], x[ok] .- o[ok]; c=colour, ms=2.5, msw=0, alpha=0.6)
        stat!("D4", "$label: FAC10 bias / RMSE / r (site-months)", NaN,
            @sprintf("%+.2f / %.2f m / %.2f (%d)", mean(x[ok] .- o[ok]), sqrt(mean((x[ok] .- o[ok]) .^ 2)), cor(x[ok], o[ok]), count(ok)))
        for (lo, hi) in ELEVATION_BANDS
            sel = ok .& (lo .<= z .< hi)
            any(sel) && stat!("D4", "$label: FAC10 bias / RMSE (n), $(lo)–$(hi) m", NaN,
                @sprintf("%+.2f / %.2f m (%d)", mean(x[sel] .- o[sel]), sqrt(mean((x[sel] .- o[sel]) .^ 2)), count(sel)))
        end
    end
    plot!(a, [0, 7], [0, 7]; c=MUTED, ls=:dash, lw=0.8); xlims!(a, (0, 7)); ylims!(a, (0, 7))
    tag!(a, 1, "Top-10 m firn air content"); tag!(b, 2, "Bias"; loc=:topright)
    entries = [(label, colour, (; seriestype=:scatter, marker=:circle, ms=4, msw=0)) for (_, label, colour) in MODES]
    layout = @layout [grid(1, 2); l{0.1h}]
    save_figure(plot(a, b, legend_strip(entries); layout, size=(720, 360), left_margin=3Plots.mm, bottom_margin=3Plots.mm),
        DOUT, "D4_fac10")
end

function firn_evaluation_main()
    mkpath(DOUT)
    models = Dict(mode => load_model(mode) for (mode, _, _) in MODES)
    length(models) == length(MODES) || return @warn "Firn profiles missing" keys(models)
    profiles = sumup_profiles()
    stat!("D1", "SUMup density profiles $(first(OBS_YEARS))–$(last(OBS_YEARS))", length(profiles), "1")
    matched, density, density_notes = figure_density(models, profiles)
    examples = figure_examples(models, matched)
    temperature, temperature_notes = figure_temperature(models)
    figure_combined(density, density_notes, examples, temperature, temperature_notes)
    AIS && figure_fac10(models, matched)
    write_stats(DOUT)
    @info "Firn evaluation figures written" DOUT
end

if abspath(PROGRAM_FILE) == @__FILE__
    firn_evaluation_main()
end
