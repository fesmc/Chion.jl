#!/usr/bin/env julia

# Evaluation figures for the model description paper (Greenland, MAR v3.14.3):
#   E1-E3  transient series, residuals and maps: prescribed (daily and native 3-hourly) and parameterized
# Reads the runs of run_paper_experiments.sh (output/greenland_prescribed_steps_8,
# output/greenland_parameterized and the native 3-hourly transient in
# output/gris_substrate5_suite/transient_native_3hourly) and the MAR daily files; writes PDF/PNG figures and stats.csv to
# output/evaluation (or CHION_EVALUATION_DIR). Captions and notes: output/evaluation/README.md.
#
#   julia --project=. --threads 8 examples/figures/evaluation_figures.jl
include(joinpath(@__DIR__, "..", "research", "run_gris_fully_prescribed_transient.jl"))
include(joinpath(@__DIR__, "compare_with_mar.jl"))
using Printf, Dates

const ROOT = normpath(joinpath(@__DIR__, "..", "..", "output"))
const EOUT = get(ENV, "CHION_EVALUATION_DIR", joinpath(ROOT, "evaluation"))
const SUITE = joinpath(ROOT, "gris_substrate5_suite")
PRESCRIBED(steps) = joinpath(ROOT, "greenland_prescribed_steps_$steps")
const PARAMETERIZED = joinpath(ROOT, "greenland_parameterized")
const MAR_DIR = DOMAIN_SETUP.forcing_dir
const EVAL_YEARS = 1980:2017
const LAST_COMPLETE = 2025
mar_file(year) = joinpath(MAR_DIR, "MARv3.14.3-10km-daily-ERA5-$year.nc")

# ---------- style (same as the verification suite) ----------
const C1, C2, C3, C4, C5 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
const INK, MUTED = "#3a3a38", "#8a8a86"
const SMB_CMAP = cgrad(:vik; rev=true)   # Crameri vik: blue accumulation, red-brown ablation
const DIFF_CMAP = cgrad(:bam)            # Crameri bam: green where Chion > MAR, magenta where Chion < MAR
Plots.default(fontfamily="DejaVu Sans", titlefontsize=9, guidefontsize=8, tickfontsize=7,
    legendfontsize=7, framestyle=:box, grid=false, linewidth=1.5, foreground_color_legend=nothing,
    background_color_legend=nothing, titlelocation=:left, dpi=300)
panel_tag(i) = "($(Char('a' + i - 1)))"
"""Text inside the axes at one corner. Call after the limits are final."""
function note!(p, label; loc=:topleft, size=8, row=1)
    (x0, x1), (y0, y1) = Plots.xlims(p), Plots.ylims(p)
    left, top = loc in (:topleft, :bottomleft), loc in (:topleft, :topright)
    annotate!(p, x0 + (left ? 0.03 : 0.97) * (x1 - x0), y0 + (top ? 1.05 - 0.08 * row : 0.08 * row - 0.05) * (y1 - y0),
        text(label, size, INK, left ? :left : :right, top ? :top : :bottom))
    return p
end
"""Panel letter above the axes (top left) and the panel's subject inside the axes."""
tag!(p, i, label=""; loc=:topleft) = (title!(p, panel_tag(i)); isempty(label) ? p : note!(p, label; loc))
"""A legend-only strip placed below multi-panel figures; entries are (label, colour, keywords)."""
function legend_strip(entries)
    p = plot(; framestyle=:none, legend=:top, legend_column=length(entries), xlims=(0, 1), ylims=(0, 1))
    for (label, colour, kw) in entries; plot!(p, [NaN], [NaN]; c=colour, label, kw...); end
    return p
end
function save_figure(fig, name)
    mkpath(EOUT)
    savefig(fig, joinpath(EOUT, "$name.pdf")); savefig(fig, joinpath(EOUT, "$name.png"))
    return name
end
const STATS = NamedTuple[]
stat!(figure, quantity, value, unit) = push!(STATS, (; figure, quantity, value, unit))

# ---------- small helpers ----------
function read_table(path)
    lines = filter(!isempty, readlines(path))
    header = Symbol.(split(lines[1], ','))
    cells = [split(l, ',') for l in lines[2:end]]
    column(i) = (raw = [c[i] for c in cells]; parsed = tryparse.(Float64, raw);
                 any(isnothing, parsed) ? String.(raw) : Float64.(parsed))
    return Dict(h => column(i) for (i, h) in enumerate(header))
end
yearly(dir) = mar_comparison(dir)          # Chion and MAR yearly totals of a run (Gt)
const VARS = (:surface_smb => "SMB", :melt => "Melt", :runoff => "Runoff", :refreezing => "Refreezing")
mar_series(t, v) = t[Symbol("mar_$v")]
chion_series(t, v) = t[Symbol("chion_$v")]
window(t, years) = findall(y -> y in years, t[:year])
mean_bias(t, v, years) = (k = window(t, years); mean(chion_series(t, v)[k] .- mar_series(t, v)[k]))
slope_per_decade(x, y) = 10 * sum((x .- mean(x)) .* (y .- mean(y))) / sum((x .- mean(x)) .^ 2)
correlation(a, b) = sum((a .- mean(a)) .* (b .- mean(b))) / sqrt(sum((a .- mean(a)) .^ 2) * sum((b .- mean(b)) .^ 2))
const BAND_EDGES = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500]
band_label(i) = i == length(BAND_EDGES) - 1 ? "≥$(BAND_EDGES[i])" : "$(BAND_EDGES[i])–$(BAND_EDGES[i + 1])"
band_of(h) = clamp(searchsortedlast(BAND_EDGES, h), 1, length(BAND_EDGES) - 1)

# ---------- shared MAR / grid data ----------
const GRID_INFO = Ref{Any}(nothing)
function grid_info()
    isnothing(GRID_INFO[]) || return GRID_INFO[]
    p = mar_file(1990); l = load_mar(p); g = l.grid
    area = read_mar_static_columns(p, "AREA", g); height = read_mar_static_columns(p, "SH", g)
    GRID_INFO[] = (; grid=g, area, height, bands=band_of.(height))
    return GRID_INFO[]
end
gt(values, area) = sum(values .* area) * 1e-6

"""MAR per-column annual-mean melt, runoff, refreezing and SMB over EVAL_YEARS, plus the
area-integrated monthly climatology of each variable (Gt per calendar month)."""
function mar_column_climatology()
    info = grid_info(); g = info.grid; ncol = length(g.js)
    sums = Dict(v => zeros(ncol) for v in (:melt, :runoff, :refreezing, :surface_smb))
    monthly = Dict(v => zeros(12) for v in keys(sums))
    names = (melt="ME", runoff="RU", refreezing="RZ", surface_smb="SMB")
    for year in EVAL_YEARS
        p = mar_file(year)
        times = NCDataset(ds -> ds["TIME"][:], p); ntime = length(times)
        for (v, name) in pairs(names)
            daily = read_mar_columns(p, name, g, ntime)          # mmWE per day
            sums[v] .+= vec(sum(daily; dims=2))
            for m in 1:12
                monthly[v][m] += gt(vec(sum(daily[:, month.(times) .== m]; dims=2)), info.area)
            end
        end
    end
    n = length(EVAL_YEARS)
    return Dict(v => s ./ n for (v, s) in sums), Dict(v => x ./ n for (v, x) in monthly)
end

"""Chion per-column annual means (EVAL_YEARS) and the monthly climatology of each variable (Gt per
calendar month) from monthly NetCDF output (one file for the whole run, or one file per year)."""
function chion_monthly_netcdf(paths)
    info = grid_info(); g = info.grid
    out = Dict(v => zeros(length(g.js)) for v in (:melt, :runoff, :refreezing, :surface_smb))
    monthly = Dict(v => zeros(12) for v in keys(out))
    for path in (paths isa AbstractString ? [paths] : paths)
        NCDataset(path) do ds
            t = ds[haskey(ds, "time") ? "time" : "t"][:]           # research outputs use `t`
            sel = findall(d -> year(d) in EVAL_YEARS, t)
            isempty(sel) && return
            for v in keys(out)
                field = coalesce.(ds[String(v)][:, :, sel], 0.0)       # (x, y, months), mmWE per month
                total = dropdims(sum(field; dims=3); dims=3)
                out[v] .+= [total[g.is[k], g.js[k]] for k in eachindex(g.js)]
                for (k, s) in enumerate(sel)
                    monthly[v][month(t[s])] += gt([field[g.is[c], g.js[c], k] for c in eachindex(g.js)], info.area)
                end
            end
        end
    end
    n = length(EVAL_YEARS)
    return Dict(v => x ./ n for (v, x) in out), Dict(v => x ./ n for (v, x) in monthly)
end
const NATIVE_MONTHLY = [joinpath(SUITE, "transient_native_3hourly", "monthly", "native_$(y)_monthly.nc") for y in EVAL_YEARS]
PRESCRIBED_NC(steps) = joinpath(PRESCRIBED(steps), "monthly.nc")
const PARAMETERIZED_NC = joinpath(PARAMETERIZED, "monthly.nc")

function to_map(values)
    g = grid_info().grid
    M = fill(NaN, length(g.x), length(g.y))
    for k in eachindex(g.js); M[g.is[k], g.js[k]] = values[k]; end
    return M
end

# ---------- yearly Greenland-integrated series of the three main runs ----------
"""(years, mar, chion) with one vector per variable; years after LAST_COMPLETE are dropped."""
function run_series(dir)
    t = yearly(dir); keep = findall(<=(LAST_COMPLETE), t[:year])
    return (; years=t[:year][keep], mar=Dict(v => mar_series(t, v)[keep] for (v, _) in VARS),
        chion=Dict(v => chion_series(t, v)[keep] for (v, _) in VARS))
end
function native_series()
    t = read_table(joinpath(SUITE, "transient_native_3hourly", "native_3hourly_transient_1980_2017_budgets.csv"))
    return (; years=t[:year], mar=Dict(v => t[Symbol("mar_$(v)_gt")] for (v, _) in VARS),
        chion=Dict(v => t[Symbol("chion_$(v)_gt")] for (v, _) in VARS))
end
series_bias(s, v, years=EVAL_YEARS) = (k = findall(in(years), s.years); mean(s.chion[v][k] .- s.mar[v][k]))
main_runs() = [("Prescribed, daily forcing", run_series(PRESCRIBED(8)), C1),
               ("Prescribed, native 3-hourly", native_series(), C2),
               ("Parameterized", run_series(PARAMETERIZED), C3)]
with_headroom!(p, f=0.14) = (yl = Plots.ylims(p); ylims!(p, (yl[1], yl[2] + f * (yl[2] - yl[1]))); p)

# ---------- E1: transient series (left) and residuals (right) ----------
function figure_transient()
    runs = main_runs()
    mar = runs[1][2]                          # the MAR series are the same in all runs
    panels = Any[]
    for (i, (v, label)) in enumerate(VARS)
        p = plot(; xlabel=i == length(VARS) ? "Year" : "", ylabel="$label (Gt yr⁻¹)", legend=false, xlims=(1940, LAST_COMPLETE))
        plot!(p, mar.years, mar.mar[v]; c=INK, lw=2)
        r = plot(; xlabel=i == length(VARS) ? "Year" : "", ylabel="Chion − MAR (Gt yr⁻¹)", legend=false, xlims=(1940, LAST_COMPLETE))
        hline!(r, [0]; c=MUTED, lw=0.8)
        rmse = String[]
        for ((name, s, colour), short) in zip(runs, ("daily", "3-hourly", "param."))
            plot!(p, s.years, s.chion[v]; c=colour, lw=1.2)
            plot!(r, s.years, s.chion[v] .- s.mar[v]; c=colour, lw=1.2)
            m, x = s.mar[v], s.chion[v]
            push!(rmse, @sprintf("%s %.0f", short, sqrt(mean((x .- m) .^ 2))))
            stat!("E1", "$label, $name: mean bias 1980–2017", series_bias(s, v), "Gt yr⁻¹")
            stat!("E1", "$label, $name: r $(Int(first(s.years)))–$(Int(last(s.years)))", correlation(m, x), "1")
            stat!("E1", "$label, $name: RMSE $(Int(first(s.years)))–$(Int(last(s.years)))", sqrt(mean((x .- m) .^ 2)), "Gt yr⁻¹")
            stat!("E1", "$label, $name: interannual std MAR / Chion", NaN, @sprintf("%.0f / %.0f Gt yr⁻¹", std(m), std(x)))
            stat!("E1", "$label, $name: trend MAR / Chion", NaN,
                @sprintf("%+.1f / %+.1f Gt yr⁻¹ per decade", slope_per_decade(s.years, m), slope_per_decade(s.years, x)))
        end
        push!(panels, tag!(with_headroom!(p, 0.2), 2i - 1, label))
        push!(panels, note!(tag!(with_headroom!(r, 0.4), 2i, label), "RMSE: " * join(rmse, ", ") * " Gt yr⁻¹"; size=7, row=2.4))
    end
    entries = vcat([("MAR", INK, (; lw=2))], [(name, colour, (; lw=1.2)) for (name, _, colour) in runs])
    layout = @layout [grid(4, 2); l{0.04h}]
    save_figure(plot(panels..., legend_strip(entries); layout, size=(760, 900), left_margin=4Plots.mm, bottom_margin=1Plots.mm),
        "E1_transient")
end

# ---------- E2: mean seasonal cycle 1980–2017 ----------
"""Annual totals as a two-column block (name, value) in the top-right or bottom-right corner,
right of the July peak, so the axis range needs no extra room."""
function annual_totals!(p, totals; top=true)
    (x0, x1), (y0, y1) = Plots.xlims(p), Plots.ylims(p)
    n = length(totals)
    for (k, (name, x)) in enumerate(totals)
        fy = top ? 0.95 - 0.08 * (k - 1) : 0.05 + 0.08 * (n - k)
        y = y0 + fy * (y1 - y0); va = top ? :top : :bottom
        annotate!(p, x0 + 0.74 * (x1 - x0), y, text(name, 7, INK, :left, va))
        annotate!(p, x0 + 0.97 * (x1 - x0), y, text(@sprintf("%.0f", x), 7, INK, :right, va))
    end
    return p
end

function figure_seasonal_cycle(mar_monthly, configs)
    short = Dict("Prescribed, daily forcing" => "daily", "Prescribed, native 3-hourly" => "3-hourly", "Parameterized" => "param.")
    panels = Any[]
    for (i, (v, label)) in enumerate(VARS)
        p = plot(; xlabel="Month", ylabel="$label (Gt month⁻¹)", xticks=(1:12, first.(monthname.(1:12), 1)), xlims=(0.5, 12.5), legend=false)
        v == :surface_smb && hline!(p, [0]; c=MUTED, lw=0.8)
        plot!(p, 1:12, mar_monthly[v]; c=INK, lw=2)
        for (name, monthly, colour) in configs
            plot!(p, 1:12, monthly[v]; c=colour, lw=1.4, marker=:circle, ms=2.5, msw=0)
            for m in (6, 7, 8)
                stat!("E2", "$label, $name: $(monthname(m)) Chion − MAR, 1980–2017", monthly[v][m] - mar_monthly[v][m], "Gt month⁻¹")
            end
        end
        # no extra headroom: label and totals go into corners the curves leave empty
        # (SMB: below zero outside summer; the others: above the winter values)
        top = v != :surface_smb
        tag!(p, i, label; loc=top ? :topleft : :bottomleft)
        # annual totals (sum of the monthly climatology), one line per dataset
        totals = vcat([("MAR", sum(mar_monthly[v]))], [(get(short, name, name), sum(m[v])) for (name, m, _) in configs])
        annual_totals!(p, totals; top)
        for (n, x) in totals
            stat!("E2", "$label annual total 1980–2017, $n", x, "Gt yr⁻¹")
        end
        push!(panels, p)
    end
    entries = vcat([("MAR", INK, (; lw=2))], [(name, colour, (; lw=1.4)) for (name, _, colour) in configs])
    layout = @layout [grid(2, 2); l{0.06h}]
    save_figure(plot(panels..., legend_strip(entries); layout, size=(720, 560), left_margin=4Plots.mm, bottom_margin=2Plots.mm),
        "E2_seasonal_cycle")
end

# ---------- E3: maps per variable (MAR and the three runs) with Chion − MAR; SMB differences as SI figure ----------
const MAP_VARS = ((:surface_smb, "SMB"), (:melt, "Melt"), (:runoff, "Runoff"), (:refreezing, "Refreezing"))
const SEQ_CMAP = cgrad(:lajolla; rev=true) # Crameri lajolla: light (low) to dark (high), for non-negative fields
nice(x) = (e = 10.0^floor(log10(x)); m = x / e; (m <= 1 ? 1 : m <= 2 ? 2 : m <= 2.5 ? 2.5 : m <= 5 ? 5 : 10) * e)
robust(v, q=0.98) = quantile_sorted(sort(abs.(filter(isfinite, v))), q)
function map_panel(values, clims, cmap, i, label)
    g = grid_info().grid
    p = heatmap(g.x, g.y, permutedims(to_map(values)); c=cmap, clims, aspect_ratio=:equal, colorbar=false,
        xlims=extrema(g.x), ylims=(g.y[1] - 0.13 * (g.y[end] - g.y[1]), g.y[end]), framestyle=:none,
        ticks=false, background_color_inside="#ffffff")
    return tag!(p, i, label; loc=:bottomleft)
end
"""One vertical colour bar for a row of maps, as its own narrow subplot."""
function colour_bar(clims, cmap, title)
    levels = range(clims..., length=256)
    return heatmap([0, 1], levels, repeat(collect(levels), 1, 2); c=cmap, clims, colorbar=false, xticks=false, ymirror=true,
        ylabel=title, framestyle=:box, tickfontsize=7, guidefontsize=8, ylims=clims, top_margin=8Plots.mm, bottom_margin=8Plots.mm,
        right_margin=6Plots.mm)
end
"""Area-weighted spatial RMSE of `x` against `ref` over the ice-sheet columns."""
spatial_rmse(x, ref) = (a = grid_info().area; sqrt(sum(a .* (x .- ref) .^ 2) / sum(a)))
function figure_maps(mar, cases)
    unit = "mm w.e. yr⁻¹"
    rmse_label(name, x, ref) = name * @sprintf("\nRMSE %.0f mm w.e. yr⁻¹", spatial_rmse(x, ref))
    for (v, label) in MAP_VARS
        # fields: one shared scale for MAR and the three runs
        lim = nice(maximum(robust(x[v]) for x in vcat([mar], [c for (_, _, c) in cases])))
        clims, cmap = v == :surface_smb ? ((-lim, lim), SMB_CMAP) : ((0, lim), SEQ_CMAP)
        dlim = nice(maximum(robust(c[v] .- mar[v], 0.95) for (_, _, c) in cases))
        diff_bar = colour_bar((-dlim, dlim), DIFF_CMAP, "Chion − MAR, $label ($unit)")
        for (_, name, chion) in cases
            stat!("E3", "$name: area-weighted spatial RMSE vs MAR, $label", spatial_rmse(chion[v], mar[v]), unit)
        end
        if v == :surface_smb
            # SMB: fields in the main text, Chion − MAR (one shared scale) as an SI figure
            panels = Any[map_panel(mar[v], clims, cmap, 1, "MAR")]
            for (k, (_, name, chion)) in enumerate(cases)
                push!(panels, map_panel(chion[v], clims, cmap, k + 1, rmse_label(name, chion[v], mar[v])))
            end
            layout = @layout [grid(1, 4) a{0.035w}]
            save_figure(plot(panels..., colour_bar(clims, cmap, "$label ($unit)"); layout, size=(1100, 400),
                left_margin=0Plots.mm), "E3_maps_$(v)")
            diffs = [map_panel(c[v] .- mar[v], (-dlim, dlim), DIFF_CMAP, k, rmse_label(name, c[v], mar[v]))
                     for (k, (_, name, c)) in enumerate(cases)]
            layout = @layout [grid(1, 3) a{0.045w}]
            save_figure(plot(diffs..., diff_bar; layout, size=(860, 400), left_margin=0Plots.mm), "S_E3_differences_$(v)")
        else
            # other variables: fields (top) and Chion − MAR below each run (bottom) in one figure
            n = length(cases)
            top = Any[map_panel(mar[v], clims, cmap, 1, "MAR")]
            for (k, (_, name, chion)) in enumerate(cases)
                push!(top, map_panel(chion[v], clims, cmap, k + 1, name))
            end
            bottom = Any[plot(; framestyle=:none, ticks=false)]
            for (k, (_, name, c)) in enumerate(cases)
                push!(bottom, map_panel(c[v] .- mar[v], (-dlim, dlim), DIFF_CMAP, n + 1 + k,
                    @sprintf("Chion − MAR\nRMSE %.0f mm w.e. yr⁻¹", spatial_rmse(c[v], mar[v]))))
            end
            widths = [fill(0.965 / (n + 1), n + 1); 0.035]
            layout = grid(2, n + 2; widths)
            save_figure(plot(top..., colour_bar(clims, cmap, "$label ($unit)"), bottom..., diff_bar; layout,
                size=(1100, 800), left_margin=0Plots.mm), "E3_maps_$(v)")
        end
        stat!("E3", "colour range $label", lim, unit); stat!("E3", "colour range Δ$label (SI)", dlim, unit)
        for (_, name, c) in cases
            stat!("E3", "$name: area-weighted mean |Chion − MAR| $label", sum(abs.(c[v] .- mar[v]) .* grid_info().area) /
                sum(grid_info().area), unit)
        end
    end
end

quantile_sorted(s, q) = isempty(s) ? 1.0 : s[clamp(round(Int, q * length(s)), 1, length(s))]

std(x) = sqrt(sum((x .- mean(x)) .^ 2) / (length(x) - 1))

# ---------- main ----------
function write_stats()
    open(joinpath(EOUT, "stats.csv"), "w") do io
        println(io, "figure,quantity,value,unit")
        for s in STATS
            println(io, join((s.figure, '"' * s.quantity * '"', s.value, '"' * s.unit * '"'), ','))
        end
    end
end

function evaluation_main()
    mkpath(EOUT)
    @info "E1"; figure_transient()
    @info "map data"
    mar, mar_monthly = mar_column_climatology()
    prescribed, monthly8 = chion_monthly_netcdf(PRESCRIBED_NC(8))
    native, monthly_native = chion_monthly_netcdf(NATIVE_MONTHLY)
    parameterized, monthly_param = chion_monthly_netcdf(PARAMETERIZED_NC)
    @info "E2"; figure_seasonal_cycle(mar_monthly, [("Prescribed, daily forcing", monthly8, C1),
        ("Prescribed, native 3-hourly", monthly_native, C2), ("Parameterized", monthly_param, C3)])
    @info "E3"; figure_maps(mar, [("prescribed_daily", "Prescribed, daily", prescribed),
                                  ("native_3hourly", "Prescribed, native 3-hourly", native),
                                  ("parameterized", "Parameterized", parameterized)])
    write_stats()
    @info "Evaluation figures written" EOUT
end

if abspath(PROGRAM_FILE) == @__FILE__
    evaluation_main()
end
