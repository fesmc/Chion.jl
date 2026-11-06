#!/usr/bin/env julia

# Evaluation figures for Antarctica (MAR v3.14.3, 27.5 km, ERA5), same style as evaluation_figures.jl:
#   A1  transient AIS-integrated SMB, melt, runoff and refreezing 1979–2025 (left) and Chion − MAR (right), as E1
#   A3  maps (1979–2025 means) per variable: MAR and both runs (top), Chion − MAR below each run
# Reads the runs of run_paper_experiments.sh (output/antarctica_{parameterized,prescribed})
# and the MAR daily files (the forcing of those runs); does not run the model. Writes PDF/PNG figures and
# stats.csv to output/antarctica/figures (or CHION_ANTARCTICA_FIGURE_DIR); captions in README.md there.
#
# Notes on the MAR reference:
#   * MAR Antarctica writes refreezing (RZ) only from 2025 on. The MAR refreezing used here is the
#     surface water budget ME + RF − RU (melt plus rain minus runoff); it ignores changes in liquid-water
#     storage and is checked against RZ for 2025 (stats.csv).
# The per-year MAR reads are cached in output/antarctica/figures/.mar_cache.jls (delete to recompute).
#
#   julia --project=. --threads 8 examples/figures/antarctica_figures.jl
ENV["CHION_DOMAIN"] = "antarctica"
include(joinpath(@__DIR__, "..", "research", "run_gris_fully_prescribed_transient.jl"))
include(joinpath(@__DIR__, "compare_with_mar.jl"))
using Printf, Dates, Serialization

const ROOT = normpath(joinpath(@__DIR__, "..", "..", "output"))
const AOUT = get(ENV, "CHION_ANTARCTICA_FIGURE_DIR", joinpath(ROOT, "antarctica", "figures"))
const PRESCRIBED_DIR = joinpath(ROOT, "antarctica_prescribed")
const PARAMETERIZED_DIR = joinpath(ROOT, "antarctica_parameterized")
const PRESCRIBED_NC = joinpath(PRESCRIBED_DIR, "monthly.nc")
const PARAMETERIZED_NC = joinpath(PARAMETERIZED_DIR, "monthly.nc")
const MAR_DIR = DOMAIN_SETUP.forcing_dir
mar_file(year) = joinpath(MAR_DIR, "MARv3.14.3-27.5km-daily-ERA5-$year.nc")
# Complete MAR years with a file (2026 is partial).
const ALL_YEARS = filter(y -> isfile(mar_file(y)), 1979:2025)
const EVAL_YEARS = ALL_YEARS                      # metrics and maps: 1979–2025, both runs
const LAST_COMPLETE = Dict(:prescribed => 2025, :parameterized => 2025)

# ---------- style (same as evaluation_figures.jl) ----------
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
    mkpath(AOUT)
    savefig(fig, joinpath(AOUT, "$name.pdf")); savefig(fig, joinpath(AOUT, "$name.png"))
    return name
end
const STATS = NamedTuple[]
stat!(figure, quantity, value, unit) = push!(STATS, (; figure, quantity, value, unit))
with_headroom!(p, f=0.14) = (yl = Plots.ylims(p); ylims!(p, (yl[1], yl[2] + f * (yl[2] - yl[1]))); p)

# ---------- small helpers ----------
function read_table(path)
    lines = filter(!isempty, readlines(path))
    header = Symbol.(split(lines[1], ','))
    cells = [split(l, ',') for l in lines[2:end]]
    column(i) = (raw = [c[i] for c in cells]; parsed = tryparse.(Float64, raw);
                 any(isnothing, parsed) ? String.(raw) : Float64.(parsed))
    return Dict(h => column(i) for (i, h) in enumerate(header))
end
const VARS = (:surface_smb => "SMB", :melt => "Melt", :runoff => "Runoff", :refreezing => "Refreezing")
slope_per_decade(x, y) = 10 * sum((x .- mean(x)) .* (y .- mean(y))) / sum((x .- mean(x)) .^ 2)
correlation(a, b) = sum((a .- mean(a)) .* (b .- mean(b))) / sqrt(sum((a .- mean(a)) .^ 2) * sum((b .- mean(b)) .^ 2))
std(x) = sqrt(sum((x .- mean(x)) .^ 2) / (length(x) - 1))
rmse(a, b) = sqrt(mean((a .- b) .^ 2))
quantile_sorted(s, q) = isempty(s) ? 1.0 : s[clamp(round(Int, q * length(s)), 1, length(s))]
nice(x) = (e = 10.0^floor(log10(x)); m = x / e; (m <= 1 ? 1 : m <= 2 ? 2 : m <= 2.5 ? 2.5 : m <= 5 ? 5 : 10) * e)

# ---------- shared MAR / grid data ----------
const GRID_INFO = Ref{Any}(nothing)
function grid_info()
    isnothing(GRID_INFO[]) || return GRID_INFO[]
    p = mar_file(1990); g = load_mar(p).grid
    area = read_mar_static_columns(p, "AREA", g)          # from the 2025 file (DOMAIN_SETUP.area_file)
    GRID_INFO[] = (; grid=g, area)
    return GRID_INFO[]
end
gt(values, area) = sum(values .* area) * 1e-6
to_map(values) = (g = grid_info().grid; M = fill(NaN, length(g.x), length(g.y));
                  for k in eachindex(g.js); M[g.is[k], g.js[k]] = values[k]; end; M)

"""Per year: AIS-integrated MAR ME, RU, RF, SMB, RZ (NaN when absent) in Gt, and per-column annual
totals (mm w.e.) of ME, RU, RF, SMB for EVAL_YEARS. Cached on disk."""
function mar_reference()
    cache = joinpath(AOUT, ".mar_cache.jls")
    isfile(cache) && return deserialize(cache)
    info = grid_info(); g = info.grid; ncol = length(g.js)
    names = (melt="ME", runoff="RU", rain="RF", surface_smb="SMB")
    totals = Dict(v => Float64[] for v in (:melt, :runoff, :rain, :surface_smb, :rz))
    columns = Dict(v => zeros(ncol) for v in keys(names))
    for year in ALL_YEARS
        p = mar_file(year); ntime = NCDataset(ds -> length(ds["TIME"][:]), p)
        for (v, name) in pairs(names)
            annual = vec(sum(read_mar_columns(p, name, g, ntime); dims=2))      # mm w.e. per day, daily files
            push!(totals[v], gt(annual, info.area))
            year in EVAL_YEARS && (columns[v] .+= annual)
        end
        push!(totals[:rz], mar_has(p, "RZ") ? gt(vec(sum(read_mar_columns(p, "RZ", g, ntime); dims=2)), info.area) : NaN)
        @info "MAR year read" year
    end
    for v in keys(columns); columns[v] ./= length(EVAL_YEARS); end
    ref = (; years=collect(ALL_YEARS), totals, columns)
    mkpath(AOUT); serialize(cache, ref)
    return ref
end
"""MAR refreezing from the surface water budget: melt + rain − runoff."""
mar_refreezing(totals) = totals[:melt] .+ totals[:rain] .- totals[:runoff]

"""Chion per-column annual means over EVAL_YEARS from the monthly NetCDF output."""
function chion_columns(path)
    g = grid_info().grid
    out = Dict(v => zeros(length(g.js)) for v in (:melt, :runoff, :refreezing, :surface_smb))
    NCDataset(path) do ds
        t = ds["time"][:]
        sel = findall(d -> year(d) in EVAL_YEARS, t)
        for v in keys(out)
            total = dropdims(sum(coalesce.(ds[String(v)][:, :, sel], 0.0); dims=3); dims=3)   # (x, y), mm w.e.
            out[v] .= [total[g.is[k], g.js[k]] for k in eachindex(g.js)] ./ length(EVAL_YEARS)
        end
    end
    return out
end

# ---------- yearly AIS-integrated series ----------
"""(years, mar, chion) per variable, up to the run's last usable year. The MAR refreezing is
ME + RF − RU from mar_reference (the run CSVs have NaN before 2025)."""
function run_series(dir, mode, ref)
    t = mar_comparison(dir)
    keep = findall(<=(LAST_COMPLETE[mode]), t[:year]); years = t[:year][keep]
    rows = [findfirst(==(y), ref.years) for y in years]
    mar = Dict(v => t[Symbol("mar_$v")][keep] for (v, _) in VARS)
    mar[:refreezing] = mar_refreezing(ref.totals)[rows]
    return (; years, mar, chion=Dict(v => t[Symbol("chion_$v")][keep] for (v, _) in VARS))
end
"""Values on a continuous year axis with NaN in missing years (2021), so lines break there."""
gapped(years, x) = (all = first(years):last(years); v = fill(NaN, length(all));
                    for (y, xi) in zip(years, x); v[Int(y - first(years)) + 1] = xi; end; (all, v))
window(s, years=EVAL_YEARS) = findall(in(years), s.years)
series_bias(s, v, years=EVAL_YEARS) = (k = window(s, years); mean(s.chion[v][k] .- s.mar[v][k]))

# ---------- A1: transient series (left) and residuals (right), as E1 ----------
function figure_transient(runs)
    mar = runs[1][2]                          # the MAR series are the same in both runs (prescribed reaches 2025)
    panels = Any[]
    span = "$(first(EVAL_YEARS))–$(last(EVAL_YEARS))"
    for (i, (v, label)) in enumerate(VARS)
        xlabel = i == length(VARS) ? "Year" : ""
        p = plot(; xlabel, ylabel="$label (Gt yr⁻¹)", legend=false, xlims=(1978, 2026))
        plot!(p, gapped(mar.years, mar.mar[v])...; c=INK, lw=2)
        r = plot(; xlabel, ylabel="Chion − MAR (Gt yr⁻¹)", legend=false, xlims=(1978, 2026))
        hline!(r, [0]; c=MUTED, lw=0.8)
        rm = String[]
        for ((name, s, colour), short) in zip(runs, ("presc.", "param."))
            plot!(p, gapped(s.years, s.chion[v])...; c=colour, lw=1.2)
            plot!(r, gapped(s.years, s.chion[v] .- s.mar[v])...; c=colour, lw=1.2)
            k = window(s); m, x = s.mar[v][k], s.chion[v][k]
            push!(rm, @sprintf("%s %.0f", short, rmse(x, m)))
            stat!("A1", "$label, $name: MAR mean $span", mean(m), "Gt yr⁻¹")
            stat!("A1", "$label, $name: Chion mean $span", mean(x), "Gt yr⁻¹")
            stat!("A1", "$label, $name: mean bias $span", mean(x .- m), "Gt yr⁻¹")
            stat!("A1", "$label, $name: relative mean bias $span", 100 * mean(x .- m) / mean(m), "%")
            stat!("A1", "$label, $name: r $span", correlation(m, x), "1")
            stat!("A1", "$label, $name: RMSE $span", rmse(x, m), "Gt yr⁻¹")
            stat!("A1", "$label, $name: interannual std MAR", std(m), "Gt yr⁻¹")
            stat!("A1", "$label, $name: interannual std Chion", std(x), "Gt yr⁻¹")
            stat!("A1", "$label, $name: trend MAR $span", slope_per_decade(s.years[k], m), "Gt yr⁻¹ per decade")
            stat!("A1", "$label, $name: trend Chion $span", slope_per_decade(s.years[k], x), "Gt yr⁻¹ per decade")
        end
        rmse_text = "RMSE: " * join(rm, ", ") * " Gt yr⁻¹"
        tag!(with_headroom!(p, v == :refreezing ? 0.3 : 0.2), 2i - 1, label)
        v == :refreezing && note!(p, "MAR refreezing: ME + RF − RU"; size=7, row=2.4)
        push!(panels, p)
        push!(panels, note!(tag!(with_headroom!(r, 0.4), 2i, label), rmse_text; size=7, row=2.4))
    end
    entries = vcat([("MAR", INK, (; lw=2))], [(name, colour, (; lw=1.2)) for (name, _, colour) in runs])
    layout = @layout [grid(4, 2); l{0.04h}]
    save_figure(plot(panels..., legend_strip(entries); layout, size=(760, 900), left_margin=4Plots.mm,
        bottom_margin=1Plots.mm), "A1_transient")
end

# ---------- A3: maps per variable (MAR and the two runs) with Chion − MAR below ----------
# Same layout as E3 in evaluation_figures.jl.
const SEQ_CMAP = cgrad(:lajolla; rev=true) # Crameri lajolla: light (low) to dark (high), for non-negative fields
function map_panel(values, clims, cmap, i, label)
    g = grid_info().grid
    p = heatmap(g.x, g.y, permutedims(to_map(values)); c=cmap, clims, aspect_ratio=:equal,
        xlims=extrema(g.x), ylims=(g.y[1] - 0.10 * (g.y[end] - g.y[1]), g.y[end]), framestyle=:none, ticks=false,
        background_color_inside="#ffffff", colorbar=false)
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
    info = grid_info()
    unit = "mm w.e. yr⁻¹"
    robust(v, q) = quantile_sorted(sort(abs.(filter(isfinite, v))), q)
    # Melt and runoff are zero over most of the ice sheet, so the ranges of the fields (98th percentile) and of
    # the differences (90th; SMB 98th, shared by the runs) use only the columns where MAR or a run exceeds
    # 1 mm w.e. yr⁻¹.
    active(v) = reduce((a, c) -> a .| (abs.(c[v]) .> 1), [c for (_, _, c) in cases]; init=abs.(mar[v]) .> 1)
    n = length(cases)
    for (v, label) in VARS
        k = active(v)
        lim = nice(maximum(robust(x[v][k], 0.98) for x in vcat([mar], [c for (_, _, c) in cases])))
        clims, cmap = v == :surface_smb ? ((-lim, lim), SMB_CMAP) : ((0, lim), SEQ_CMAP)
        dlim = nice(maximum(robust((c[v] .- mar[v])[k], v == :surface_smb ? 0.98 : 0.90) for (_, _, c) in cases))
        diff_bar = colour_bar((-dlim, dlim), DIFF_CMAP, "Chion − MAR, $label ($unit)")
        top = Any[map_panel(mar[v], clims, cmap, 1, "MAR")]
        # fields (top) and Chion − MAR below each run (bottom) in one figure
        for (j, (_, name, c)) in enumerate(cases)
            push!(top, map_panel(c[v], clims, cmap, j + 1, name))
        end
        bottom = Any[plot(; framestyle=:none, ticks=false)]
        for (j, (_, name, c)) in enumerate(cases)
            push!(bottom, map_panel(c[v] .- mar[v], (-dlim, dlim), DIFF_CMAP, n + 1 + j,
                @sprintf("Chion − MAR\nRMSE %.0f %s", spatial_rmse(c[v], mar[v]), unit)))
        end
        layout = grid(2, n + 2; widths=[fill(0.955 / (n + 1), n + 1); 0.045])
        save_figure(plot(top..., colour_bar(clims, cmap, "$label ($unit)"), bottom..., diff_bar; layout,
            size=(1000, 760), left_margin=0Plots.mm), "A3_maps_$(v)")
        stat!("A3", "colour range $label", lim, unit); stat!("A3", "colour range Δ$label", dlim, unit)
        for (_, name, c) in cases
            d = c[v] .- mar[v]
            stat!("A3", "$name: area-weighted mean Chion − MAR $label", sum(d .* info.area) / sum(info.area), unit)
            stat!("A3", "$name: area-weighted mean |Chion − MAR| $label", sum(abs.(d) .* info.area) / sum(info.area), unit)
            stat!("A3", "$name: area-weighted spatial RMSE vs MAR, $label", spatial_rmse(c[v], mar[v]), unit)
        end
    end
    for (_, name, c) in cases
        stat!("A3", "$name: spatial r of mean SMB (columns)", correlation(mar[:surface_smb], c[:surface_smb]), "1")
    end
    stat!("A3", "MAR area-weighted mean SMB", sum(mar[:surface_smb] .* info.area) / sum(info.area), unit)
    stat!("A3", "ice-sheet area (selected columns)", sum(info.area), "km²")
    stat!("A3", "number of columns", length(info.area), "1")
end

# ---------- main ----------
function write_stats()
    open(joinpath(AOUT, "stats.csv"), "w") do io
        println(io, "figure,quantity,value,unit")
        for s in STATS
            println(io, join((s.figure, '"' * s.quantity * '"', s.value, '"' * s.unit * '"'), ','))
        end
    end
end

function antarctica_main()
    mkpath(AOUT)
    @info "MAR reference"
    ref = mar_reference()
    # Consistency of the recomputed MAR totals with the run CSVs, and of the refreezing estimate with RZ.
    t = mar_comparison(PRESCRIBED_DIR)
    rows = [findfirst(==(y), t[:year]) for y in ref.years]
    for (v, _) in VARS
        v == :refreezing && continue
        stat!("check", "max |recomputed − CSV| MAR $v $(first(ALL_YEARS))–$(last(ALL_YEARS))",
            maximum(abs.(ref.totals[v] .- t[Symbol("mar_$v")][rows])), "Gt yr⁻¹")
    end
    k = findfirst(==(2025), ref.years)
    stat!("check", "MAR 2025 refreezing RZ", ref.totals[:rz][k], "Gt yr⁻¹")
    stat!("check", "MAR 2025 refreezing ME + RF − RU", mar_refreezing(ref.totals)[k], "Gt yr⁻¹")
    stat!("check", "MAR rainfall mean $(first(EVAL_YEARS))–$(last(EVAL_YEARS))",
        mean(ref.totals[:rain][findall(in(EVAL_YEARS), ref.years)]), "Gt yr⁻¹")
    runs = [("Chion prescribed", run_series(PRESCRIBED_DIR, :prescribed, ref), C1),
            ("Chion parameterized", run_series(PARAMETERIZED_DIR, :parameterized, ref), C3)]
    tp = mar_comparison(PARAMETERIZED_DIR); kp = findfirst(==(2025), tp[:year])
    stat!("check", "parameterized 2025: Chion vapor mass", tp[:chion_vapor_mass][kp], "Gt yr⁻¹")
    stat!("check", "parameterized 2025: Chion − MAR SMB", tp[:chion_surface_smb][kp] - tp[:mar_surface_smb][kp], "Gt yr⁻¹")
    # Vapor mass (deposition − sublimation) bias: with the runoff bias it closes the SMB bias.
    for (name, dir) in (("Chion prescribed", PRESCRIBED_DIR), ("Chion parameterized", PARAMETERIZED_DIR))
        tv = mar_comparison(dir); kv = findall(in(EVAL_YEARS), tv[:year])
        stat!("A1", "Vapor mass, $name: mean bias $(first(EVAL_YEARS))–$(last(EVAL_YEARS))",
            mean(tv[:chion_vapor_mass][kv] .- tv[:mar_vapor_mass][kv]), "Gt yr⁻¹")
    end
    @info "A1"; figure_transient(runs)
    @info "A3"
    mar = Dict(:surface_smb => ref.columns[:surface_smb], :melt => ref.columns[:melt], :runoff => ref.columns[:runoff],
               :refreezing => ref.columns[:melt] .+ ref.columns[:rain] .- ref.columns[:runoff])
    figure_maps(mar, [("prescribed", "Prescribed", chion_columns(PRESCRIBED_NC)),
                      ("parameterized", "Parameterized", chion_columns(PARAMETERIZED_NC))])
    write_stats()
    @info "Antarctica figures written" AOUT
end

if abspath(PROGRAM_FILE) == @__FILE__
    antarctica_main()
end
