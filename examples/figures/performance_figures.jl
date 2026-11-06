#!/usr/bin/env julia

# Performance figures for the model description paper (replaces plot_gpu_benchmark.jl, whose
# input files no longer exist):
#   P1_scaling    throughput and wall time per simulated year vs number of columns, GPU vs CPU threads
#   P_performance P1 (top) and P2 (bottom) in one figure, panels (a)–(d)
#   P2_greenland  wall time of one Greenland year (18,620 columns, daily MAR forcing) on GPU and CPU;
#                 panel (b), when present, against CPU threads for Chion and Chion with one daily step
#                 (greenland_raw_substeps1.csv). The original Fortran BESSI
#                 (greenland_raw_bessi_fortran.csv, research/bessi_fortran/run_bessi_fortran_benchmark.sh)
#                 is reported in stats.csv only.
# Single-precision runs (number_type=Float32; *_float32.csv, written by the launchers with
# NUMBER_TYPE=Float32) are added when present: dotted lines in P1, light bars in P2.
# Reads scaling_raw.csv, time_block_raw.csv and greenland_raw.csv written by
# performance_benchmark.jl and performance_greenland.jl (launchers run_performance_benchmark_*.sh)
# and writes PDF/PNG figures and stats.csv to output/performance (or CHION_PERFORMANCE_DIR).
#
#   julia --project=. examples/figures/performance_figures.jl
import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using Plots, Printf
using Statistics: mean

const POUT = get(ENV, "CHION_PERFORMANCE_DIR", normpath(joinpath(@__DIR__, "..", "..", "output", "performance")))
const DEFAULT_NTOT = 15

# ---------- style (same as evaluation_figures.jl) ----------
const C1, C2, C3, C4, C5 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
const INK, MUTED = "#3a3a38", "#8a8a86"
Plots.default(fontfamily="DejaVu Sans", titlefontsize=9, guidefontsize=8, tickfontsize=7,
    legendfontsize=7, framestyle=:box, grid=false, linewidth=1.5, foreground_color_legend=nothing,
    background_color_legend=nothing, titlelocation=:left, dpi=300)
panel_tag(i) = "($(Char('a' + i - 1)))"
"""Text inside the axes at one corner. Call after the limits are final (log axes supported)."""
function note!(p, label; loc=:topleft, size=8, row=1)
    (x0, x1), (y0, y1) = Plots.xlims(p), Plots.ylims(p)
    sp = p[1]
    tx(f) = sp[:xaxis][:scale] == :log10 ? 10^(log10(x0) + f * (log10(x1) - log10(x0))) : x0 + f * (x1 - x0)
    ty(f) = sp[:yaxis][:scale] == :log10 ? 10^(log10(y0) + f * (log10(y1) - log10(y0))) : y0 + f * (y1 - y0)
    left, top = loc in (:topleft, :bottomleft), loc in (:topleft, :topright)
    annotate!(p, tx(left ? 0.03 : 0.97), ty(top ? 1.05 - 0.08 * row : 0.08 * row - 0.05),
        text(label, size, INK, left ? :left : :right, top ? :top : :bottom))
    return p
end
tag!(p, i, label=""; loc=:topleft) = (title!(p, panel_tag(i)); isempty(label) ? p : note!(p, label; loc))
function save_figure(fig, name)
    mkpath(POUT)
    savefig(fig, joinpath(POUT, "$name.pdf")); savefig(fig, joinpath(POUT, "$name.png"))
    return name
end
const STATS = NamedTuple[]
stat!(figure, quantity, value, unit) = push!(STATS, (; figure, quantity, value, unit))

function read_table(path)
    lines = filter(!isempty, readlines(path))
    header = lines[1]
    cells = [split(l, ',') for l in lines[2:end] if l != header]     # concurrent jobs may repeat the header
    names = Symbol.(split(header, ','))
    column(i) = (raw = [c[i] for c in cells]; parsed = tryparse.(Float64, raw);
                 any(isnothing, parsed) ? String.(raw) : Float64.(parsed))
    return Dict(h => column(i) for (i, h) in enumerate(names))
end
rows_where(t, pred) = findall(pred, eachindex(first(values(t))))

"""Per (backend, threads, columns, ntot) statistics over the timed repetitions (probe excluded)."""
function summarize(t)
    keys_ = unique([(t[:backend][i], Int(t[:threads][i]), Int(t[:columns][i]), Int(t[:ntot][i]))
                    for i in eachindex(t[:backend]) if t[:repetition][i] > 0])
    out = NamedTuple[]
    for (backend, threads, columns, ntot) in keys_
        k = rows_where(t, i -> t[:repetition][i] > 0 && t[:backend][i] == backend &&
            t[:threads][i] == threads && t[:columns][i] == columns && t[:ntot][i] == ntot)
        spy = t[:seconds_per_year][k]
        push!(out, (; backend, threads, columns, ntot, n=length(k), years=Int(t[:years][k[1]]),
            mean_spy=mean(spy), min_spy=minimum(spy), mean_cyps=columns / mean(spy), max_cyps=columns / minimum(spy),
            setup=mean(t[:setup_s][k]), melt=mean(t[:mean_melt_kg_m2_yr][k]),
            refreezing=mean(t[:mean_refreezing_kg_m2_yr][k]), runoff=mean(t[:mean_runoff_kg_m2_yr][k]),
            substeps=t[:substep_fraction][k[1]], device=t[:device][k[1]]))
    end
    return sort(out; by=r -> (r.backend, r.threads, r.ntot, r.columns))
end
series(s, backend, threads, ntot) = filter(r -> r.backend == backend && r.threads == threads && r.ntot == ntot, s)
label_of(backend, threads) = backend == "gpu" ? "GPU" : threads == 1 ? "CPU, 1 thread" : "CPU, $threads threads"
const SERIES = [("gpu", 8, C1), ("threads", 1, C2), ("threads", 8, C3), ("threads", 16, C4), ("threads", 32, C5)]

function fig_scaling(s, s32=NamedTuple[])
    xt = [10.0^k for k in 2:6]
    xtl = ["10²", "10³", "10⁴", "10⁵", "10⁶"]
    pa = plot(; xscale=:log10, yscale=:log10, xlabel="Number of columns", ylabel="Column-years per second",
        legend=false, xticks=(xt, xtl), xlims=(60, 5e6))
    pb = plot(; xscale=:log10, yscale=:log10, xlabel="Number of columns", ylabel="Wall time per simulated year (s)",
        legend=:bottomright, xticks=(xt, xtl), xlims=(60, 5e6))
    for (backend, threads, colour) in SERIES
        r = series(s, backend, threads, DEFAULT_NTOT)
        isempty(r) && continue
        x = [q.columns for q in r]
        plot!(pa, x, [q.mean_cyps for q in r]; c=colour, marker=:circle, ms=3, msw=0, label="")
        plot!(pb, x, [q.mean_spy for q in r]; c=colour, marker=:circle, ms=3, msw=0, label=label_of(backend, threads))
        r30 = series(s, backend, threads, 30)
        isempty(r30) && continue
        plot!(pa, [q.columns for q in r30], [q.mean_cyps for q in r30]; c=colour, ls=:dash, marker=:utriangle,
            ms=3, msw=0, label="")
        plot!(pb, [q.columns for q in r30], [q.mean_spy for q in r30]; c=colour, ls=:dash, marker=:utriangle,
            ms=3, msw=0, label="")
    end
    for (backend, threads, colour) in SERIES
        r = series(s32, backend, threads, DEFAULT_NTOT)
        isempty(r) && continue
        x = [q.columns for q in r]
        plot!(pa, x, [q.mean_cyps for q in r]; c=colour, ls=:dot, marker=:diamond, ms=3, msw=0, label="")
        plot!(pb, x, [q.mean_spy for q in r]; c=colour, ls=:dot, marker=:diamond, ms=3, msw=0, label="")
    end
    plot!(pb, [NaN], [NaN]; c=MUTED, ls=:dash, marker=:utriangle, ms=3, msw=0, label="Ntot = 30")
    isempty(s32) || plot!(pb, [NaN], [NaN]; c=MUTED, ls=:dot, marker=:diamond, ms=3, msw=0, label="Float32")
    all_ = [s; s32]
    ylims!(pa, (10.0^floor(log10(minimum(q.mean_cyps for q in all_)) - 0.2), 10.0^(ceil(log10(maximum(q.mean_cyps for q in all_))) + 0.4)))
    ylims!(pb, (10.0^(floor(log10(minimum(q.mean_spy for q in all_))) - 0.6), 10.0^(ceil(log10(maximum(q.mean_spy for q in all_))) + 0.3)))
    tag!(pa, 1, "Throughput, daily steps, Ntot = $DEFAULT_NTOT")
    tag!(pb, 2, "Time per simulated year")
    fig = plot(pa, pb; layout=(1, 2), size=(900, 360), left_margin=4Plots.mm, bottom_margin=4Plots.mm)
    save_figure(fig, "P1_scaling")
    return [pa, pb]
end

thousands(x::Integer) = replace(string(x), r"(?<=\d)(?=(\d{3})+$)" => ",")

function greenland_summary(g; precision="Float64", phase="one_year_from_warm_state")
    keys_ = unique([(g[:backend][i], Int(g[:threads][i])) for i in eachindex(g[:backend])])
    out = NamedTuple[]
    for (backend, threads) in keys_
        k = rows_where(g, i -> g[:backend][i] == backend && g[:threads][i] == threads &&
            g[:phase][i] == phase)
        w = rows_where(g, i -> g[:backend][i] == backend && g[:threads][i] == threads &&
            g[:phase][i] == "warmup_from_initial_state")
        total = g[:setup_s][k] .+ g[:integration_s][k] .+ g[:finalize_s][k]
        push!(out, (; backend, threads, n=length(k), mean_integration=mean(g[:integration_s][k]),
            min_integration=minimum(g[:integration_s][k]), mean_total=mean(total), min_total=minimum(total),
            setup=mean(g[:setup_s][k]), io=mean(g[:forcing_io_s][[k; w]]), warmup=g[:integration_s][w[1]],
            columns=Int(g[:columns][k[1]]), substeps=g[:substep_fraction][k[1]], device=g[:device][k[1]], precision))
    end
    return sort(out; by=r -> (r.backend != "gpu", r.threads, r.precision != "Float64"))
end

function fig_greenland(all_rows, comparison=())
    # with panel (b), the single-thread time (15 s) is shown there only, so that (a) stays readable
    gs = isempty(comparison) ? all_rows : filter(r -> r.threads != 1 || r.backend == "gpu", all_rows)
    colour(r) = r.backend == "gpu" ? C1 : r.threads == 16 ? C4 : r.threads == 32 ? C5 : r.threads == 8 ? C3 : C2
    labels = [label_of(r.backend, r.threads) * (r.precision == "Float64" ? "" : ", $(r.precision)") for r in gs]
    vals = [r.mean_total for r in gs]
    cols = [colour(r) for r in gs]
    single = [r.precision != "Float64" for r in gs]
    n = length(vals)
    ypos = collect(n:-1:1)             # first entry on top
    p = plot(; legend=false, yticks=(ypos, labels), ytick_direction=:none, xlabel="Wall time (s)", xlims=(0, 1.5 * maximum(vals)),
        ylims=(0.4, n + 1.0))
    for (sel, alpha) in ((.!single, 1.0), (single, 0.45))   # single precision as lighter bars
        any(sel) && bar!(p, ypos[sel], vals[sel]; orientation=:horizontal, c=cols[sel], lc=cols[sel],
            fillalpha=alpha, bar_width=0.6)
    end
    for (i, v) in enumerate(vals)
        annotate!(p, v + 0.015 * maximum(vals), ypos[i],
            text(@sprintf("%.2f s  (%s col-yr/s)", v, thousands(round(Int, round(gs[i].columns / v; sigdigits=3)))), 7, INK, :left, :vcenter))
    end
    note!(p, "Greenland, one simulated year, $(thousands(gs[1].columns)) columns"; loc=:topright)
    if isempty(comparison)
        save_figure(plot(p; size=(620, 160 + 30 * n), left_margin=3Plots.mm, bottom_margin=3Plots.mm), "P2_greenland")
        return [p]
    end
    title!(p, panel_tag(1))
    pt = fig_greenland_threads(all_rows, comparison...)
    save_figure(plot(p, pt; layout=grid(1, 2; widths=[0.58, 0.42]), size=(1050, max(360, 160 + 30 * n)),
        left_margin=3Plots.mm, bottom_margin=5Plots.mm), "P2_greenland")
    return [p, pt]
end

"""P1 (top row) and P2 (bottom row) as one figure, panels relabelled (a)–(d)."""
function fig_performance(scaling, greenland)
    panels = [scaling; greenland]
    for (i, p) in enumerate(panels); title!(p, panel_tag(i)); end
    layout = grid(2, 1)                  # rows nested, so the bottom widths do not apply to the top row
    layout[1, 1] = grid(1, 2)
    layout[2, 1] = length(greenland) == 2 ? grid(1, 2; widths=[0.58, 0.42]) : grid(1, 1)
    return save_figure(plot(panels...; layout, size=(1050, 800), left_margin=4Plots.mm, bottom_margin=4Plots.mm),
        "P_performance")
end

"""Wall time against CPU threads: Chion (default and one daily step), optionally the original Fortran BESSI."""
function fig_greenland_threads(gs, one_step, fortran)
    cpu(rows) = sort(filter(r -> r.backend != "gpu" && r.precision == "Float64", rows); by=r -> r.threads)
    p = plot(; xscale=:log10, yscale=:log10, xlabel="CPU threads", ylabel="Wall time per simulated year (s)",
        xticks=([1, 2, 4, 8, 16, 32], ["1", "2", "4", "8", "16", "32"]), xlims=(0.8, 45), legend=:topright)
    for (rows, colour, ls, marker, label) in ((cpu(gs), INK, :solid, :circle, "Chion"),
            (cpu(one_step), INK, :dash, :utriangle, "Chion, one daily step"),
            (fortran, MUTED, :solid, :square, "BESSI (Fortran)"))
        isempty(rows) && continue
        plot!(p, [r.threads for r in rows], [r.mean_total for r in rows]; c=colour, ls, marker, ms=4, msw=0, label)
    end
    gpu = filter(r -> r.backend == "gpu" && r.precision == "Float64", gs)
    isempty(gpu) || hline!(p, [only(gpu).mean_total]; c=C1, ls=:dot, label="Chion, GPU")
    title!(p, panel_tag(2))
    return p
end

function main()
    s = summarize(read_table(joinpath(POUT, "scaling_raw.csv")))
    open(joinpath(POUT, "scaling_summary.csv"), "w") do io
        println(io, "backend,threads,columns,ntot,repetitions,years_per_repetition,mean_seconds_per_year," *
            "min_seconds_per_year,mean_column_years_per_s,max_column_years_per_s,mean_setup_s,mean_melt_kg_m2_yr," *
            "mean_refreezing_kg_m2_yr,mean_runoff_kg_m2_yr,substep_fraction")
        for r in s
            println(io, join((r.backend, r.threads, r.columns, r.ntot, r.n, r.years, r.mean_spy, r.min_spy, r.mean_cyps,
                r.max_cyps, r.setup, r.melt, r.refreezing, r.runoff, r.substeps), ','))
        end
    end
    path32 = joinpath(POUT, "scaling_raw_float32.csv")
    s32 = isfile(path32) ? summarize(read_table(path32)) : NamedTuple[]
    scaling_panels = fig_scaling(s, s32)
    for (backend, threads, _) in SERIES
        r32 = series(s32, backend, threads, DEFAULT_NTOT); r64 = series(s, backend, threads, DEFAULT_NTOT)
        (isempty(r32) || isempty(r64)) && continue
        tag = label_of(backend, threads)
        stat!("P1", "peak mean throughput, $tag, Float32, Ntot=$DEFAULT_NTOT", round(maximum(q.mean_cyps for q in r32); sigdigits=4), "column-years/s")
        stat!("P1", "Float32 / Float64 peak throughput, $tag", round(maximum(q.mean_cyps for q in r32) /
            maximum(q.mean_cyps for q in r64); sigdigits=3), "ratio")
        for c in intersect([q.columns for q in r32], [q.columns for q in r64])
            a = only(filter(q -> q.columns == c, r64)); b = only(filter(q -> q.columns == c, r32))
            stat!("P1", "Float32 speed-up at $c columns, $tag", round(a.mean_spy / b.mean_spy; sigdigits=3), "ratio")
            a.years == b.years && stat!("P1", "Float32 - Float64 mean melt at $c columns, $tag", round(b.melt - a.melt; sigdigits=3), "kg m-2 yr-1")
        end
    end
    for (backend, threads, _) in SERIES, ntot in (DEFAULT_NTOT, 30)
        r = series(s, backend, threads, ntot)
        isempty(r) && continue
        best = r[argmax([q.mean_cyps for q in r])]
        tag = "$(label_of(backend, threads)), Ntot=$ntot"
        stat!("P1", "peak mean throughput, $tag", round(best.mean_cyps; sigdigits=4), "column-years/s")
        stat!("P1", "columns at peak, $tag", best.columns, "columns")
    end
    gpu = series(s, "gpu", 8, DEFAULT_NTOT)
    for (_, threads, _) in SERIES[2:end]
        cpu = series(s, "threads", threads, DEFAULT_NTOT)
        (isempty(cpu) || isempty(gpu)) && continue
        stat!("P1", "GPU peak / CPU $threads-thread peak", round(maximum(q.mean_cyps for q in gpu) /
            maximum(q.mean_cyps for q in cpu); sigdigits=3), "ratio")
        for c in intersect([q.columns for q in gpu], [q.columns for q in cpu])
            g = only(filter(q -> q.columns == c, gpu)); q = only(filter(q -> q.columns == c, cpu))
            stat!("P1", "GPU / CPU $threads-thread speed-up at $c columns", round(q.mean_spy / g.mean_spy; sigdigits=3), "ratio")
        end
    end
    block_path = joinpath(POUT, "time_block_raw.csv")
    if isfile(block_path)
        b = summarize(read_table(block_path))
        for r in b
            base = filter(q -> q.backend == r.backend && q.threads == r.threads && q.columns == r.columns &&
                q.ntot == r.ntot, s)
            isempty(base) && continue
            stat!("P1", "speed-up from 365-step kernel blocking, $(label_of(r.backend, r.threads)), $(r.columns) columns",
                round(only(base).mean_spy / r.mean_spy; sigdigits=3), "ratio")
        end
    end

    # CPU runs at column counts that are multiples of the KA CPU workgroup size (16), which trigger
    # heap allocation in the threaded kernel (README, caveats); compared with the +1 counts.
    m16_path = joinpath(POUT, "scaling_raw_cpu_multiple_of_16.csv")
    if isfile(m16_path)
        for r in summarize(read_table(m16_path))
            base = filter(q -> q.backend == r.backend && q.threads == r.threads && q.columns == r.columns + 1 &&
                q.ntot == r.ntot, s)
            isempty(base) && continue
            stat!("P1", "wall time per year at $(r.columns + 1) columns / at $(r.columns) columns (multiple of 16), " *
                label_of(r.backend, r.threads), round(only(base).mean_spy / r.mean_spy; sigdigits=3), "ratio")
        end
    end

    gs = greenland_summary(read_table(joinpath(POUT, "greenland_raw.csv")))
    g32_path = joinpath(POUT, "greenland_raw_float32.csv")
    gs32 = isfile(g32_path) ? greenland_summary(read_table(g32_path); precision="Float32") : NamedTuple[]
    # Chion with one daily step and the original Fortran BESSI (bessi_fortran/run_bessi_fortran_benchmark.sh)
    s1_path = joinpath(POUT, "greenland_raw_substeps1.csv")
    gs1 = isfile(s1_path) ? greenland_summary(read_table(s1_path)) : NamedTuple[]
    f_path = joinpath(POUT, "greenland_raw_bessi_fortran.csv")
    ft = isfile(f_path) ? read_table(f_path) : nothing
    # serial build for 1 thread, OpenMP build for more
    fortran = isnothing(ft) ? NamedTuple[] : sort(filter(r -> (r.backend == "bessi_fortran_serial") == (r.threads == 1),
        greenland_summary(ft)); by=r -> r.threads)
    # the Fortran BESSI timings go to stats.csv only, not into the figure
    comparison = isempty(gs1) ? () : (gs1, NamedTuple[])
    greenland_panels = fig_greenland(sort([gs; gs32]; by=r -> (r.backend != "gpu", r.threads, r.precision != "Float64")), comparison)
    fig_performance(scaling_panels, greenland_panels)
    for r in gs1
        tag = label_of(r.backend, r.threads)
        stat!("P2", "Greenland one-year wall time, one daily step, mean, $tag", round(r.mean_total; digits=3), "s")
        base = filter(q -> q.backend == r.backend && q.threads == r.threads, gs)
        isempty(base) || stat!("P2", "Greenland wall time default / one daily step, $tag",
            round(only(base).mean_total / r.mean_total; sigdigits=3), "ratio")
    end
    isempty(gs1) || stat!("P2", "Greenland fraction of column-days with diurnal substeps, one daily step", round(gs1[1].substeps; digits=3), "fraction")
    for r in fortran
        tag = "BESSI Fortran, $(label_of("threads", r.threads))"
        stat!("P2", "Greenland one-year wall time, mean, $tag", round(r.mean_total; digits=3), "s")
        k = rows_where(ft, i -> ft[:backend][i] == r.backend && ft[:threads][i] == r.threads && ft[:phase][i] == "one_year_spun_up")
        isempty(k) || stat!("P2", "Greenland one-year wall time after 30 years, mean, $tag", round(mean(ft[:integration_s][k]); digits=3), "s")
        stat!("P2", "Greenland forcing I/O for one year, $tag", round(r.io; digits=3), "s")
        for (name, rows) in (("Chion", gs), ("Chion one daily step", gs1))
            q = filter(q -> q.backend == "threads" && q.threads == r.threads && q.precision == "Float64", rows)
            isempty(q) || stat!("P2", "$name / BESSI Fortran wall time, $(label_of("threads", r.threads))",
                round(only(q).mean_total / r.mean_total; sigdigits=3), "ratio")
        end
        g = filter(q -> q.backend == "gpu" && q.precision == "Float64", gs)
        isempty(g) || stat!("P2", "$tag / Chion GPU wall time", round(r.mean_total / only(g).mean_total; sigdigits=3), "ratio")
    end
    # after 30 years from an empty column (filled firn column), Chion and BESSI
    spun(path) = isfile(path) ? greenland_summary(read_table(path); phase="one_year_spun_up") : NamedTuple[]
    fs = isnothing(ft) ? NamedTuple[] : filter(r -> (r.backend == "bessi_fortran_serial") == (r.threads == 1),
        greenland_summary(ft; phase="one_year_spun_up"))
    for (name, rows) in (("Chion", spun(joinpath(POUT, "greenland_raw_spunup.csv"))),
                         ("Chion one daily step", spun(joinpath(POUT, "greenland_raw_spunup_substeps1.csv"))))
        for r in rows
            tag = label_of(r.backend, r.threads)
            stat!("P2", "Greenland one-year wall time after 30 years, mean, $name, $tag", round(r.mean_total; digits=3), "s")
            f = filter(q -> q.threads == r.threads, fs)
            r.backend == "gpu" || isempty(f) || stat!("P2", "$name / BESSI Fortran wall time after 30 years, $tag",
                round(r.mean_total / only(f).mean_total; sigdigits=3), "ratio")
        end
    end
    if !isnothing(ft)
        for b in unique(ft[:backend]), t in unique(ft[:threads])
            k = rows_where(ft, i -> ft[:backend][i] == b && ft[:threads][i] == t)
            isempty(k) && continue
            ref = rows_where(ft, i -> ft[:backend][i] == "bessi_fortran_serial")
            same = ft[:smb_checksum][k] == ft[:smb_checksum][ref] && ft[:snowmass_checksum][k] == ft[:snowmass_checksum][ref]
            stat!("P2", "BESSI Fortran $b $(Int(t)) threads identical to serial (SMB and mass checksums)", same ? 1 : 0, "bool")
        end
    end
    for r in gs32
        tag = label_of(r.backend, r.threads)
        stat!("P2", "Greenland one-year wall time (setup+integration), mean, $tag, Float32", round(r.mean_total; digits=3), "s")
        base = filter(q -> q.backend == r.backend && q.threads == r.threads, gs)
        isempty(base) || stat!("P2", "Greenland Float32 speed-up, $tag", round(only(base).mean_total / r.mean_total; sigdigits=3), "ratio")
    end
    for r in gs
        tag = label_of(r.backend, r.threads)
        stat!("P2", "Greenland one-year wall time (setup+integration), mean, $tag", round(r.mean_total; digits=3), "s")
        stat!("P2", "Greenland one-year wall time (setup+integration), min, $tag", round(r.min_total; digits=3), "s")
        stat!("P2", "Greenland one-year integration time, mean, $tag", round(r.mean_integration; digits=3), "s")
        stat!("P2", "Greenland first year incl. compilation, $tag", round(r.warmup; digits=2), "s")
        stat!("P2", "Greenland forcing I/O for one year, $tag job", round(r.io; digits=2), "s")
    end
    g = filter(r -> r.backend == "gpu", gs)
    for r in filter(r -> r.backend != "gpu", gs)
        isempty(g) || stat!("P2", "Greenland GPU speed-up vs $(label_of(r.backend, r.threads))",
            round(r.mean_total / only(g).mean_total; sigdigits=3), "ratio")
    end
    stat!("P2", "Greenland fraction of column-days with diurnal substeps", round(gs[1].substeps; digits=3), "fraction")
    open(joinpath(POUT, "stats.csv"), "w") do io
        println(io, "figure,quantity,value,unit")
        for q in STATS
            println(io, join((q.figure, '"' * q.quantity * '"', q.value, '"' * q.unit * '"'), ','))
        end
    end
    println("Wrote figures and stats to $POUT")
end

main()
