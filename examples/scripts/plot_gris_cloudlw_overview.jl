#!/usr/bin/env julia

# Overview of the calibrated fully parameterized setup (cloud-proxy longwave,
# SEMIX sensible heat 2.5/40, dynamic albedo 0.81/0.70/0.40) against MAR, the
# prescribed-flux runs and the uncalibrated parameterized runs.
# Writes figures and README.md to output/gris_cloudlw_parameterized/overview.
using Plots, Printf, Statistics

const ROOT = joinpath(@__DIR__, "..", "..", "output")
const NEW = joinpath(ROOT, "gris_cloudlw_parameterized")
const SUITE = joinpath(ROOT, "gris_substrate5_suite")
const OUT = joinpath(NEW, "overview")
const VARS = (:surface_smb => "Surface SMB", :melt => "Melt", :runoff => "Runoff", :refreezing => "Refreezing")
const LAST_COMPLETE_YEAR = 2025

const TRANSIENTS = [
    ("Prescribed fluxes", joinpath(SUITE, "transient_daily_prescribed", "steps_8"), "#2a78d6"),
    ("Parameterized, uncalibrated", joinpath(SUITE, "transient_daily_parameterized", "steps_8"), "#9a9a9a"),
    ("Parameterized, calibrated", joinpath(NEW, "transient_daily", "steps_8"), "#eb6834"),
]
const CALIBRATED_STEPS = [(s, joinpath(NEW, "transient_daily", "steps_$s")) for s in (4, 8, 24)]
const EQUILIBRIA = [
    ("Cycled 1940–1979, prescribed", joinpath(NEW, "equilibrium_daily_cycled_prescribed", "steps_8", "equilibrium_annual_integrated_budgets.csv")),
    ("Cycled 1940–1979, parameterized (calibrated)", joinpath(NEW, "equilibrium_daily_cycled_parameterized", "steps_8", "equilibrium_annual_integrated_budgets.csv")),
    ("Climatology, parameterized (calibrated)", joinpath(NEW, "equilibrium_daily", "steps_8", "parameterized_constant_2k_steps_8_ntot_15_annual_integrated_budgets.csv")),
]

"""Minimal CSV reader: Dict of column name => vector (Float64 where parseable)."""
function read_table(path)
    lines = filter(!isempty, readlines(path))
    header = Symbol.(split(lines[1], ','))
    cells = [split(l, ',') for l in lines[2:end]]
    column(i) = (raw = [c[i] for c in cells]; parsed = tryparse.(Float64, raw); any(isnothing, parsed) ? String.(raw) : Float64.(parsed))
    return Dict(h => column(i) for (i, h) in enumerate(header))
end

function yearly(dir)
    t = read_table(joinpath(dir, "yearly_integrated_rates.csv"))
    keep = t[:year] .<= LAST_COMPLETE_YEAR
    return Dict(k => v[keep] for (k, v) in t)
end
mar_col(v) = Symbol("mar_$(v)_gt_per_year")
chion_col(v) = Symbol("chion_prescribed_$(v)_gt_per_year")
window(df, a, b) = (keep = a .<= df[:year] .<= b; Dict(k => v[keep] for (k, v) in df))
mean_bias(df, v) = mean(df[chion_col(v)] .- df[mar_col(v)])

function transient_figure(path)
    panels = Any[]
    reference = yearly(TRANSIENTS[1][2])
    for (v, label) in VARS
        p = plot(reference[:year], reference[mar_col(v)]; c=:black, lw=2.5, label="MAR", title=label,
            xlabel="Year", ylabel="Gt yr⁻¹", framestyle=:box, grid=:y, gridalpha=0.15,
            legend=v == :surface_smb ? :bottomleft : false, xlims=(1940, LAST_COMPLETE_YEAR))
        for (name, dir, colour) in TRANSIENTS
            df = yearly(dir)
            plot!(p, df[:year], df[chion_col(v)]; c=colour, lw=1.6, label=name)
        end
        push!(panels, p)
    end
    savefig(plot(panels...; layout=(2, 2), size=(1300, 850), left_margin=5Plots.mm, bottom_margin=4Plots.mm,
            plot_title="Daily MAR transient 1940–2025, 8 diurnal substeps", plot_titlefontsize=12), path)
end

function calibrated_bias_figure(path)
    panels = Any[]
    for (v, label) in VARS
        p = plot(; title=label, xlabel="Year", ylabel="Chion − MAR [Gt yr⁻¹]", framestyle=:box, grid=:y, gridalpha=0.15,
            legend=v == :surface_smb ? :topleft : false, xlims=(1940, LAST_COMPLETE_YEAR))
        for (s, dir) in CALIBRATED_STEPS
            df = yearly(dir)
            plot!(p, df[:year], df[chion_col(v)] .- df[mar_col(v)]; lw=1.6, label="$s substeps")
        end
        df = yearly(TRANSIENTS[1][2])
        plot!(p, df[:year], df[chion_col(v)] .- df[mar_col(v)]; c=:black, lw=1.2, ls=:dash, label="prescribed, 8 substeps")
        plot!(p, [1940, LAST_COMPLETE_YEAR], [0, 0]; c=:gray, lw=1, label="", primary=false)
        push!(panels, p)
    end
    savefig(plot(panels...; layout=(2, 2), size=(1300, 850), left_margin=5Plots.mm, bottom_margin=4Plots.mm,
            plot_title="Calibrated parameterized transient: yearly bias by diurnal substeps", plot_titlefontsize=12), path)
end

function equilibrium_rows()
    rows = []
    for (name, file) in EQUILIBRIA
        isfile(file) || continue
        df = read_table(file)
        get_bias(v) = (i = findfirst(==(String(v)), df[:variable]); isnothing(i) ? NaN : df[:Chion_minus_MAR_Gt][i])
        push!(rows, (name, Dict(v => get_bias(v) for (v, _) in VARS)))
    end
    return rows
end

function equilibrium_figure(path, rows)
    names = first.(rows)
    panels = Any[]
    for (v, label) in VARS
        values = [r[2][v] for r in rows]
        p = bar(1:length(rows), values; c=["#2a78d6", "#eb6834", "#c5a028"][1:length(rows)], label="", title=label,
            ylabel="Chion − MAR [Gt yr⁻¹]", xticks=(1:length(rows), ["cycled\nprescribed", "cycled\nparameterized", "climatology\nparameterized"][1:length(rows)]),
            framestyle=:box, grid=:y, gridalpha=0.15)
        hline!(p, [0]; c=:gray, lw=1, label="")
        push!(panels, p)
    end
    savefig(plot(panels...; layout=(2, 2), size=(1200, 850), left_margin=5Plots.mm,
            plot_title="Equilibrium bias against MAR (final cycle / final year), 8 substeps", plot_titlefontsize=12), path)
end

function stats_line(df, v)
    m, c = df[mar_col(v)], df[chion_col(v)]
    x = df[:year] .- mean(df[:year])
    slope(y) = sum(x .* (y .- mean(y))) / sum(x .^ 2) * 10
    return @sprintf("%.2f | %.1f | %.0f / %.0f | %+.1f / %+.1f", cor(m, c), sqrt(mean((c .- m) .^ 2)), std(m), std(c), slope(m), slope(c))
end

function main()
    mkpath(OUT)
    transient_figure(joinpath(OUT, "transient_yearly_budgets.pdf"))
    calibrated_bias_figure(joinpath(OUT, "calibrated_yearly_bias_by_substeps.pdf"))
    eq = equilibrium_rows()
    equilibrium_figure(joinpath(OUT, "equilibrium_bias.pdf"), eq)
    open(joinpath(OUT, "README.md"), "w") do io
        println(io, """# Calibrated fully parameterized surface setup

All runs: graybody-snow SEMIX SEB, 5-layer ice substrate, fine near-surface layers, 15 snow layers, daily MAR forcing with 1 K reconstructed diurnal temperature cycle. The calibrated parameterized setup adds no forcing variable:

- longwave down `longwave_scheme=:cloud_proxy`: ε = 0.624 + 0.0032 (Tₐ − T₀) + 0.613 n, with cloudiness n = 1 − τ/τ_clear from the daily shortwave transmissivity τ = SW↓/SW_TOA, τ_clear = 0.85 + 0.075 km⁻¹ × elevation; n = 0.389 in the polar night;
- SEMIX sensible heat: exchange factor 2.5, stable coefficient 40, constant 5 m s⁻¹ wind;
- dynamic albedo: dry 0.81, wet 0.70, ice 0.40, aging scale 1.0.

Calibrated in free coupled runs 1980–1989 (`calibrate_gris_surface_parameterization.jl`, `CHION_CALIBRATION=diagnose`, sets `sensible`, `sensible_wide`, `albedo`); launcher `run_gris_cloudlw_parameterized_gpu.sh`. Transient metrics exclude the incomplete 2026.

## Transient, mean Chion − MAR (Gt/yr)
""")
        println(io, "| run | period | melt | runoff | refreezing | surface SMB |\n|---|---|---:|---:|---:|---:|")
        runs = vcat([(n, d) for (n, d, _) in TRANSIENTS], [("Parameterized, calibrated, $s substeps", d) for (s, d) in CALIBRATED_STEPS if s != 8])
        for (name, dir) in runs, (a, b) in ((1940, 1979), (1980, 2017))
            df = window(yearly(dir), a, b)
            println(io, @sprintf("| %s | %d–%d | %+.1f | %+.1f | %+.1f | %+.1f |", name, a, b,
                mean_bias(df, :melt), mean_bias(df, :runoff), mean_bias(df, :refreezing), mean_bias(df, :surface_smb)))
        end
        println(io, "\n## Interannual agreement 1940–2025, 8 substeps\n\ncorrelation | RMSE (Gt/yr) | std MAR / Chion | trend MAR / Chion (Gt/yr per decade)\n")
        println(io, "| run | variable | r | RMSE | std | trend |\n|---|---|---:|---:|---|---|")
        for (name, dir, _) in TRANSIENTS[[1, 3]], (v, label) in VARS[1:3]
            parts = split(stats_line(yearly(dir), v), " | ")
            println(io, "| $name | $label | ", join(parts, " | "), " |")
        end
        println(io, "\n## Equilibrium, Chion − MAR (Gt/yr)\n\nCycled: real daily MAR 1940–1979 repeated 5× from the initial state, final cycle against MAR's 1940–1979 mean (`run_gris_daily_cycled_equilibrium.jl`). Climatology: final year after 200 y on the 1940–1980 daily climatology, whose near-daily snowfall keeps the dynamic albedo fresh.\n")
        println(io, "| run | melt | runoff | refreezing | surface SMB |\n|---|---:|---:|---:|---:|")
        for (name, b) in eq
            cell(x) = isnan(x) ? "n/a" : @sprintf("%+.1f", x)
            println(io, "| $name | ", join(cell.((b[:melt], b[:runoff], b[:refreezing], b[:surface_smb])), " | "), " |")
        end
        println(io, """

## Caveats

- Above 1000 m the parameterized summer sensible heat stays 3–5 W/m² below MAR for every SEMIX setting; the calibrated albedo compensates with 10–18 W/m² more absorbed shortwave there, so integrated budgets agree but individual fluxes do not.
- The climatology equilibrium under-melts with any dynamic albedo; use the cycled equilibrium.

## Figures

- `transient_yearly_budgets.pdf` — MAR vs prescribed, uncalibrated and calibrated parameterized (8 substeps)
- `calibrated_yearly_bias_by_substeps.pdf` — yearly bias of the calibrated setup at 4/8/24 substeps
- `equilibrium_bias.pdf` — cycled and climatology equilibria
- per run: `transient_daily/steps_*/yearly_integrated_rates.pdf`, `equilibrium_daily_cycled_*/steps_8/{cycle_drift,final_cycle_budgets}.pdf`
""")
    end
    println("Overview written to ", OUT)
end

main()
