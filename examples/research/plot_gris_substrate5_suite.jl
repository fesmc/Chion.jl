#!/usr/bin/env julia
# Overview of the 5-layer ice-substrate experiment suite (run_gris_substrate5_suite_gpu.sh):
# equilibrium and transient runs with daily prescribed, daily parameterized,
# daily parameterized with MAR albedo, and native 3-hourly MAR forcing.
import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using NCDatasets, Plots, Printf, Statistics, Dates

const SUITE = get(ENV, "CHION_SUITE_DIR", joinpath(@__DIR__, "..", "..", "output", "gris_substrate5_suite"))
const DEST = joinpath(SUITE, "overview")
const DAILY_DIR = "/p/projects/ou/labs/ai/Nils/MAR_daily/mar_daily"
const STEPS = (1, 4, 8, 24)
const COLORS = Dict(:prescribed => "#2a78d6", :native => "#eb6834", :parameterized => "#1baf7a",
    :parameterized_mar_albedo => "#eda100")
const STEP_RAMP = Dict(1 => "#86b6ef", 4 => "#3987e5", 8 => "#1c5cab", 24 => "#0d366b")
const PARAM_RAMP = Dict(1 => "#8fd9bd", 4 => "#3fbf8f", 8 => "#1a8f63", 24 => "#0f5e40")
const ALBEDO_RAMP = Dict(1 => "#f7cf6e", 4 => "#eda100", 8 => "#b07800", 24 => "#6e4b00")
# Daily forcing types: (kind, equilibrium case family, label, yearly-plot ramp).
const DAILY_KINDS = (
    (:prescribed, "fully_prescribed", "prescribed fluxes", STEP_RAMP),
    (:parameterized, "parameterized_constant_2k", "parameterized fluxes", PARAM_RAMP),
    (:parameterized_mar_albedo, "parameterized_mar_albedo", "parameterized fluxes, MAR albedo", ALBEDO_RAMP),
)
const KIND_INFO = Dict(k[1] => k for k in DAILY_KINDS)
daily_kinds() = first.(DAILY_KINDS)
const VARS = (("melt", "Melt"), ("runoff", "Runoff"), ("refreezing", "Refreezing"), ("surface_smb", "Surface SMB"))
const WINDOW = (1980, 2017)
const CATEGORIES = [string.(STEPS)..., "native 3h"]

function read_csv(path)
    lines = filter(!isempty, readlines(path))
    header = split(lines[1], ',')
    [Dict(zip(header, split(l, ','))) for l in lines[2:end]]
end
num(x) = parse(Float64, x)
style(; kw...) = (framestyle=:box, grid=:y, gridalpha=0.15, titlefontsize=11, kw...)
zero_line!(p, xs) = plot!(p, collect(xs), [0, 0]; c=:gray, lw=1, label="", primary=false)

# ---------- equilibrium ----------
function equilibrium_file(kind, step=nothing)
    kind == :native && return joinpath(SUITE, "equilibrium_native_3hourly_cycled", "equilibrium_annual_integrated_budgets.csv")
    case = KIND_INFO[kind][2]
    joinpath(SUITE, "equilibrium_daily_$(kind)", "steps_$step", "$(case)_steps_$(step)_ntot_15_annual_integrated_budgets.csv")
end
function equilibrium_bias(kind, step, var)
    file = equilibrium_file(kind, step)
    isfile(file) || return NaN
    rows = filter(r -> r["variable"] == var, read_csv(file))
    isempty(rows) ? NaN : num(only(rows)["Chion_minus_MAR_Gt"])
end

# ---------- transients: yearly (years, MAR, Chion) ----------
function daily_series(kind, step)
    file = joinpath(SUITE, "transient_daily_$(kind)", "steps_$step", "yearly_integrated_rates.csv")
    isfile(file) || return nothing
    # The 2026 forcing file is incomplete, so its annual budget is not comparable.
    rows = filter(r -> parse(Int, r["year"]) <= 2025, read_csv(file))
    years = [parse(Int, r["year"]) for r in rows]
    Dict(v => (years, num.(getindex.(rows, "mar_$(v)_gt_per_year")), num.(getindex.(rows, "chion_prescribed_$(v)_gt_per_year"))) for (v, _) in VARS)
end
function native_series()
    file = joinpath(SUITE, "transient_native_3hourly", "native_3hourly_transient_1980_2017_budgets.csv")
    isfile(file) || return nothing
    rows = read_csv(file)
    years = [parse(Int, r["year"]) for r in rows]
    Dict(v => (years, num.(getindex.(rows, "mar_$(v)_gt")), num.(getindex.(rows, "chion_$(v)_gt"))) for (v, _) in VARS)
end
function window_bias(series, var, window=WINDOW)
    series === nothing && return NaN
    years, mar, chion = series[var]
    k = window[1] .<= years .<= window[2]
    any(k) ? mean(chion[k] .- mar[k]) : NaN
end

function category_panel(title, values_by_kind; ylabel="Chion − MAR [Gt yr⁻¹]", legend=false)
    x = 1:length(CATEGORIES)
    p = plot(; xticks=(x, CATEGORIES), xlabel="Daily substeps / native", ylabel, title, legend,
        xlims=(0.5, length(x) + 0.5), legendfontsize=8, style()...)
    zero_line!(p, (0.5, length(x) + 0.5))
    vline!(p, [length(STEPS) + 0.5]; c=:gray, ls=:dot, lw=1, label="", primary=false)
    for kind in daily_kinds()
        vals = values_by_kind[kind]
        all(isnan, vals) && continue
        label = "daily, " * KIND_INFO[kind][3]
        plot!(p, x[1:length(STEPS)], vals; c=COLORS[kind], lw=2, marker=:circle, ms=6, msw=0, label)
    end
    scatter!(p, [x[end]], [values_by_kind[:native]]; c=COLORS[:native], marker=:diamond, ms=9, msw=0,
        label="native 3-hourly, prescribed fluxes")
    return p
end

function equilibrium_figure(io)
    println(io, "\n## Equilibrium (Chion − MAR, Gt/yr)\n")
    println(io, "Daily: final year after 200-y spin-up on the 1940-1980 daily climatology. ",
        "Native: mean of the final cycle of native 1980-1989 (after the same climatology spin-up), against MAR's 1980-1989 mean.\n")
    println(io, "| run | melt | runoff | refreezing |\n|---|---:|---:|---:|")
    for kind in daily_kinds(), s in STEPS
        println(io, @sprintf("| daily %s, %d substep%s | %+.1f | %+.1f | %+.1f |", KIND_INFO[kind][3], s, s == 1 ? "" : "s",
            (equilibrium_bias(kind, s, v) for v in ("melt", "runoff", "refreezing"))...))
    end
    println(io, @sprintf("| native 3-hourly, cycled 1980-1989 | %+.1f | %+.1f | %+.1f |",
        (equilibrium_bias(:native, nothing, v) for v in ("melt", "runoff", "refreezing"))...))
    panels = Any[]
    for (var, title) in VARS[1:3]
        vals = Dict{Symbol,Any}(kind => [equilibrium_bias(kind, s, var) for s in STEPS] for kind in daily_kinds())
        vals[:native] = equilibrium_bias(:native, nothing, var)
        push!(panels, category_panel(title, vals; legend=var == "melt" ? :bottomright : false))
    end
    savefig(plot(panels...; layout=(1, 3), size=(1650, 520), bottom_margin=8Plots.mm, left_margin=6Plots.mm, top_margin=6Plots.mm,
            plot_title="Equilibrium bias", plot_titlefontsize=12), joinpath(DEST, "equilibrium_bias.pdf"))
end

function transient_mean_figure(io)
    series = Dict((kind, s) => daily_series(kind, s) for kind in daily_kinds() for s in STEPS)
    native = native_series()
    for window in (WINDOW, (1940, 1980))
        println(io, "\n## Transients: mean Chion − MAR over $(window[1])-$(window[2]) (Gt/yr)\n")
        println(io, "| run | melt | runoff | refreezing | surface SMB |\n|---|---:|---:|---:|---:|")
        for kind in daily_kinds(), s in STEPS
            println(io, @sprintf("| daily %s, %d substep%s | %+.1f | %+.1f | %+.1f | %+.1f |", KIND_INFO[kind][3], s, s == 1 ? "" : "s",
                (window_bias(series[(kind, s)], v, window) for (v, _) in VARS)...))
        end
        window == WINDOW && println(io, @sprintf("| native 3-hourly | %+.1f | %+.1f | %+.1f | %+.1f |",
            (window_bias(native, v, window) for (v, _) in VARS)...))
    end
    panels = Any[]
    for (var, title) in VARS
        vals = Dict{Symbol,Any}(kind => [window_bias(series[(kind, s)], var) for s in STEPS] for kind in daily_kinds())
        vals[:native] = window_bias(native, var)
        push!(panels, category_panel(title, vals; legend=var == "melt" ? :bottomright : false))
    end
    savefig(plot(panels...; layout=(2, 2), size=(1400, 950), left_margin=5Plots.mm, top_margin=4Plots.mm,
            plot_title="Transients: mean bias $(WINDOW[1])–$(WINDOW[2])", plot_titlefontsize=12),
        joinpath(DEST, "transient_mean_bias_$(WINDOW[1])_$(WINDOW[2]).pdf"))
    return series, native
end

function transient_yearly_figures(series, native)
    for (kind, _, label, ramp) in DAILY_KINDS
        panels = Any[]
        for (var, title) in VARS
            p = plot(; xlabel="Year", ylabel="Gt yr⁻¹", title, legend=var == "melt" ? :topleft : false, legendfontsize=7, style()...)
            ref = series[(kind, 8)] === nothing ? nothing : series[(kind, 8)][var]
            ref === nothing || plot!(p, ref[1], ref[2]; c=:black, lw=2.5, label="MAR")
            for s in STEPS
                d = series[(kind, s)]
                d === nothing && continue
                plot!(p, d[var][1], d[var][3]; c=ramp[s], lw=1.5, label="Chion daily $label, $s substep$(s == 1 ? "" : "s")")
            end
            native === nothing || plot!(p, native[var][1], native[var][3]; c=COLORS[:native], lw=1.8,
                label="Chion native 3-hourly, prescribed")
            push!(panels, p)
        end
        savefig(plot(panels...; layout=(2, 2), size=(1500, 950), left_margin=5Plots.mm, top_margin=4Plots.mm,
                plot_title="Transients, daily forcing, $label", plot_titlefontsize=12),
            joinpath(DEST, "transient_yearly_budgets_$(kind).pdf"))
    end
end

# ---------- mean seasonal cycle 1980-2017 ----------
function geometry()
    NCDataset(joinpath(DAILY_DIR, "MARv3.14.3-10km-daily-ERA5-1980.nc")) do d
        Float64.(d["MSK"][:, :]) .>= 50, Float64.(d["AREA"][:, :])
    end
end
gt(field, mask, area) = sum(coalesce.(field, NaN)[mask] .* area[mask]) * 1e-6

function mar_seasonal(mask, area)
    names = Dict("melt" => "ME", "runoff" => "RU", "refreezing" => "RZ", "surface_smb" => "SMB")
    out = Dict(v => zeros(12) for (v, _) in VARS)
    for year in WINDOW[1]:WINDOW[2]
        NCDataset(joinpath(DAILY_DIR, "MARv3.14.3-10km-daily-ERA5-$(year).nc")) do d
            for (v, _) in VARS
                var = d[names[v]].var
                A = Float64.(ndims(var) == 4 ? var[:, :, 1, :] : var[:, :, :])
                for k in axes(A, 3)
                    out[v][month(Date(year, 1, 1) + Day(k - 1))] += gt(A[:, :, k], mask, area)
                end
            end
        end
    end
    Dict(v => x ./ (WINDOW[2] - WINDOW[1] + 1) for (v, x) in out)
end

function chion_seasonal(files_and_records, mask, area)
    out = Dict(v => zeros(12) for (v, _) in VARS)
    n = 0
    for (file, records) in files_and_records
        NCDataset(file) do d
            times = d["t"][:]
            for r in records
                m = month(times[r])
                for (v, _) in VARS
                    out[v][m] += gt(Float64.(coalesce.(d[v][:, :, r], NaN32)), mask, area)
                end
            end
        end
    end
    Dict(v => x ./ (WINDOW[2] - WINDOW[1] + 1) for (v, x) in out)
end

function daily_monthly_records(kind, step)
    file = joinpath(SUITE, "transient_daily_$(kind)", "steps_$step", "gris_$(kind)_monthly_1940_2026.nc")
    isfile(file) || return nothing
    records = NCDataset(file) do d
        findall(t -> WINDOW[1] <= year(t) <= WINDOW[2], d["t"][:])
    end
    [(file, records)]
end
function native_monthly_records()
    dir = joinpath(SUITE, "transient_native_3hourly", "monthly")
    files = [joinpath(dir, "native_$(y)_monthly.nc") for y in WINDOW[1]:WINDOW[2]]
    all(isfile, files) || return nothing
    [(f, 1:12) for f in files]
end

function seasonal_figure()
    mask, area = geometry()
    mar = mar_seasonal(mask, area)
    runs = Any[]
    for kind in daily_kinds()
        label = "daily, $(KIND_INFO[kind][3]), 8 substeps"
        fr = daily_monthly_records(kind, 8)
        fr === nothing || push!(runs, (label, chion_seasonal(fr, mask, area), COLORS[kind]))
    end
    nr = native_monthly_records()
    nr === nothing || push!(runs, ("native 3-hourly, prescribed", chion_seasonal(nr, mask, area), COLORS[:native]))
    panels = Any[]
    for (var, title) in VARS
        p = plot(; xlabel="Month", ylabel="Gt month⁻¹", xticks=1:12, title, legend=var == "melt" ? :topleft : false,
            legendfontsize=8, style()...)
        plot!(p, 1:12, mar[var]; c=:black, lw=2.5, ls=:dash, label="MAR")
        for (label, values, color) in runs
            plot!(p, 1:12, values[var]; c=color, lw=2, marker=:circle, ms=3, msw=0, label="Chion: $label")
        end
        push!(panels, p)
    end
    savefig(plot(panels...; layout=(2, 2), size=(1300, 900), left_margin=5Plots.mm, top_margin=4Plots.mm,
            plot_title="Mean seasonal cycle $(WINDOW[1])–$(WINDOW[2])", plot_titlefontsize=12),
        joinpath(DEST, "seasonal_cycle_$(WINDOW[1])_$(WINDOW[2]).pdf"))
end

function main()
    mkpath(DEST)
    open(joinpath(DEST, "README.md"), "w") do io
        println(io, "# Ice-substrate experiment suite\n")
        println(io, "All runs: graybody LW (`seb_scheme=:semix`), 5-layer ice substrate (0.05 m top, 1.55 m total), ",
            "fine near-surface layers (0.02/0.05/0.10/0.30 m), no deep-layer cap, 15 snow layers, 200-y spin-up, ",
            "1 K reconstructed diurnal air-temperature cycle for daily forcing. Prescribed runs use MAR SW net, LWD, SHF, LHF; ",
            "parameterized runs compute them (dynamic albedo, constant 5 m/s wind); the MAR-albedo runs use the same ",
            "parameterized SEB with MAR's AL2 albedo prescribed. Launcher: `run_gris_substrate5_suite_gpu.sh`.")
        equilibrium_figure(io)
        series, native = transient_mean_figure(io)
        transient_yearly_figures(series, native)
        println(io, "\n## Figures\n")
        for f in ("equilibrium_bias.pdf", "transient_mean_bias_$(WINDOW[1])_$(WINDOW[2]).pdf",
                  ("transient_yearly_budgets_$(k).pdf" for k in daily_kinds())..., "seasonal_cycle_$(WINDOW[1])_$(WINDOW[2]).pdf")
            println(io, "- `$f`")
        end
        println(io, "\nPer-run plots: `equilibrium_daily_*/steps_*/` (full equilibrium diagnostics), `transient_daily_*/steps_*/yearly_integrated_rates.pdf`, ",
            "`transient_native_3hourly/native_3hourly_transient_1980_2017_budgets.pdf`, `equilibrium_native_3hourly_cycled/{cycle_drift,final_cycle_budgets}.pdf`.")
    end
    seasonal_figure()
    println("Overview written to ", DEST)
end

main()
