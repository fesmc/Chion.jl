#!/usr/bin/env julia
# Fit of the cloud-proxy longwave parameterisation to daily MAR over the Greenland ice sheet,
# and the figure documenting it (L1_cloud_proxy_longwave).
#
# Three stages, each fitted on 1980–1989 and evaluated on 2010–2019 (every 7th ice cell):
#   1. clear-sky transmissivity  τcs = τ0 + τz z + τS ln(SW_TOA / 300 W m⁻²)     against MAR SWDC/SW_TOA
#   2. cloud transmission        1 − τ/τcs = a nᵖ exp(−b z)                      against MAR cloud cover
#      (maximum overlap of low, middle and high cloud), inverted to a cloud fraction n̂ ∈ [0, 1];
#      in the polar night n̂ is the constant night cloud fraction (fitted in stage 3)
#   3. emissivity                ε = εcs(T)(1 − n̂^q) + εoc(T) n̂^q, εcs/oc linear in T − T0
#                                                                               against MAR LWD
# The functional forms follow Konzelmann et al. (1994). The sample is cached in OUT/mar_sample.jls.
#
# Usage: julia --project=. examples/research/fit_cloud_proxy_longwave.jl
#   CHION_FORCING_DIR  MAR daily files (default /p/projects/ou/labs/ai/Nils/MAR_daily/mar_daily)
#   CHION_OUTPUT_DIR   output directory (default output/longwave_fit)
using NCDatasets, Serialization, Statistics, Printf, LinearAlgebra, Chion
include(joinpath(@__DIR__, "..", "figures", "figure_style.jl"))

const FORCING_DIR = get(ENV, "CHION_FORCING_DIR", "/p/projects/ou/labs/ai/Nils/MAR_daily/mar_daily")
const OUT = get(ENV, "CHION_OUTPUT_DIR", joinpath(@__DIR__, "..", "..", "output", "longwave_fit"))
const CAL_YEARS, VAL_YEARS = 1980:1989, 2010:2019
const THIN = 7
const σSB = 5.670373e-8
mkpath(OUT)

# ---------- MAR sample ----------
function extract_sample()
    path = joinpath(OUT, "mar_sample.jls")
    isfile(path) && return deserialize(path)
    keys = (:year, :doy, :month, :lat, :z, :T, :lwd, :swd, :swdc, :cloud, :melt, :toa)
    d = Dict(k => Float64[] for k in keys)
    for year in vcat(CAL_YEARS, VAL_YEARS)
        ds = NCDataset(joinpath(FORCING_DIR, "MARv3.14.3-10km-daily-ERA5-$year.nc"))
        msk = ds["MSK"][:, :]; lat = ds["LAT"][:, :]; sh = ds["SH"][:, :]
        ice = findall(x -> !ismissing(x) && x > 50, msk)[1:THIN:end]
        lwd, swd, swdc = ds["LWD"][:, :, :], ds["SWD"][:, :, :], ds["SWDC"][:, :, :]
        cloud = max.(ds["CU"][:, :, :], ds["CM"][:, :, :], ds["CD"][:, :, :])
        tt = ds["TTZ"][:, :, 1, :]; me = ds["ME"][:, :, 1, :]; mm = ds["MM"][:]
        for t in axes(tt, 3)
            λ = Chion._solar_longitude_deg_from_calendar_day(t)
            for I in ice
                push!(d[:year], year); push!(d[:doy], t); push!(d[:month], mm[t]); push!(d[:lat], lat[I]); push!(d[:z], sh[I])
                push!(d[:T], tt[I, t] + 273.15); push!(d[:lwd], lwd[I, t]); push!(d[:swd], swd[I, t]); push!(d[:swdc], swdc[I, t])
                push!(d[:cloud], cloud[I, t]); push!(d[:melt], me[I, t])
                push!(d[:toa], Chion._daily_toa_shortwave(Float64(lat[I]), λ, t))
            end
        end
        close(ds)
        @info "MAR sample" year
    end
    serialize(path, d)
    return d
end

# ---------- fit ----------
rmse(x) = sqrt(mean(x .^ 2))
clear_sky(t, z_km, toa) = t[1] .+ t[2] .* z_km .+ t[3] .* log.(max.(toa, 50) ./ 300)
cloud_fraction(deficit, z_km, a, p, b) = clamp.(clamp.(deficit ./ (a .* exp.(-b .* z_km)), 0, 1) .^ (1 / p), 0, 1)
emissivity(c, ΔT, n, q) = (c[1] .+ c[2] .* ΔT) .* (1 .- n .^ q) .+ (c[3] .+ c[4] .* ΔT) .* n .^ q
function fit_emissivity(lwd, B, ΔT, n, q)
    w = n .^ q
    return hcat(1 .- w, ΔT .* (1 .- w), w, ΔT .* w) .* B \ lwd
end

d = extract_sample()
cal = in.(d[:year], Ref(CAL_YEARS)); val = in.(d[:year], Ref(VAL_YEARS))
sun = d[:toa] .> 50
z_km = d[:z] ./ 1000; ΔT = d[:T] .- 273.15; B = σSB .* d[:T] .^ 4
τ = d[:swd] ./ max.(d[:toa], 1e-6); τcs_mar = d[:swdc] ./ max.(d[:toa], 1e-6)
jja = in.(d[:month], Ref((6, 7, 8))); melt_day = d[:melt] .> 1

# stage 1
k1 = cal .& sun .& (τcs_mar .> 0.3)
t = hcat(ones(count(k1)), z_km[k1], log.(d[:toa][k1] ./ 300)) \ τcs_mar[k1]
τcs = min.(clear_sky(t, z_km, d[:toa]), 1)
# stage 2
deficit = clamp.(1 .- τ ./ τcs, 0, 1)
k2 = cal .& sun .& (d[:toa] .> 100)
best = (Inf, 0.0, 0.0, 0.0)
for p in 1.0:0.25:5.0, b in 0.0:0.025:0.6
    w = d[:cloud][k2] .^ p .* exp.(-b .* z_km[k2])
    a = dot(w, deficit[k2]) / dot(w, w)
    e = rmse(a .* w .- deficit[k2])
    e < best[1] && (global best = (e, a, p, b))
end
_, a, p, b = best
n_day = cloud_fraction(deficit, z_km, a, p, b)
# stage 3 (q and the night cloud fraction by grid search on a thinned calibration sample)
idx3 = findall(cal)[1:5:end]
best3 = (Inf, 0.0, 0.0)
for q in 1.0:0.5:5.0, nn in 0.5:0.05:1.0
    n = ifelse.(sun[idx3], n_day[idx3], nn)
    c = fit_emissivity(d[:lwd][idx3], B[idx3], ΔT[idx3], n, q)
    e = rmse(emissivity(c, ΔT[idx3], n, q) .* B[idx3] .- d[:lwd][idx3])
    e < best3[1] && (global best3 = (e, q, nn))
end
_, q, n_night = best3
n_hat = ifelse.(sun, n_day, n_night)
c = fit_emissivity(d[:lwd][cal], B[cal], ΔT[cal], n_hat[cal], q)
ε_new = clamp.(emissivity(c, ΔT, n_hat, q), 0.4, 1.0)

# previous linear scheme (Chion before the recalibration), uncapped and capped at one
n_old = ifelse.(sun, clamp.(1 .- τ ./ (0.85 .+ 0.075 .* z_km), 0, 1), 0.389)
ε_old = clamp.(0.624 .+ 0.0032 .* ΔT .+ 0.613 .* n_old, 0.4, 1.3)
ε_old_cap = min.(ε_old, 1.0)

params = [("lw_clear_sky_transmissivity", t[1]), ("lw_clear_sky_transmissivity_per_km", t[2]),
    ("lw_clear_sky_transmissivity_insolation", t[3]), ("lw_cloud_transmission_coefficient", a),
    ("lw_cloud_transmission_exponent", p), ("lw_cloud_transmission_decay_per_km", b),
    ("lw_clear_sky_emissivity", c[1]), ("lw_clear_sky_emissivity_temperature_slope", c[2]),
    ("lw_overcast_emissivity", c[3]), ("lw_overcast_emissivity_temperature_slope", c[4]),
    ("lw_cloud_emissivity_exponent", q), ("lw_night_cloud_fraction", n_night)]
open(joinpath(OUT, "fitted_parameters.csv"), "w") do io
    println(io, "parameter,value")
    for (k, v) in params; println(io, "$k,$v"); end
end
for (k, v) in params; @printf("%-44s %.5g\n", k, v); end
# compare with the model defaults
cm = Chion.SnowpackPhysicalConstants()
for (k, v) in params
    abs(getproperty(cm, Symbol(k)) - v) > 0.01 * max(abs(v), 1e-3) &&
        @warn "Model default differs from the fit by more than 1 %" k default=getproperty(cm, Symbol(k)) fit=v
end

# ---------- statistics ----------
function scheme_stats(name, ε)
    for (per, P) in (("1980–1989", cal), ("2010–2019", val)), (sub, S) in (("all days", trues(length(ε))), ("JJA", jja), ("JJA melt days", jja .& melt_day))
        k = P .& S; err = ε[k] .* B[k] .- d[:lwd][k]
        stat!("L1", "$name, $per, $sub: LWD bias", mean(err), "W m⁻²")
        stat!("L1", "$name, $per, $sub: LWD RMSE", rmse(err), "W m⁻²")
    end
    stat!("L1", "$name: share of days with ε > 1", 100mean(ε .> 1), "%")
end
scheme_stats("previous", ε_old); scheme_stats("previous capped", ε_old_cap); scheme_stats("physical", ε_new)
for (lab, nn) in (("physical", n_hat), ("previous", n_old))
    k = val .& sun
    stat!("L1", "$lab cloud fraction vs MAR (2010–2019, sunlit): bias", mean(nn[k] .- d[:cloud][k]), "")
    stat!("L1", "$lab cloud fraction vs MAR (2010–2019, sunlit): r", cor(nn[k], d[:cloud][k]), "")
end
stat!("L1", "clear-sky transmissivity RMSE (2010–2019)", rmse(τcs[val .& sun .& (τcs_mar .> 0.3)] .- τcs_mar[val .& sun .& (τcs_mar .> 0.3)]), "")

# ---------- figure ----------
binned(x, y, edges, sel; f=median) = [(k = sel .& (x .>= lo) .& (x .< hi); count(k) > 200 ? f(y[k]) : NaN) for (lo, hi) in zip(edges[1:end-1], edges[2:end])]
mids(e) = (e[1:end-1] .+ e[2:end]) ./ 2
q10(v) = quantile(v, 0.1); q90(v) = quantile(v, 0.9); q25(v) = quantile(v, 0.25); q75(v) = quantile(v, 0.75)

fmt(x, d=3) = replace(@sprintf("%.*f", d, x), "-" => "−")
MARK, LINE = "MAR, median per bin", "fit"
legend_keys!(p) = (scatter!(p, [NaN], [NaN]; c=INK, ms=3, msw=0, label=MARK); plot!(p, [NaN], [NaN]; c=INK, label=LINE))

# (a) clear-sky transmissivity: MAR clear-sky SW↓ / TOA SW↓ (median and 10–90 % range per elevation bin)
ze = 0:0.25:3.5; zm = mids(ze); s1 = cal .& sun .& (τcs_mar .> 0.3)
med = binned(z_km, τcs_mar, ze, s1)
pa = plot(zm, med; ribbon=(med .- binned(z_km, τcs_mar, ze, s1; f=q10), binned(z_km, τcs_mar, ze, s1; f=q90) .- med), c=INK, fillalpha=0.15,
    label="MAR, median and 10–90 %", xlabel="Surface elevation z (km)", ylabel="Clear-sky transmissivity τcs", legend=:topleft,
    ylims=(0.5, 1.05), xlims=(0, 3.5))
for (toa, ls) in ((150, :dot), (300, :solid), (450, :dash))
    plot!(pa, zm, clear_sky(t, zm, fill(toa, length(zm))); c=C1, ls, lw=toa == 300 ? 1.5 : 1.0, label="fit, S = $toa W m⁻²")
end
tag!(pa, 1)
note!(pa, "τcs = τ₀ + τ₁ z + τ₂ ln(S / 300 W m⁻²)\nτ₀ = $(fmt(t[1])), τ₁ = $(fmt(t[2])) km⁻¹, τ₂ = $(fmt(t[3]))\nS: daily mean top-of-atmosphere SW↓"; loc=:bottomright, size=6)

# (b) cloud transmission: τ/τcs against MAR cloud cover, median per 0.1 bin, by elevation band
ce = 0:0.1:1.0; cmid = mids(ce)
pb = plot(; xlabel="MAR cloud cover n", ylabel="Cloud transmission τ/τcs", legend=:bottomleft, ylims=(0.4, 1.17), xlims=(0, 1))
for ((lo, hi), col) in zip(((0, 1), (1, 2), (2, 4)), (C2, C4, C1))
    k = cal .& sun .& (d[:toa] .> 100) .& (z_km .>= lo) .& (z_km .< hi)
    scatter!(pb, cmid, binned(d[:cloud], τ ./ τcs, ce, k); c=col, ms=3, msw=0, label=@sprintf("z = %.0f–%.0f km", lo, hi))
    xs = range(0, 1, length=50)
    plot!(pb, xs, 1 .- a .* xs .^ p .* exp(-b * (lo + hi) / 2); c=col, label="")
end
legend_keys!(pb)
tag!(pb, 2)
note!(pb, "τ/τcs = 1 − a nᵖ exp(−b z)\na = $(fmt(a)), p = $(fmt(p, 2)), b = $(fmt(b)) km⁻¹"; loc=:topright, size=6)

# (c) diagnosed cloud fraction (inversion of b) against MAR cloud cover, validation years, sunlit days
kc = val .& sun; m = binned(d[:cloud], n_hat, ce, kc)
pc = plot([0, 1], [0, 1]; c=MUTED, lw=0.75, label="", xlabel="MAR cloud cover", ylabel="Diagnosed cloud fraction n̂", legend=:topleft,
    xlims=(0, 1), ylims=(0, 1.05))
plot!(pc, cmid, m; ribbon=(m .- binned(d[:cloud], n_hat, ce, kc; f=q25), binned(d[:cloud], n_hat, ce, kc; f=q75) .- m), c=C1, fillalpha=0.15,
    label="median and 25–75 % per bin")
tag!(pc, 3)
note!(pc, "2010–2019, sunlit days\nr = $(fmt(cor(n_hat[kc], d[:cloud][kc]), 2)), bias $(fmt(mean(n_hat[kc] .- d[:cloud][kc]), 2))\npolar night: n̂ = $(fmt(n_night, 2))"; loc=:bottomright, size=6)

# (d) effective emissivity MAR LW↓/σT⁴ against diagnosed cloud fraction, median per 0.1 bin, by temperature band
ne = 0:0.1:1.0; nm = mids(ne)
pd = plot(; xlabel="Diagnosed cloud fraction n̂", ylabel="Effective emissivity LW↓/σT⁴", legend=:topleft, ylims=(0.5, 1.12), xlims=(0, 1))
for ((lo, hi), col) in zip(((-35, -20), (-20, -10), (-10, -3), (-3, 3)), (C1, C3, C4, C2))
    k = cal .& sun .& (ΔT .>= lo) .& (ΔT .< hi)
    scatter!(pd, nm, binned(n_hat, d[:lwd] ./ B, ne, k); c=col, ms=3, msw=0, label=replace("T = $lo to $hi °C", "-" => "−"))
    xs = range(0, 1, length=50)
    plot!(pd, xs, emissivity(c, (lo + hi) / 2, xs, q); c=col, label="")
end
legend_keys!(pd)
hline!(pd, [1]; c=MUTED, lw=0.75, label="")
tag!(pd, 4)
sup = Dict('0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴', '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹', '.' => '·')
nq = "n̂" * map(ch -> sup[ch], q == round(q) ? string(Int(q)) : string(q))
note!(pd, "ε = εcs (1 − $nq) + εoc $nq ≤ 1\nεcs = $(fmt(c[1])) + $(fmt(c[2], 4)) ΔT\nεoc = $(fmt(c[3])) + $(fmt(c[4], 4)) ΔT  (ΔT in °C)"; loc=:bottomright, size=6)

# (e) daily LW↓, validation years
edges = 100:5:350
H = zeros(length(edges) - 1, length(edges) - 1)
model_lwd = ε_new .* B
for i in findall(val)
    ix = searchsortedlast(edges, d[:lwd][i]); iy = searchsortedlast(edges, model_lwd[i])
    (1 <= ix < length(edges)) && (1 <= iy < length(edges)) && (H[iy, ix] += 1)
end
hmax = ceil(maximum(log10.(H[H .> 0])); digits=1)
pe = heatmap(mids(edges), mids(edges), log10.(replace(H, 0 => NaN)); c=cgrad(:batlow), clims=(0, hmax), colorbar=false,
    xlabel="MAR LW↓ (W m⁻²)", ylabel="Parameterised LW↓ (W m⁻²)", xlims=(100, 350), ylims=(100, 350), right_margin=-1Plots.mm)
# narrow colour bar as its own subplot next to (e)
levels = range(0, hmax, length=200)
pcb = heatmap([0, 1], levels, repeat(collect(levels), 1, 2); c=cgrad(:batlow), clims=(0, hmax), colorbar=false, xticks=false,
    ymirror=true, ylabel="log₁₀ number of cell-days", framestyle=:box, tickfontsize=6, guidefontsize=7, ylims=(0, hmax),
    left_margin=-2Plots.mm, right_margin=3Plots.mm)
plot!(pe, [100, 350], [100, 350]; c=INK, lw=0.75, label="")
err = model_lwd[val] .- d[:lwd][val]
tag!(pe, 5)
note!(pe, "2010–2019, daily\nbias $(fmt(mean(err), 1))\nRMSE $(fmt(rmse(err), 1))"; loc=:topleft, size=6)

# (f) LW↓ bias by elevation, validation years
be = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5]; bm = mids(be)
pf = plot(; xlabel="Surface elevation z (km)", ylabel="LW↓ bias, parameterised − MAR (W m⁻²)", legend=:topright, xlims=(0, 3.5))
hline!(pf, [0]; c=MUTED, lw=0.75, label="")
for (sel, col, lab) in ((trues(length(val)), INK, "all days"), (jja, C1, "JJA"), (jja .& melt_day, C2, "JJA melt days"))
    k = val .& sel
    plot!(pf, bm, binned(z_km, model_lwd .- d[:lwd], be, k; f=mean); c=col, marker=:circle, ms=3, msw=0, label=lab)
end
with_headroom!(pf, 0.25)
tag!(pf, 6)
note!(pf, "2010–2019, mean per bin"; loc=:topleft, size=6)

layout = @layout [a b c; d [e{0.955w} cb] f]
fig = plot(pa, pb, pc, pd, pe, pcb, pf; layout, size=(960, 600), left_margin=4Plots.mm, bottom_margin=4Plots.mm, right_margin=2Plots.mm)
save_figure(fig, OUT, "L1_cloud_proxy_longwave")
write_stats(OUT)
println("Figure and statistics written to ", OUT)
