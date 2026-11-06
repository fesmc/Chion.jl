#!/usr/bin/env julia

# Verification suite for the model description paper. Idealised single-column
# tests against analytical solutions (heat conduction, phase change), exact
# conservation checks on a full Greenland year, and discretisation convergence.
# Writes paper-ready figures (PDF + PNG), metrics.csv and README.md to
# output/verification (or CHION_VERIFICATION_DIR).
#
#   julia --project=. --threads 8 examples/figures/verification_suite.jl
include(joinpath(@__DIR__, "..", "research", "run_gris_fully_prescribed_transient.jl"))
using Printf, Random

const VOUT = get(ENV, "CHION_VERIFICATION_DIR", joinpath(@__DIR__, "..", "..", "output", "verification"))
const SECONDS_PER_DAY = 86_400.0
const RHO_SNOW = 350.0

# ---------- figure style ----------
# Validated reference palette (categorical slots in fixed order; one-hue ramp for
# ordered sequences such as resolution or time).
const C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
const RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]
const INK = "#3a3a38"
Plots.default(fontfamily="DejaVu Sans", titlefontsize=9, guidefontsize=8, tickfontsize=7,
    legendfontsize=7, framestyle=:box, grid=false, linewidth=1.5, foreground_color_legend=nothing,
    background_color_legend=nothing, titlelocation=:left, dpi=300)
panel_title(i) = "($(Char('a' + i - 1)))"
"""Panel description inside the axes at relative position `at` as seen on the page, (0, 0)
bottom left, (1, 1) top right (also for flipped y axes). Call after the last series so the
axis limits are final."""
function panel_text!(p, s; at=(0.04, 0.95), halign=:left, valign=:top)
    y = p[1][:yaxis][:flip] ? 1 - at[2] : at[2]
    annotate!(p, ((at[1], y), text(s, 8, halign, valign, INK)))
end
function save_figure(fig, name)
    mkpath(VOUT)
    savefig(fig, joinpath(VOUT, "$name.pdf")); savefig(fig, joinpath(VOUT, "$name.png"))
    return name
end

const METRICS = NamedTuple[]
record!(test, metric, value, unit, expectation) = push!(METRICS, (; test, metric, value, unit, expectation))

# ---------- analytical helpers ----------
"""erf(x) by its Taylor series (|x| < 3, accurate to ~1e-12)."""
function erf_series(x)
    term, total, n = x, x, 0
    while abs(term) > 1e-17 * abs(total)
        n += 1
        term *= -x^2 / n
        total += term / (2n + 1)
    end
    return 2 / sqrt(pi) * total
end

"""Scaled complementary error function exp(x²)·erfc(x) for x ≥ 0."""
function erfcx(x)
    x < 3 && return exp(x^2) * (1 - erf_series(x))
    # Continued fraction (modified Lentz) for large x.
    tiny = 1e-300
    f = x; C = x; D = 0.0
    for n in 1:200
        a = n / 2
        D = x + a * D; D = abs(D) < tiny ? tiny : D; D = 1 / D
        C = x + a / C; C = abs(C) < tiny ? tiny : C
        delta = C * D; f *= delta
        abs(delta - 1) < 1e-16 && break
    end
    return 1 / (sqrt(pi) * f)
end
erfc_(x) = x < 0 ? 2 - erfc_(-x) : exp(-x^2) * erfcx(x)

"""Step change of the ambient temperature over a semi-infinite solid with
surface heat-transfer coefficient h (h → ∞: fixed surface temperature);
Carslaw & Jaeger (1959, ch. 2), Incropera et al. (2007, sect. 5.7)."""
function step_response(z, t, κ, K, h)
    η = z / (2sqrt(κ * t))
    isinf(h) && return erfc_(η)
    return erfc_(η) - exp(-η^2) * erfcx(η + h * sqrt(κ * t) / K)
end

"""Complex amplitude of a periodic surface temperature wave in snow (thickness H,
conductivity K1, diffusivity κ1) over semi-infinite ice (K2, κ2); unit surface amplitude."""
function two_layer_wave(z, ω, H, K1, κ1, K2, κ2)
    k1 = (1 + im) * sqrt(ω / (2κ1)); k2 = (1 + im) * sqrt(ω / (2κ2))
    # θ1 = a e^{-k1 z} + b e^{k1 z}, θ2 = c e^{-k2 (z - H)}
    M = [1 1 0;
         exp(-k1 * H) exp(k1 * H) -1;
         -K1 * k1 * exp(-k1 * H) K1 * k1 * exp(k1 * H) K2 * k2]
    a, b, c = M \ ComplexF64[1, 0, 0]
    return z <= H ? a * exp(-k1 * z) + b * exp(k1 * z) : c * exp(-k2 * (z - H))
end

"""Amplitude and phase lag (rad) of a sampled periodic series over one full period."""
function harmonic(values, phases)
    s = 2 * sum(values .* sin.(phases)) / length(values)
    c = 2 * sum(values .* cos.(phases)) / length(values)
    return hypot(s, c), mod(-atan(c, s), 2pi)
end

# ---------- single-column harness ----------
"""Column model for idealised energy tests: BESSI turbulent sensible heat with
coefficient `h` (a Robin condition towards the air temperature), no emission,
optional ice substrate."""
column_model(nlayer; h=1e6, nice=0, ice_top=0.05) = BESSIModel(SnowpackGrid(1); Ntot=max(nlayer, 1),
    turbulent_flux_scheme=:bessi, D_sh=h, ϵ_snow=0.0, eps_ice=0.0, ice_substrate_layers=nice,
    ice_substrate_top_thickness_m=ice_top, diurnal_shortwave_substeps=false, longwave_scheme=:graybody)

function column_state(model, dz, T; rho=RHO_SNOW)
    state = BESSIState(model)
    n = length(dz)
    state.N[1] = n
    for k in 1:n
        state.mass[k, 1] = rho * dz[k]; state.density[k, 1] = rho
        state.temperature[k, 1] = T; state.mass_w[k, 1] = 0.0
    end
    state.ice_temperature .= T
    state.Tsrf[1] = T
    return state
end

"""One energy solve with only the ambient-air exchange (q_sh not prescribed)."""
solve!(state, scratch, Ta, dt) = Chion.go_energy_flux!(state, 1, Ta, 0.0, 0.0, 0.0, dt;
    scratch, q_sw_net=0.0, q_lw_down=0.0, q_lh=0.0)

"""One energy solve with a prescribed net surface flux and no ambient exchange."""
solve_flux!(state, scratch, Q, dt) = Chion.go_energy_flux!(state, 1, 273.15, 0.0, 0.0, 0.0, dt;
    scratch, q_sw_net=Q, q_lw_down=0.0, q_sh=0.0, q_lh=0.0)

conductivity(rho, T) = Chion._snow_thermal_conductivity(rho, T, 917.0)
centres(dz) = cumsum(dz) .- dz ./ 2
substrate_dz(nice, top) = [top * 2.0^(k - 1) for k in 1:nice]

"""Heat content (J m⁻², relative to T0) of the snow layers and the ice substrate."""
function heat_content(state, c, ice_dz)
    H = sum(state.mass[k, 1] * c.ci * (state.temperature[k, 1] - c.T0) for k in 1:state.N[1]; init=0.0)
    H += sum(c.rho_i * ice_dz[k] * c.ci * (state.ice_temperature[k, 1] - c.T0) for k in eachindex(ice_dz); init=0.0)
    return H
end

# ---------- 1. heat conduction ----------
function run_step_response()
    ci = 2110.0; Ti, Ta = 263.15, 264.15
    dz = fill(0.02, 300); z = centres(dz)
    K = conductivity(RHO_SNOW, (Ti + Ta) / 2); κ = K / (RHO_SNOW * ci)
    days = (1, 10, 50); dt = 3600.0
    panels = Any[]
    for (i, (label, h)) in enumerate((("Fixed surface\ntemperature", Inf), ("Surface heat\nexchange,\nh = 5 W m⁻² K⁻¹", 5.0)))
        model = column_model(length(dz); h=isinf(h) ? 1e6 : h)
        state = column_state(model, dz, Ti); scratch = Chion.EnergyWorkspace(state)
        p = plot(; title=panel_title(i), xlabel="Normalised temperature (T − Tᵢ)/(Tₐ − Tᵢ)",
            ylabel=i == 1 ? "Depth (m)" : "", yflip=true, ylims=(0, 2.5), xlims=(-0.02, 1.02), legend=i == 1 ? :bottomright : false)
        step = 0
        for (j, day) in enumerate(days)
            while step * dt < day * SECONDS_PER_DAY - 1e-6
                solve!(state, scratch, Ta, dt); step += 1
            end
            model_theta = (state.temperature[1:length(dz), 1] .- Ti) ./ (Ta - Ti)
            exact = step_response.(z, day * SECONDS_PER_DAY, κ, K, h)
            colour = RAMP[2j - 1]
            plot!(p, exact, z; c=colour, lw=1.5, label="$day d, analytical")
            sel = 1:3:length(z)
            scatter!(p, model_theta[sel], z[sel]; c=colour, ms=2.5, msw=0, label="$day d, Chion")
            record!("step response, $(isinf(h) ? "fixed surface" : "h = 5")", "max |θ error| at $day d", maximum(abs.(model_theta .- exact)), "1", "≪ 1 (dt = 1 h, Δz = 2 cm)")
        end
        panel_text!(p, label; at=(0.96, 0.55), halign=:right)
        push!(panels, p)
    end
    save_figure(plot(panels...; layout=(1, 2), size=(720, 330), left_margin=3Plots.mm, bottom_margin=3Plots.mm), "fig1_step_response")
end

"""Periodic surface forcing: returns depth, model amplitude and lag for snow and substrate rows."""
function periodic_column(dz_snow, rho_snow, nice, ice_top, period_days, dt, ncycles; Tm=258.15, amplitude=1.0)
    model = column_model(length(dz_snow); nice, ice_top)
    state = column_state(model, dz_snow, Tm; rho=rho_snow)
    isempty(dz_snow) && (state.N[1] = 0)
    scratch = Chion.EnergyWorkspace(state)
    ω = 2pi / (period_days * SECONDS_PER_DAY)
    nsteps = round(Int, period_days * SECONDS_PER_DAY / dt)
    ice_dz = substrate_dz(nice, ice_top)
    series = zeros(length(dz_snow) + nice, nsteps); phases = zeros(nsteps)
    for cycle in 1:ncycles, s in 1:nsteps
        t = ((cycle - 1) * nsteps + s) * dt
        solve!(state, scratch, Tm + amplitude * sin(ω * t), dt)
        if cycle == ncycles
            series[:, s] .= vcat(state.temperature[1:length(dz_snow), 1], state.ice_temperature[:, 1]) .- Tm
            phases[s] = ω * t
        end
    end
    widths = vcat(dz_snow, ice_dz)
    bottoms = cumsum(widths); tops = bottoms .- widths
    fits = [harmonic(series[k, :], phases) for k in axes(series, 1)]
    return (tops .+ bottoms) ./ 2, tops, bottoms, first.(fits) ./ amplitude, last.(fits), ω
end

"""Mean of `f` over the cell [a, b] (finite-volume cell average), trapezoidal rule."""
cell_mean(f, a, b; n=200) = (zs = range(a, b; length=n + 1); v = f.(zs); (sum(v) - (v[1] + v[end]) / 2) / n)

function run_periodic()
    ci = 2110.0; Tm = 258.15
    K_snow = conductivity(RHO_SNOW, Tm); κ_snow = K_snow / (RHO_SNOW * ci)
    K_ice = conductivity(917.0, Tm); κ_ice = K_ice / (917.0 * ci)
    cases = [
        ("Snow, annual", fill(0.05, 240), 0, 0.05, 365.0, SECONDS_PER_DAY, 6),
        ("1 m snow over ice, annual", fill(0.05, 20), 8, 0.05, 365.0, SECONDS_PER_DAY, 6),
        ("Bare ice, daily", Float64[], 10, 0.005, 1.0, 600.0, 12),
    ]
    amp_panels = Any[]; lag_panels = Any[]
    for (i, (label, dz, nice, top, period, dt, ncycles)) in enumerate(cases)
        z, tops, bottoms, amp, lag, ω = periodic_column(dz, RHO_SNOW, nice, top, period, dt, ncycles; Tm)
        H = sum(dz; init=0.0)
        wave(zz) = isempty(dz) ? exp(-(1 + im) * sqrt(ω / (2κ_ice)) * zz) :
                   nice == 0 ? exp(-(1 + im) * sqrt(ω / (2κ_snow)) * zz) :
                   two_layer_wave(zz, ω, H, K_snow, κ_snow, K_ice, κ_ice)
        # The model holds finite-volume cell means, so compare with analytical cell means.
        exact = [cell_mean(wave, a, b) for (a, b) in zip(tops, bottoms)]
        damping_depth = isempty(dz) ? sqrt(2κ_ice / ω) : nice == 0 ? sqrt(2κ_snow / ω) : sqrt(2κ_ice / ω)
        # Stay well above the insulated base, whose reflection the analytical
        # semi-infinite solution does not contain.
        zmax = min(H + 4damping_depth, 0.6 * maximum(bottoms))
        zfine = range(0, zmax; length=400); exact_fine = wave.(zfine)
        keep = bottoms .<= zmax; keep_fine = trues(length(zfine))
        a = plot(; title=panel_title(i), xlabel="Relative amplitude", ylabel=i == 1 ? "Depth (m)" : "",
            yflip=true, xscale=:log10, legend=i == 1 ? :bottomright : false, ylims=(0, zmax))
        plot!(a, abs.(exact_fine[keep_fine]), zfine[keep_fine]; c=C1, label="Analytical")
        scatter!(a, abs.(exact[keep]), z[keep]; c=:white, msc=C1, msw=1, ms=3.5, label="Analytical cell mean")
        scatter!(a, amp[keep], z[keep]; c=C2, ms=2.5, msw=0, label="Chion")
        nice > 0 && !isempty(dz) && hline!(a, [H]; c=INK, ls=:dot, lw=0.8, label="")
        l = plot(; title=panel_title(i + 3), xlabel="Phase lag (rad)", ylabel=i == 1 ? "Depth (m)" : "",
            yflip=true, legend=false, ylims=(0, zmax))
        plot!(l, unwrap(-angle.(exact_fine[keep_fine])), zfine[keep_fine]; c=C1)
        scatter!(l, unwrap(-angle.(exact[keep])), z[keep]; c=:white, msc=C1, msw=1, ms=3.5)
        scatter!(l, unwrap(lag[keep]), z[keep]; c=C2, ms=2.5, msw=0)
        nice > 0 && !isempty(dz) && hline!(l, [H]; c=INK, ls=:dot, lw=0.8)
        i == 1 ? panel_text!(a, label) : panel_text!(a, label; at=(0.96, 0.05), halign=:right, valign=:bottom)
        panel_text!(l, label; at=(0.04, 0.05), valign=:bottom)
        push!(amp_panels, a); push!(lag_panels, l)
        relevant = keep .& (abs.(exact) .> 0.05)
        record!("periodic: $label", "max relative amplitude error vs cell mean (|A| > 0.05)",
            maximum(abs.(amp[relevant] .- abs.(exact[relevant])) ./ abs.(exact[relevant])), "1", "≪ 1")
        record!("periodic: $label", "max phase-lag error vs cell mean (|A| > 0.05)",
            maximum(abs.(unwrap(lag[relevant]) .- unwrap(-angle.(exact[relevant])))), "rad", "≪ 1")
    end
    save_figure(plot(amp_panels..., lag_panels...; layout=(2, 3), size=(780, 560), left_margin=3Plots.mm, bottom_margin=3Plots.mm),
        "fig2_periodic_forcing")
end

"""Remove 2π jumps from a phase profile."""
function unwrap(phase)
    out = copy(phase)
    for k in 2:length(out)
        while out[k] - out[k - 1] > pi; out[k] -= 2pi; end
        while out[k] - out[k - 1] < -pi; out[k] += 2pi; end
    end
    return out
end

function run_convergence()
    ci = 2110.0; Ti = 263.15; ΔT = 0.01; depth = 4.0; t_end = 10 * SECONDS_PER_DAY
    K = conductivity(RHO_SNOW, Ti + ΔT / 2); κ = K / (RHO_SNOW * ci)
    function integrate(dz, dt)
        n = round(Int, depth / dz)
        model = column_model(n); state = column_state(model, fill(dz, n), Ti); scratch = Chion.EnergyWorkspace(state)
        for _ in 1:round(Int, t_end / dt); solve!(state, scratch, Ti + ΔT, dt); end
        return centres(fill(dz, n)), (state.temperature[1:n, 1] .- Ti) ./ ΔT
    end
    rms(x) = sqrt(sum(abs2, x) / length(x))
    dzs = [0.2, 0.1, 0.05, 0.025, 0.0125]
    space_errors = map(dzs) do dz
        z, θ = integrate(dz, 10.0); sel = z .<= 2.0
        rms(θ[sel] .- step_response.(z[sel], t_end, κ, K, Inf))
    end
    dts = [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625] .* SECONDS_PER_DAY
    z_ref, θ_ref = integrate(0.0125, 60.0)
    time_errors = map(dts) do dt
        z, θ = integrate(0.0125, dt); sel = z .<= 2.0
        rms(θ[sel] .- θ_ref[sel])
    end
    # Least-squares slope of log(error) against log(resolution).
    order(err, h) = (x = log.(h); y = log.(err); sum((x .- mean(x)) .* (y .- mean(y))) / sum((x .- mean(x)) .^ 2))
    p_space, p_time = order(space_errors, dzs), order(time_errors, dts)
    for (h, err) in zip(dzs, space_errors); record!("convergence", @sprintf("RMS error, Δz = %.4f m (Δt = 10 s)", h), err, "1", ""); end
    for (h, err) in zip(dts, time_errors); record!("convergence", @sprintf("RMS difference, Δt = %.4f d (Δz = 1.25 cm)", h / SECONDS_PER_DAY), err, "1", ""); end
    record!("convergence", "observed spatial order (least-squares fit)", p_space, "1", "2 (second-order finite volume)")
    record!("convergence", "observed temporal order (least-squares fit)", p_time, "1", "1 (implicit Euler)")
    guide(h, err, p) = exp(mean(log.(err) .- p .* log.(h))) .* h .^ p
    a = plot(dzs, space_errors; xscale=:log10, yscale=:log10, marker=:circle, ms=4, msw=0, c=C1,
        label=@sprintf("Chion (order %.2f)", p_space), title=panel_title(1),
        xlabel="Layer thickness Δz (m)", ylabel="RMS error of normalised temperature", legend=:topleft)
    plot!(a, dzs, guide(dzs, space_errors, 2); c=INK, ls=:dash, lw=1, label="Second order")
    b = plot(dts ./ SECONDS_PER_DAY, time_errors; xscale=:log10, yscale=:log10, marker=:circle, ms=4, msw=0, c=C1,
        label=@sprintf("Chion (order %.2f)", p_time), title=panel_title(2),
        xlabel="Time step Δt (days)", ylabel="RMS difference to Δt = 60 s", legend=:topleft)
    plot!(b, dts ./ SECONDS_PER_DAY, guide(dts, time_errors, 1); c=INK, ls=:dash, lw=1, label="First order")
    panel_text!(a, "Spatial convergence"; at=(0.96, 0.05), halign=:right, valign=:bottom)
    panel_text!(b, "Temporal convergence"; at=(0.96, 0.05), halign=:right, valign=:bottom)
    save_figure(plot(a, b; layout=(1, 2), size=(720, 320), left_margin=4Plots.mm, bottom_margin=4Plots.mm), "fig3_convergence")
end

# ---------- 2. conservation ----------
"""Exact column budgets at every step of a full Greenland year (all processes, all columns,
default parameterized setup).
Mass: Δ(snow solid + liquid) + Δsmb_ice = precipitation − runoff + vapour exchange.
Energy: Δ`column_enthalpy` = Δ`boundary_enthalpy`, i.e. the change of the enthalpy of the
snow/firn column and ice substrate (relative to liquid water at T₀) equals the surface energy
applied by the solver plus the enthalpy of the mass that crossed the column boundary."""
function run_full_step_conservation()
    p = joinpath(DOMAIN_SETUP.forcing_dir, "MARv3.14.3-10km-daily-ERA5-1990.nc")
    l = load_mar(p); g = l.grid
    forcing = parameterized_forcing(l.forcing)
    model = BESSIModel(g)
    spin = Simulation(model; forcing, years=5, backend=:threads, write_netcdf=false, compute_year_metrics=false)
    run!(spin; io=devnull)
    sim = Simulation(model; forcing, state=spin.now, years=1, backend=:threads, write_netcdf=false, compute_year_metrics=false)
    integ = init_integrator(sim; io=devnull)
    s = integ.model_runtime.data.state
    w = integ.model_runtime.data.workspace
    storage() = [sum(s.mass[k, j] + s.mass_w[k, j] for k in 1:s.N[j]; init=0.0) for j in 1:s.ncol]
    enthalpy() = [Chion.column_enthalpy(s, j) for j in 1:s.ncol]
    snap() = (S=storage(), smb=copy(s.smb_ice), R=copy(s.runoff), V=copy(s.vapor_mass), H=enthalpy(), B=copy(w.boundary_enthalpy))
    prev = snap(); mass_res = Float64[]; energy_res = Float64[]; precip_total = 0.0; energy_scale = 0.0
    snowfall, rainfall = Array(forcing.snowfall_rate), Array(forcing.rainfall_rate)
    for k in eachindex(forcing.dt_days)
        step!(integ)
        now = snap()
        dt = forcing.dt_days[k] * SECONDS_PER_DAY
        P = (snowfall[:, k] .+ rainfall[:, k]) .* dt
        append!(mass_res, (now.S .- prev.S) .+ (now.smb .- prev.smb) .- (P .- (now.R .- prev.R) .+ (now.V .- prev.V)))
        append!(energy_res, ((now.H .- prev.H) .- (now.B .- prev.B)) ./ dt)
        precip_total += sum(P); energy_scale += sum(abs, now.B .- prev.B) / dt
        prev = now
    end
    finalize!(integ)
    n = length(mass_res)
    record!("mass conservation (GrIS 1990, all columns and days)", "max |residual|", maximum(abs, mass_res), "kg m⁻²", "round-off")
    record!("mass conservation (GrIS 1990, all columns and days)", "mean daily precipitation", precip_total / n, "kg m⁻²", "scale")
    record!("energy conservation, full step (GrIS 1990, all columns and days)", "max |residual|", maximum(abs, energy_res), "W m⁻²", "round-off")
    record!("energy conservation, full step (GrIS 1990, all columns and days)", "mean |boundary enthalpy flux|", energy_scale / n, "W m⁻²", "scale")
    return mass_res, energy_res, length(g.js)
end

"""Energy closure of the implicit energy solve alone on real columns (prescribed MAR fluxes):
surface energy input = change in heat content of snow and ice substrate + melt energy.
All columns, threaded over columns."""
function run_energy_closure()
    p = joinpath(DOMAIN_SETUP.forcing_dir, "MARv3.14.3-10km-daily-ERA5-1990.nc")
    l = load_mar(p); g, f = l.grid, l.forcing
    pf = fully_prescribed_forcing(p, l)
    model = BESSIModel(g)
    spin = Simulation(model; forcing=pf, years=5, backend=:threads, write_netcdf=false, compute_year_metrics=false)
    run!(spin; io=devnull)
    c = model.c; ice_dz = substrate_dz(model.ice_substrate_layers, model.ice_substrate_top_thickness_m)
    swn, lwd, shf, lhf, ta = Array(pf.q_sw_net), Array(pf.q_lw_down), Array(pf.q_sh), Array(pf.q_lh), Array(pf.air_temperature)
    ndays, ncol = length(f.dt_days), length(g.js)
    residuals = zeros(ndays, ncol)
    Threads.@threads for j in 1:ncol
        single = BESSIModel(SnowpackGrid(1))
        state = BESSIState(single)
        for name in Chion._BESSI_ARRAY_FIELD_NAMES
            src, dst = getfield(spin.now, name), getfield(state, name)
            dst isa AbstractMatrix ? (dst[:, 1] .= src[:, j]) : (dst[1] = src[j])
        end
        scratch = Chion.EnergyWorkspace(state)
        dt = SECONDS_PER_DAY
        for k in 1:ndays
            H0 = heat_content(state, c, ice_dz)
            e = Chion.go_energy_flux!(state, 1, ta[j, k], 0.0, 0.0, 0.0, dt; scratch,
                q_sw_net=swn[j, k], q_lw_down=lwd[j, k], q_sh=shf[j, k], q_lh=lhf[j, k])
            E_in = (e.surface_flux_constant - e.surface_flux_linear * state.Tsrf[1]) * dt
            residuals[k, j] = (E_in - (heat_content(state, c, ice_dz) - H0) - e.melt_energy_available) / dt
        end
    end
    residuals = vec(residuals)
    record!("energy closure of the energy solve ($ncol GrIS columns × $ndays d)", "max |residual|", maximum(abs, residuals), "W m⁻²", "round-off")
    return residuals
end

function plot_conservation(mass_res, solve_energy_res)
    hist_panel(i, title, values, unit) = begin
        absval = abs.(values)
        top = max(-6, ceil(Int, log10(max(maximum(absval), 1e-300))))   # widen beyond 1e-6 if needed
        nonzero = clamp.(absval[absval .> 0], 1e-18, 10.0^top)
        zeros_frac = 1 - count(>(0), absval) / length(absval)
        edges = 10.0 .^ range(-18, top; length=4 * (top + 18) + 1)
        counts = [count(x -> edges[b] <= x < edges[b + 1], nonzero) for b in 1:length(edges) - 1]
        ymax = max(maximum(counts; init=0), 1)
        share = 100 .* counts ./ length(absval); smax = 100ymax / length(absval)
        h = plot(edges, vcat(share, share[end]); seriestype=:steppost, fillrange=0, c=C1, fillcolor=C1, fillalpha=0.85,
            xscale=:log10, label="", title=panel_title(i), xlabel="|residual| ($unit)", ylabel="Share of samples (%)",
            xlims=(1e-18, 10.0^top), ylims=(0, 1.6smax))
        panel_text!(h, title * @sprintf("\nmax %.1e %s\nN = %d, %.0f %% exactly zero",
            maximum(absval), unit, length(absval), 100zeros_frac))
        h
    end
    a = hist_panel(1, "Column mass balance, every step", mass_res, "kg m⁻²")
    b = hist_panel(2, "Energy-solve closure", solve_energy_res, "W m⁻²")
    save_figure(plot(a, b; layout=(1, 2), size=(720, 300), left_margin=4Plots.mm, bottom_margin=4Plots.mm), "fig4_conservation")
end

# ---------- 3. phase change ----------
function run_phase_change()
    c = BESSIModel(SnowpackGrid(1)).c
    panels = Any[]

    # (a, b) Constant net flux on cold snow: surface warming until melt onset
    # (semi-infinite constant-flux solution) and melt afterwards.
    F, Ti, dt = 5.0, c.T0 - 2.0, 30.0
    dz = fill(0.005, 200)
    model = column_model(length(dz))
    state = column_state(model, dz, Ti); scratch = Chion.EnergyWorkspace(state)
    K = conductivity(RHO_SNOW, Ti + 1.0); κ = K / (RHO_SNOW * c.ci)
    t_onset_exact = pi * (K * (c.T0 - Ti) / (2F))^2 / κ
    nsteps = round(Int, 4 * t_onset_exact / dt)
    times = (1:nsteps) .* dt; Ts = zeros(nsteps); melt_rate = zeros(nsteps)
    H_start = heat_content(state, c, Float64[]); melt_energy = 0.0
    for s in 1:nsteps
        e = solve_flux!(state, scratch, F, dt)
        Ts[s] = state.Tsrf[1]; melt_rate[s] = e.melt_energy_available / (c.Lm * dt)
        melt_energy += e.melt_energy_available
    end
    t_onset = times[findfirst(>(0), melt_rate)]
    closure = F * times[end] - (heat_content(state, c, Float64[]) - H_start) - melt_energy
    record!("constant flux on cold snow", "melt onset (model / analytical)", t_onset / t_onset_exact, "1", "≈ 1")
    record!("constant flux on cold snow", "energy closure F·t − ΔH − E_melt", closure, "J m⁻²", "round-off")
    hours = times ./ 3600
    a = plot(hours, Ts .- c.T0; c=C2, label="Chion", title=panel_title(1),
        xlabel="Time (h)", ylabel="Surface temperature (°C)", legend=:bottomright)
    t_an = range(0, t_onset_exact; length=200)
    plot!(a, t_an ./ 3600, (Ti - c.T0) .+ 2F .* sqrt.(κ .* t_an ./ pi) ./ K; c=C1, ls=:dash, label="Analytical (until onset)")
    vline!(a, [t_onset_exact / 3600]; c=INK, ls=:dot, lw=0.8, label="Analytical onset")
    b = plot(hours, melt_rate .* SECONDS_PER_DAY; c=C2, label="Chion", title=panel_title(2),
        xlabel="Time (h)", ylabel="Melt rate (kg m⁻² d⁻¹)", legend=:bottomright)
    hline!(b, [F / c.Lm * SECONDS_PER_DAY]; c=C1, ls=:dash, label="F / Lₘ (limit)")
    panel_text!(a, "Surface warming under\nconstant flux"; at=(0.96, 0.85), halign=:right)
    panel_text!(b, "Melt after onset"; at=(0.04, 0.5))
    push!(panels, a, b)

    # (c) Refreezing of added liquid water in a cold layer (cold-content limit).
    m, T_layer = 100.0, c.T0 - 10.0
    added = range(0, 15; length=61)
    refrozen = map(added) do W
        liquid, solid, rho, temp = [W], [m], [RHO_SNOW], [T_layer]
        Chion.go_refreezing!(liquid, solid, rho, temp, c.T0, c.ci, c.Lm, c.rho_i)
        solid[1] - m
    end
    exact = min.(added, c.ci * m * (c.T0 - T_layer) / c.Lm)
    record!("refreezing", "max |refrozen − min(W, cₚ m ΔT / Lₘ)|", maximum(abs.(refrozen .- exact)), "kg m⁻²", "round-off")
    cc = plot(added, exact; c=C1, label="min(W, cᵢ m (T₀ − T) / Lₘ)", title=panel_title(3),
        xlabel="Added liquid water W (kg m⁻²)", ylabel="Refrozen mass (kg m⁻²)", legend=:bottomright)
    scatter!(cc, added[1:3:end], refrozen[1:3:end]; c=C2, ms=3, msw=0, label="Chion")
    panel_text!(cc, "Refreezing in a cold layer"; at=(0.96, 0.6), halign=:right)
    push!(panels, cc)

    # (d) Cold-content payback in the ice substrate: cooling, then warming flux.
    Q1, Q2, days1, days2 = -20.0, 20.0, 30, 60
    dt = 3600.0; schedule = vcat(fill(Q1, days1 * 24), fill(Q2, days2 * 24))
    cum(model) = begin
        st = column_state(model, Float64[], c.T0); st.N[1] = 0
        sc = Chion.EnergyWorkspace(st); total = 0.0; out = zeros(length(schedule))
        for (s, Q) in enumerate(schedule)
            e = solve_flux!(st, sc, Q, dt); total += e.melt_energy_available / c.Lm; out[s] = total
        end
        out
    end
    with_substrate = cum(column_model(1; nice=5))
    # Without a substrate the model step holds bare ice at T0 and melts with max(Q, 0).
    without = cumsum(max.(schedule, 0.0) .* dt ./ c.Lm)
    final_exact = (Q1 * days1 + Q2 * days2) * SECONDS_PER_DAY / c.Lm
    record!("ice-substrate cold content", "final melt (model / energy budget)", with_substrate[end] / final_exact, "1", "→ 1 once the ice is back at T0")
    t_days = (1:length(schedule)) .* dt ./ SECONDS_PER_DAY
    d = plot(t_days, with_substrate; c=C2, label="With ice substrate", title=panel_title(4),
        xlabel="Time (days)", ylabel="Cumulative melt (kg m⁻²)", legend=:topleft)
    plot!(d, t_days, without; c=C3, label="Ice held at T₀ (no substrate)")
    hline!(d, [final_exact]; c=C1, ls=:dash, label="Energy budget (Q₁t₁ + Q₂t₂)/Lₘ")
    vline!(d, [days1]; c=INK, ls=:dot, lw=0.8, label="Flux reversal")
    panel_text!(d, "Cold content\nof bare ice"; at=(0.04, 0.45))
    push!(panels, d)
    save_figure(plot(panels...; layout=(2, 2), size=(720, 560), left_margin=4Plots.mm, bottom_margin=3Plots.mm), "fig5_phase_change")
end

# ---------- 4. discretisation ----------
"""Synthetic seasonal forcing (daily values) for one column, repeated each year."""
function synthetic_daily_forcing()
    d = 1:365
    season = sin.(2pi .* (d .- 105) ./ 365)
    Ta = 263.15 .+ 11.0 .* season
    swn = max.(0.0, 140.0 .* season)
    lwd = 0.80 * 5.670373e-8 .* Ta .^ 4 .+ 15.0
    shf = 3.0 .+ 6.0 .* max.(season, 0.0)
    snow = ifelse.(Ta .< 273.15, 2.5, 0.0); rain = ifelse.(Ta .>= 273.15, 1.0, 0.0)
    return (; Ta, swn, lwd, shf, snow, rain)
end

function column_run(substeps_per_day; nearsurface_scale=1.0, years=4)
    f = synthetic_daily_forcing()
    rep(v) = repeat(v; inner=substeps_per_day)
    forcing = SnowpackForcing(dt_days=fill(1.0 / substeps_per_day, 365 * substeps_per_day),
        air_temperature=rep(f.Ta), snowfall_rate=rep(f.snow) ./ SECONDS_PER_DAY, rainfall_rate=rep(f.rain) ./ SECONDS_PER_DAY,
        shortwave_down=rep(f.swn), q_sw_net=rep(f.swn), q_lw_down=rep(f.lwd), q_sh=rep(f.shf), q_lh=zeros(365 * substeps_per_day))
    model = BESSIModel(SnowpackGrid(1); diurnal_shortwave_substeps=false, diurnal_temperature_cycle=false,
        near_surface_layer_max_thicknesses_m=Tuple(nearsurface_scale .* (0.02, 0.05, 0.10, 0.30)))
    spin = Simulation(model; forcing, years=years - 1, backend=:threads, write_netcdf=false, compute_year_metrics=false)
    run!(spin; io=devnull)
    before = deepcopy(spin.now)
    sim = Simulation(model; forcing, state=spin.now, years=1, backend=:threads, write_netcdf=false, compute_year_metrics=false)
    run!(sim; io=devnull)
    delta(n) = getproperty(sim.now, n)[1] - getproperty(before, n)[1]
    return (melt=delta(:melt), runoff=delta(:runoff), refreezing=delta(:refreezing))
end

function run_discretisation()
    per_day = [1, 2, 4, 8, 24, 96]
    time_runs = [column_run(n) for n in per_day]
    scales = [4.0, 2.0, 1.0, 0.5, 0.25, 0.125]
    space_runs = [column_run(1; nearsurface_scale=s) for s in scales]
    terms = ((:melt, "Melt", C1), (:runoff, "Runoff", C2), (:refreezing, "Refreezing", C3))
    a = plot(; xscale=:log10, yscale=:log10, title=panel_title(1), xlabel="Time step (days)",
        ylabel="|Δ annual total| vs finest (%)", legend=:bottomright)
    b = plot(; xscale=:log10, yscale=:log10, title=panel_title(2),
        xlabel="Top-layer thickness (m)", ylabel="", legend=false)
    anchor = Ref(NaN)
    for (name, label, colour) in terms
        ref_t = getproperty(time_runs[end], name); ref_s = getproperty(space_runs[end], name)
        (ref_t > 0 && ref_s > 0) || (@warn "No $label in the reference run; skipped"; continue)
        dt_err = [abs(getproperty(r, name) - ref_t) / ref_t * 100 for r in time_runs[1:end-1]]
        sp_err = [abs(getproperty(r, name) - ref_s) / ref_s * 100 for r in space_runs[1:end-1]]
        plot!(a, 1 ./ per_day[1:end-1], max.(dt_err, 1e-6); c=colour, marker=:circle, ms=3.5, msw=0, label=label)
        plot!(b, 0.02 .* scales[1:end-1], max.(sp_err, 1e-6); c=colour, marker=:circle, ms=3.5, msw=0, label=label)
        name == :melt && (anchor[] = dt_err[end])
        record!("time-step convergence (column, daily forcing)", "$label: reference annual total (Δt = 1/96 d)", ref_t, "kg m⁻²", "")
        record!("time-step convergence (column, daily forcing)", "$label: |Δ| at daily step vs 1/96 d", dt_err[1], "%", "small, decreasing with Δt")
        record!("near-surface resolution", "$label: reference annual total (top layer 0.0025 m)", ref_s, "kg m⁻²", "")
        record!("near-surface resolution", "$label: |Δ| at default (0.02 m) vs 0.0025 m", sp_err[3], "%", "small, decreasing with Δz")
    end
    hs = 1 ./ per_day[1:end-1]
    plot!(a, hs, anchor[] .* hs ./ hs[end]; c=INK, ls=:dash, lw=1, label="First order")
    panel_text!(a, "Time step"); panel_text!(b, "Near-surface layer\nthickness"; at=(0.96, 0.4), halign=:right)
    save_figure(plot(a, b; layout=(1, 2), size=(720, 320), left_margin=4Plots.mm, bottom_margin=4Plots.mm), "fig6_discretisation")
    return time_runs, space_runs
end

# ---------- outputs ----------
function write_metrics()
    open(joinpath(VOUT, "metrics.csv"), "w") do io
        println(io, "test,metric,value,unit,expectation")
        for r in METRICS
            println(io, join(('"' * r.test * '"', '"' * r.metric * '"', r.value, r.unit, '"' * r.expectation * '"'), ','))
        end
    end
end

function verification_main()
    mkpath(VOUT)
    @info "1. heat conduction"; run_step_response(); run_periodic(); run_convergence()
    @info "2. conservation"; mass_res, step_energy_res, ncol = run_full_step_conservation(); solve_energy_res = run_energy_closure()
    plot_conservation(mass_res, solve_energy_res)
    @info "3. phase change"; run_phase_change()
    @info "4. discretisation"; run_discretisation()
    write_metrics()
    @info "Verification written" VOUT
end

if abspath(PROGRAM_FILE) == @__FILE__
    verification_main()
end
