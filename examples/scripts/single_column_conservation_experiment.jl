#!/usr/bin/env julia

"""Single-column BESSI mass and energy conservation experiment."""

using Chion, Plots, Printf

const MASS_OUTPUT = joinpath(@__DIR__, "..", "plots", "single_column_mass_conservation.pdf")
const ENERGY_OUTPUT = joinpath(@__DIR__, "..", "plots", "single_column_energy_conservation.pdf")
const COLUMN_OUTPUT = joinpath(@__DIR__, "..", "plots", "single_column_snowpack_evolution.pdf")

# Reduced enthalpy relative to ice at 0 °C. BESSI assigns liquid water Lm.
enthalpy(s) = sum((s.mass[k, 1] * s.c.ci * (s.temperature[k, 1] - s.c.T0) +
                   s.mass_w[k, 1] * s.c.Lm for k in 1:s.N[1]); init=0.0)
stored_mass(s) = sum(@view s.mass[:, 1]) + sum(@view s.mass_w[:, 1])

function audit!(f, records, step, process, state)
    before = enthalpy(state)
    expected, result = f()
    push!(records, (; step, process, expected, result,
        residual=enthalpy(state) - before - expected))
    return result
end

"One BESSI step, split into processes so each can be audited."
function audited_step!(s, forcing, step, workspace, records)
    i, c = 1, s.c
    f = Chion._step_forcing_at(forcing, i, step)
    dt = f.dt_days * c.seconds_per_day
    empty = s.N[i] == 0 || s.mass[1, i] <= Chion.EPS_EMPTY_LAYER

    audit!(records, step, :accumulation, s) do
        H0 = enthalpy(s)
        runoff0, base0 = s.runoff[i], s.mass_base[i]
        Tin = empty ? f.air_temperature : s.temperature[1, i]
        Chion._apply_accumulation_resolved!(s.N, s.mass, s.mass_w, s.density,
            s.temperature, s.mass_base, s.smb_ice, s.runoff, s.Tsrf, s.albedo,
            i, c, s.Ntot, s.mass_max, s.mass_split, s.mass_min,
            f.snowfall_rate, f.rainfall_rate, dt, f.air_temperature, f.wind_speed)
        empty && f.snowfall_rate > 0 && (s.temperature[1, i] = f.air_temperature)
        Qprecip = dt * (f.snowfall_rate * c.ci * (Tin - c.T0) + f.rainfall_rate * c.Lm)
        expected = Qprecip - (s.runoff[i] - runoff0) * c.Lm
        # Diagnose the enthalpy carried by basal ice from the exact discrete
        # storage change when the standard depth cap exports solid mass.
        if s.mass_base[i] > base0
            Qbase = H0 + expected - enthalpy(s)
            expected -= Qbase
        end
        (expected, nothing)
    end

    Chion._update_surface_albedo_arrays!(s.N, s.mass, s.mass_w, s.density,
        s.temperature, s.albedo, i, c)
    audit!(records, step, :densification, s) do
        accumulation = max(f.snowfall_rate, 0.0) + f.rainfall_rate
        Chion._go_densification!(s.N, s.mass, s.density, s.temperature,
            i, c, accumulation, dt)
        (0.0, nothing)
    end

    linear, constant = Chion._diagnose_latent_heat_flux_coefficients(
        true, c, f.air_temperature, f.snowfall_rate, f.rainfall_rate)
    energy = audit!(records, step, :energy_solve, s) do
        result = Chion._go_energy_flux_resolved!(s.N, s.mass, s.mass_w,
            s.density, s.temperature, s.Tsrf, s.albedo, i, c, workspace.energy,
            f.air_temperature, f.shortwave_down, linear, constant, dt,
            f.has_q_sw_net, f.q_sw_net, f.has_q_lw_down, f.q_lw_down,
            f.has_q_sh, f.q_sh, f.has_q_lh, f.q_lh,
            f.has_relative_humidity, f.relative_humidity, f.air_pressure)
        (result.heating, result)
    end

    audit!(records, step, :vapor_exchange, s) do
        surface_temperature = s.temperature[1, i]
        solid_phase = surface_temperature < c.T0
        runoff0 = s.runoff[i]
        fluxes = Chion._apply_snow_surface_vapor_mass_flux!(s.N, s.mass,
            s.mass_w, s.density, s.temperature, s.runoff, s.Tsrf, s.albedo,
            i, c, f, dt, s.mass_split, s.mass_min)
        s.vapor_mass[i] += fluxes.vapor_mass
        s.sublimation[i] += fluxes.sublimation_mass
        s.latent_heat_flux_sum[i] += fluxes.latent_heat_flux * f.dt_days
        solid_vapor = solid_phase ? fluxes.vapor_mass : 0.0
        liquid_vapor = solid_phase ? 0.0 : fluxes.vapor_mass
        carrier_heat = solid_vapor * c.ci * (surface_temperature - c.T0) +
                       liquid_vapor * c.Lm - (s.runoff[i] - runoff0) * c.Lm
        (carrier_heat, (; solid_vapor, liquid_vapor))
    end

    if energy.needs_melt
        audit!(records, step, :melt, s) do
            runoff0 = s.runoff[i]
            dm = energy.melt_energy_available / c.Lm
            melted = Chion._apply_melt!(s.N, s.mass, s.mass_w, s.density,
                s.temperature, s.runoff, s.Tsrf, s.albedo, i,
                s.mass_split, s.mass_min, dm, c)
            melted == dm || error("audit unexpectedly reached bare ice")
            s.melt[i] += dm
            (energy.melt_energy_available - (s.runoff[i] - runoff0) * c.Lm, nothing)
        end
    end

    if Chion._column_has_liquid_water(s.N, s.mass_w, i)
        audit!(records, step, :percolation, s) do
            runoff0 = s.runoff[i]
            s.runoff[i] += Chion._go_percolation!(s.N, s.mass, s.mass_w,
                s.density, i, c.rho_i, c.rho_w)
            (-(s.runoff[i] - runoff0) * c.Lm, nothing)
        end
    end
    if Chion._column_has_liquid_water(s.N, s.mass_w, i)
        audit!(records, step, :refreezing, s) do
            frozen = Chion._go_refreezing!(s.N, s.mass_w, s.mass, s.density,
                s.temperature, i, c.T0, c.ci, c.Lm, c.rho_i)
            s.refreezing[i] += frozen
            (0.0, nothing)
        end
    end
end

function save_state!(h, s)
    push!(h.mass, stored_mass(s)); push!(h.heat, enthalpy(s))
    push!(h.runoff, s.runoff[1]); push!(h.base, s.mass_base[1])
    push!(h.vapor, s.vapor_mass[1]); push!(h.sublimation, s.sublimation[1])
    push!(h.melt, s.melt[1])
    push!(h.freeze, s.refreezing[1]); push!(h.solid, copy(s.mass[:, 1]))
    push!(h.liquid, copy(s.mass_w[:, 1])); push!(h.density, copy(s.density[:, 1]))
end

function panel!(p, letter)
    x, y = xlims(p), ylims(p)
    annotate!(p, x[1] - 0.06(x[2] - x[1]), y[2] + 0.03(y[2] - y[1]),
        text(letter, 14, :black, :left))
end

function reserve_legend_space!(p; rows=1)
    lo, hi = ylims(p)
    ylims!(p, lo - (0.18 + 0.12(rows - 1)) * (hi - lo), hi)
end

function print_audit(records)
    println("\nProcess energy audit")
    @printf("%-14s %14s %14s\n", "process", "sum residual", "max |residual|")
    for p in unique(r.process for r in records)
        selected = [r for r in records if r.process == p]
        @printf("%-14s %+14.6e %14.6e\n", String(p),
            sum(r.residual for r in selected), maximum(abs(r.residual) for r in selected))
    end
    println("\nNegative boundary-energy contributions")
    for p in unique(r.process for r in records)
        q = -sum(min(r.expected, 0.0) for r in records if r.process == p)
        q > 0 && @printf("%-14s %14.6e J m^-2\n", String(p), q)
    end
end

function main()
    nt, nspinup, nplot = 365, 365*100, 2*365
    nsteps = nspinup + nplot
    day = 1:nt
    season = cos.(2π .* (day .- 200) ./ nt)
    Tcycle = -5 .+ 7season
    precipitation = 6 .+ 2season
    snow_fraction = ifelse.(Tcycle .<= 0, 1.0,
        ifelse.(Tcycle .>= 2, 0.0,
            1 .- Tcycle ./ 2))
    snowcycle = precipitation .* snow_fraction
    raincycle = precipitation .* (1 .- snow_fraction)*0.05
    swcycle = 20 .+ 120max.(season, 0.0)
    lwcycle = 250 .+ 45season
    shcycle = -12 .+ 14season
    lhcycle = -3 .+ 2season
    ncycles = cld(nsteps, nt)
    Tair = repeat(Tcycle, ncycles)[1:nsteps]
    snow = repeat(snowcycle, ncycles)[1:nsteps]
    rain = repeat(raincycle, ncycles)[1:nsteps]
    shortwave = repeat(swcycle, ncycles)[1:nsteps]
    lw = repeat(lwcycle, ncycles)[1:nsteps]
    sh = repeat(shcycle, ncycles)[1:nsteps]
    lh = repeat(lhcycle, ncycles)[1:nsteps]
    dt_days = fill(1.0, nsteps)
    forcing = SnowpackForcing(; dt_days, ncol=1, air_temperature_c=Tair,
        snowfall_mm_day=snow, rainfall_mm_day=rain, shortwave_down=shortwave,
        q_lw_down=lw, q_sh=sh, q_lh=lh)
    model = BESSIModel(SnowpackGrid(1))
    sim = Simulation(model; forcing, years=1, backend=:threads,
        write_netcdf=false, history_year_stride=0, compute_year_metrics=false)
    workspace, records = Chion.ColumnarStepWorkspace(sim.now), NamedTuple[]
    h = (mass=Float64[], heat=Float64[], runoff=Float64[], base=Float64[],
        vapor=Float64[], sublimation=Float64[], melt=Float64[], freeze=Float64[],
        solid=Vector{Vector{Float64}}(), liquid=Vector{Vector{Float64}}(),
        density=Vector{Vector{Float64}}())
    save_state!(h, sim.now)
    for step in 1:nsteps
        audited_step!(sim.now, forcing, step, workspace, records)
        save_state!(h, sim.now)
    end

    steps = 0:nsteps
    mass_in = [0.0; cumsum((snow + rain) .* dt_days)]
    mass_error = h.mass + h.runoff + h.base - mass_in - h.vapor
    qin_step, qout_step = zeros(nsteps), zeros(nsteps)
    for r in records
        qin_step[r.step] += max(r.expected, 0.0)
        qout_step[r.step] += max(-r.expected, 0.0)
    end
    Qin, Qout = [0.0; cumsum(qin_step)], [0.0; cumsum(qout_step)]
    Qnet, heat_error = Qin - Qout, h.heat - (Qin - Qout)
    scale = max(maximum(Qin), maximum(Qout), 1.0)

    @printf("Mass residual: %+.3e kg m^-2\n", mass_error[end])
    @printf("Heat residual: %+.3e J m^-2 (%.3e %%)\n", heat_error[end],
        100heat_error[end] / scale)
    @printf("Melt %.1f, refreezing %.1f, runoff %.1f, basal export %.1f kg m^-2\n",
        h.melt[end], h.freeze[end], h.runoff[end], h.base[end])
    @printf("Sublimation %.1f kg m^-2, net vapor mass %+.1f kg m^-2\n",
        h.sublimation[end], h.vapor[end])
    print_audit(records)

    shown = (nspinup + 1):nsteps
    before, after = (nspinup + 1):nsteps, (nspinup + 2):(nsteps + 1)
    dt = dt_days[shown]
    solid = sum.(h.solid)
    liquid = sum.(h.liquid)
    effective_solid = (solid[after] - solid[before]) ./ dt
    effective_liquid = (liquid[after] - liquid[before]) ./ dt
    refreeze = (h.freeze[after] - h.freeze[before]) ./ dt
    melt = (h.melt[after] - h.melt[before]) ./ dt
    basal = (h.base[after] - h.base[before]) ./ dt
    runoff = (h.runoff[after] - h.runoff[before]) ./ dt
    vapor_solid = zeros(length(shown)); vapor_liquid = similar(vapor_solid)
    for r in records
        if r.process == :vapor_exchange && r.step in shown
            j = r.step - nspinup
            vapor_solid[j] += r.result.solid_vapor / dt[j]
            vapor_liquid[j] += r.result.liquid_vapor / dt[j]
        end
    end
    diagnosed_solid = snow[shown] + refreeze - melt - basal + vapor_solid
    diagnosed_liquid = rain[shown] + melt - refreeze - runoff + vapor_liquid

    tol = palette(:tol_light)
    blue, coral, yellow, green, grey = tol[1], tol[2], tol[3], tol[6], tol[9]
    common = (; grid=false, linewidth=2, framestyle=:box, palette=tol,
        foreground_color_axis=:black, background_color_legend=RGBA(1, 1, 1, 0.78),
        legendfontsize=10,
        tickfontsize=11, guidefontsize=12)
    blank_x = _ -> ""

    m1 = plot(shown, effective_solid; label="Effective solid mass flux", color=:black,
        ylabel="Effective Solid Mass Flux\n(kg m⁻² d⁻¹)", legend=:bottomright,
        xformatter=blank_x, common...)
    plot!(m1, shown, effective_liquid; label="Effective liquid mass flux",
        color=grey, linewidth=3); reserve_legend_space!(m1); panel!(m1, "(a)")

    m2 = plot(shown, snow[shown]; label="Snowfall", color=coral,
        ylabel="Solid Mass Flux\n(kg m⁻² d⁻¹)", legend=:bottomright,
        xformatter=blank_x, common...)
    plot!(m2, shown, refreeze; label="Refreezing", color=blue)
    plot!(m2, shown, -melt; label="Melt", color=green)
    plot!(m2, shown, -basal; label="Basal export", color=yellow)
    plot!(m2, shown, vapor_solid; label="Vapor", color=grey)
    plot!(m2, shown, diagnosed_solid; label="Diagnosed solid flux", color=:black,
        linewidth=2.5)
    plot!(m2; legend_column=3); reserve_legend_space!(m2; rows=2); panel!(m2, "(b)")

    solid_error = effective_solid - diagnosed_solid
    m3 = plot(shown, 1e12 .* solid_error; label="Effective − diagnosed solid flux",
        color=:black, ylabel="Solid Residual\n(×10⁻¹² kg m⁻² d⁻¹)", legend=:bottomright,
        xformatter=blank_x, common...)
    reserve_legend_space!(m3); panel!(m3, "(c)")

    m4 = plot(shown, rain[shown]; label="Rainfall", color=coral,
        ylabel="Liquid Mass Flux\n(kg m⁻² d⁻¹)", legend=:bottomright,
        xformatter=blank_x, common...)
    plot!(m4, shown, -refreeze; label="Refreezing", color=blue)
    plot!(m4, shown, melt; label="Melt", color=green)
    plot!(m4, shown, -runoff; label="Runoff", color=yellow)
    plot!(m4, shown, vapor_liquid; label="Vapor", color=grey)
    plot!(m4, shown, diagnosed_liquid; label="Diagnosed liquid flux", color=:black,
        linewidth=2.5)
    plot!(m4; legend_column=3); reserve_legend_space!(m4; rows=2); panel!(m4, "(d)")

    liquid_error = effective_liquid - diagnosed_liquid
    m5 = plot(shown, 1e12 .* liquid_error; label="Effective − diagnosed liquid flux",
        color=:black, ylabel="Liquid Residual\n(×10⁻¹² kg m⁻² d⁻¹)", xlabel="forcing step",
        legend=:bottomright, common...)
    reserve_legend_space!(m5); panel!(m5, "(e)")

    mass_figure = plot(m1, m2, m3, m4, m5; layout=(5, 1), link=:x,
        size=(1200, 1800), left_margin=16Plots.mm, right_margin=4Plots.mm,
        top_margin=6Plots.mm, bottom_margin=2Plots.mm)

    diagnosed_energy = (qin_step[shown] - qout_step[shown]) ./ dt
    effective_energy = (h.heat[after] - h.heat[before]) ./ dt
    energy_error = effective_energy - diagnosed_energy
    e1 = plot(shown, 1e-6 .* diagnosed_energy; label="Diagnosed total energy flux", color=coral,
        ylabel="Total Energy Flux\n(×10⁶ J m⁻² d⁻¹)", legend=:bottomright,
        xformatter=blank_x, common...)
    plot!(e1, shown, 1e-6 .* effective_energy; label="Effective total energy flux",
        color=:black, linewidth=2.5)
    reserve_legend_space!(e1); panel!(e1, "(a)")
    e2 = plot(shown, 1e6 .* energy_error; label="Effective − diagnosed energy flux",
        color=:black, ylabel="Energy Residual\n(×10⁻⁶ J m⁻² d⁻¹)", xlabel="forcing step",
        legend=:bottomright, common...)
    reserve_legend_space!(e2); panel!(e2, "(b)")
    energy_figure = plot(e1, e2; layout=(2, 1), link=:x, size=(1200, 750),
        left_margin=16Plots.mm, right_margin=4Plots.mm,
        top_margin=6Plots.mm, bottom_margin=2Plots.mm)

    # Column geometry: the thick line is the surface and the thin lines are
    # boundaries between individual snow layers.
    nmax, nshow = length(h.solid[after[1]]), length(shown)
    interfaces = fill(NaN, nmax, nshow)
    for (j, state_index) in enumerate(after)
        nactive = something(findlast(>(Chion.EPS_EMPTY_LAYER), h.solid[state_index]), 0)
        height = 0.0
        for k in nactive:-1:1
            height += h.solid[state_index][k] / h.density[state_index][k]
            interfaces[k, j] = height
        end
    end
    surface = interfaces[1, :]
    z = range(0, maximum(surface); length=160)
    liquid_fraction = fill(NaN, length(z), nshow)
    for (j, state_index) in enumerate(after)
        nactive = something(findlast(>(Chion.EPS_EMPTY_LAYER), h.solid[state_index]), 0)
        for k in 1:nactive
            bottom = k == nactive ? 0.0 : interfaces[k + 1, j]
            fraction = h.liquid[state_index][k] /
                       (h.solid[state_index][k] + h.liquid[state_index][k])
            liquid_fraction[(bottom .<= z) .& (z .<= interfaces[k, j]), j] .= fraction
        end
    end
    wetmax = maximum(filter(isfinite, liquid_fraction))
    c1 = heatmap(shown, z, liquid_fraction; color=cgrad([RGBA(0.88, 0.94, 0.98, 1), coral]),
        clims=(0, wetmax), colorbar_title="Liquid Mass Fraction", ylabel="Snowpack Height (m)",
        xformatter=blank_x, colorbar=true, common...)
    for k in 2:nmax
        plot!(c1, shown, interfaces[k, :]; color=blue, linewidth=0.7,
            alpha=0.65, label="")
    end
    plot!(c1, shown, surface; color=blue, linewidth=3, label="")
    panel!(c1, "(a)")

    c2 = plot(shown, solid[after]; color=blue, label="Solid mass",
        ylabel="Solid Mass (kg m⁻²)", xformatter=blank_x,
        legend=:bottomright, common...)
    reserve_legend_space!(c2); panel!(c2, "(b)")
    c3 = plot(shown, liquid[after]; color=coral, linewidth=2.5,
        label="Liquid mass", ylabel="Liquid Mass (kg m⁻²)", xlabel="forcing step",
        legend=:bottomright, common...)
    reserve_legend_space!(c3); panel!(c3, "(c)")
    column_figure = plot(c1, c2, c3; layout=(3, 1), link=:x, size=(1200, 1200),
        left_margin=16Plots.mm, right_margin=18Plots.mm,
        top_margin=6Plots.mm, bottom_margin=2Plots.mm)

    mkpath(dirname(MASS_OUTPUT))
    savefig(mass_figure, MASS_OUTPUT); savefig(energy_figure, ENERGY_OUTPUT)
    savefig(column_figure, COLUMN_OUTPUT)
    println("Wrote $MASS_OUTPUT\nWrote $ENERGY_OUTPUT\nWrote $COLUMN_OUTPUT")
end

main()
