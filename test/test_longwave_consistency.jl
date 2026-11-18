using Test, Dates, Chion

_one_column_forcing() = SnowpackForcing(dt_days = [1.0], ncol = 1, air_temperature_c = -10.0, snowfall_mm_day = 0.0,
    rainfall_mm_day = 0.0, shortwave_down = 0.0, latitude_deg = 70.0)

@testset "Thin-layer Robin boundary preserves insulated equilibrium" begin
    params = Params(ϵ_snow = 0.0)    # the fluxes are prescribed; no longwave emission
    m = BESSIModel(Grid(1; layering = Layering(Ntot = 1, ice_substrate_layers = 0)); forcing = _one_column_forcing(),
        params)
    for mass in (100., 1e-8)
        s = m.state; s.N[1] = 1; s.mass[1, 1] = mass; s.density[1, 1] = 300.
        s.temperature[1, 1] = 260.; s.Tsrf[1] = 260.
        e = Chion.go_energy_flux!(m, 1, 260., 0., 0., 0., 10800.;
            q_sw_net = 0., q_lw_down = 0., q_sh = 0., q_lh = 0.)
        @test s.Tsrf[1] ≈ 260. atol = 1e-10
        @test s.temperature[1, 1] ≈ 260. atol = 1e-10
        # Zero flux derivative: the surface-row diagonal is exactly one, so the
        # resolved profile reproduces the old temperature bit for bit.
        @test m.cache.workspace.energy.rhs[1, 1] == 260.
        @test e.melt_energy_available == 0.
    end
end

@testset "Graybody longwave equilibrium with prescribed incident radiation" begin
    grid = Grid(1; layering = Layering(Ntot = 2, ice_substrate_layers = 0))
    gray = BESSIModel(grid; forcing = _one_column_forcing(), params = Params(ϵ_snow = .98, eps_ice = .98))
    legacy = Params(surface_energy_balance_scheme = BESSISurfaceEnergyBalance(), ϵ_snow = .98, eps_ice = .98)
    c = gray.params; dt = 10800.; incident = c.σ * c.T0^4
    f = Chion.SnowpackStepForcing(c.T0, .125, 0., 0., 0., 0.;
        q_sw_net = 0., has_q_sw_net = true, q_lw_down = incident, has_q_lw_down = true,
        q_sh = 0., has_q_sh = true, q_lh = 0., has_q_lh = true)
    # An isothermal blackbody radiation environment cannot melt gray ice:
    # absorbed incident LW and emitted LW must cancel, regardless of epsilon.
    bare = Chion._bare_ice_ablation_mass(c, f, dt)
    @test abs(bare.longwave_flux) < 1e-10
    @test abs(bare.melt_mass) < 1e-10
    old = Chion._bare_ice_ablation_mass(legacy, f, dt)
    @test old.longwave_flux ≈ (1 - c.ϵ_snow) * incident
    @test old.melt_mass > 0
    s = gray.state; s.N[1] = 2; s.mass[:, 1] .= (100., 300.)
    s.density[:, 1] .= (300., 500.); s.temperature[:, 1] .= c.T0; s.Tsrf[1] = c.T0
    energy = Chion.go_energy_flux!(gray, 1, c.T0, 0., 0., 0., dt;
        q_sw_net = 0., q_lw_down = incident, q_sh = 0., q_lh = 0.)
    @test s.Tsrf[1] ≈ c.T0 atol = 1e-10
    @test all(abs.(s.temperature[:, 1] .- c.T0) .< 1e-10)
    @test abs(energy.melt_energy_available) < 1e-5
end

@testset "Clear-sky transmissivity of the cloud proxy" begin
    c = Params()
    lw = c.longwave_scheme
    τ(z, toa) = min(1.0, lw.clear_sky_transmissivity_sea_level + lw.clear_sky_transmissivity_per_km * z / 1000 +
        lw.clear_sky_transmissivity_insolation * log(toa / 300))
    @test clear_sky_transmissivity(c, 1500.0, 450.0) ≈ τ(1500.0, 450.0) rtol = 1e-12
    @test clear_sky_transmissivity(c, NaN, 450.0) == clear_sky_transmissivity(c, 0.0, 450.0)
    @test clear_sky_transmissivity(c, -200.0, 450.0) == clear_sky_transmissivity(c, 0.0, 450.0)
    @test clear_sky_transmissivity(c, 1500.0, 10.0) == clear_sky_transmissivity(c, 1500.0, Chion.LONGWAVE_CLOUD_PROXY_MIN_TOA)
    @test clear_sky_transmissivity(c, 1.0e5, 450.0) == 1.0

    # A host's SW↓ from a cloud fraction n through the proxy's cloud transmission
    # gives the proxy back n.
    n = 0.6
    z = 1500.0
    forcing = SnowpackForcing(dt_days = [1.0], air_temperature = [268.0], snowfall_rate = [0.0],
        rainfall_rate = [0.0], shortwave_down = [0.0], latitude_deg = 70.0, surface_height = z,
        time_values = [DateTime(2001, 6, 15, 12)])
    step_forcing = Chion._step_forcing_at(forcing, 1, 1)
    terms = Chion._diurnal_shortwave_integral_terms(step_forcing.latitude_deg, step_forcing.solar_longitude_deg)
    toa = Chion._daily_toa_shortwave(terms, step_forcing.day_of_year)
    @test toa > Chion.LONGWAVE_CLOUD_PROXY_MIN_TOA
    cloud_effect = lw.cloud_transmission_coefficient * exp(-lw.cloud_transmission_decay_per_km * z / 1000)
    shortwave_down = toa * clear_sky_transmissivity(c, z, toa) * (1 - cloud_effect * n^lw.cloud_transmission_exponent)
    step_forcing = Chion._with_step_forcing(step_forcing, (; shortwave_down))
    @test Chion._cloud_proxy_cloud_fraction(c, step_forcing, terms) ≈ n rtol = 1e-10
end
