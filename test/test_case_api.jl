using Test
using Dates
using NCDatasets
using Chion

function _captured_exception(f::Function)
    try
        f()
        return nothing
    catch err
        return err
    end
end

# A one-step forcing for models that are only used for column-level tests.
_placeholder_forcing(ncol::Integer = 1) = SnowpackForcing(
    dt_days           = [1.0],
    ncol              = ncol,
    air_temperature_c = -10.0,
    snowfall_mm_day   = 0.0,
    rainfall_mm_day   = 0.0,
    shortwave_down    = 0.0,
    latitude_deg      = 70.0,
)

_bessi_model(grid::Grid; params = Params(; NF = Chion.number_type(grid)), forcing = _placeholder_forcing(grid.ncol)) =
    BESSIModel(grid; forcing, params)

function _sample_forcing()
    return SnowpackForcing(
        dt_days           = [1.0, 1.0, 1.0],
        ncol              = 4,
        air_temperature_c = [-15.0, -14.0, -13.0],
        snowfall_mm_day   = 1.0,
        rainfall_mm_day   = 0.0,
        shortwave_down    = [120.0, 140.0, 160.0],
        wind_speed        = 3.0,
        latitude_deg      = 70.0,
        time_values       = [DateTime(2001, 1, 15, 12), DateTime(2001, 2, 15, 12), DateTime(2001, 3, 15, 12)],
    )
end

function _sample_grid()
    return Grid(4;
        layering = Layering(Ntot = 4),
        x    = [0.0, 10_000.0],
        y    = [0.0, 10_000.0],
        js   = [1, 1, 2, 2],
        is   = [1, 2, 1, 2],
        mask = ones(2, 2),
    )
end

function _sample_model_forcing_grid()
    grid = _sample_grid()
    forcing = _sample_forcing()
    model = BESSIModel(grid; forcing)
    return model, forcing, grid
end

function _write_sample_forcing_file(path::AbstractString)
    nx = 2; ny = 2; ntime = 2
    ds = NCDataset(path, "c")
    try
        defDim(ds, "x", nx)
        defDim(ds, "y", ny)
        defDim(ds, "time", ntime)

        defVar(ds, "x", Float64, ("x",))[:] = [0.0, 10_000.0]
        defVar(ds, "y", Float64, ("y",))[:] = [0.0, 10_000.0]
        defVar(ds, "MSK", Float64, ("x", "y"))[:, :] = [100.0 0.0; 75.0 25.0]

        time = defVar(ds, "time", Float64, ("time",))
        time.attrib["units"] = "days since 2001-01-01 12:00:00"
        time[:] = [0.0, 1.0]

        tt = fill(-15.0, ntime, ny, nx); tt[2, :, :] .= -14.0
        defVar(ds, "TT", Float64, ("time", "y", "x"))[:, :, :] = tt
        defVar(ds, "SF", Float64, ("time", "y", "x"))[:, :, :] = fill(1.0, ntime, ny, nx)
        defVar(ds, "RF", Float64, ("time", "y", "x"))[:, :, :] = fill(0.0, ntime, ny, nx)
        defVar(ds, "SWD", Float64, ("time", "y", "x"))[:, :, :] = fill(100.0, ntime, ny, nx)
        defVar(ds, "LWD", Float64, ("time", "y", "x"))[:, :, :] = fill(250.0, ntime, ny, nx)
        defVar(ds, "SHF", Float64, ("time", "y", "x"))[:, :, :] = fill(0.0, ntime, ny, nx)
        defVar(ds, "LHF", Float64, ("time", "y", "x"))[:, :, :] = fill(0.0, ntime, ny, nx)
    finally
        close(ds)
    end
    return path
end

function _write_sample_forcing_file_4d(path::AbstractString)
    nx = 2; ny = 2; ntime = 2; nlayer = 1
    ds = NCDataset(path, "c")
    try
        defDim(ds, "TIME", ntime)
        defDim(ds, "layer", nlayer)
        defDim(ds, "y", ny)
        defDim(ds, "x", nx)

        defVar(ds, "x", Float64, ("x",))[:] = [0.0, 10_000.0]
        defVar(ds, "y", Float64, ("y",))[:] = [0.0, 10_000.0]
        defVar(ds, "TIME", Float64, ("TIME",))[:] = [0.0, 1.0]
        defVar(ds, "YYYY", Float64, ("TIME",))[:] = [2001.0, 2001.0]
        defVar(ds, "MM", Float64, ("TIME",))[:] = [1.0, 1.0]
        defVar(ds, "DD", Float64, ("TIME",))[:] = [1.0, 2.0]
        defVar(ds, "HH", Float64, ("TIME",))[:] = [12.0, 12.0]

        tt = fill(-15.0, nx, ny, nlayer, ntime); tt[:, :, :, 2] .= -14.0
        defVar(ds, "TT", Float32, ("x", "y", "layer", "TIME"))[:, :, :, :] = tt
        defVar(ds, "SF", Float32, ("x", "y", "TIME"))[:, :, :] = fill(1.0, nx, ny, ntime)
        defVar(ds, "RF", Float32, ("x", "y", "TIME"))[:, :, :] = fill(0.0, nx, ny, ntime)
        defVar(ds, "SWD", Float32, ("x", "y", "TIME"))[:, :, :] = fill(100.0, nx, ny, ntime)
        defVar(ds, "LWD", Float32, ("x", "y", "TIME"))[:, :, :] = fill(250.0, nx, ny, ntime)
    finally
        close(ds)
    end
    return path
end

# A PDD model on `forcing` (one column) whose state starts from `state`, after all forcing steps.
function _pdd_after(params::PDDParams, forcing::SnowpackForcing; state = nothing, grid = Grid(1))
    model = PDDModel(grid; forcing, params, state)
    integrate!(model, 1:length(forcing.time_values))
    return model
end

@testset "Run API" begin
    @testset "Refreezing correction scales cold-content capacity" begin
        grid = Grid(1)
        @test _bessi_model(grid).refreezing_correction == 1.0
        @test _bessi_model(grid; params = Params(refreezing_correction=1.5)).refreezing_correction == 1.5
        @test_throws ArgumentError Params(refreezing_correction=0.0)

        function corrected_refreezing(factor)
            liquid = [20.0]
            solid = [100.0]
            density = [400.0]
            temperature = [263.15]
            result = Chion.go_refreezing!(
                liquid, solid, density, temperature,
                273.15, 2110.0, 334000.0, 917.0;
                refreezing_correction=factor,
            )
            return (; result..., liquid, solid, density, temperature)
        end

        baseline = corrected_refreezing(1.0)
        corrected = corrected_refreezing(1.5)
        @test corrected.refrozen_mass ≈ 1.5 * baseline.refrozen_mass
        @test corrected.refrozen_mass < 20.0
        @test corrected.liquid[1] + corrected.solid[1] ≈ 120.0
        @test corrected.temperature[1] == 273.15
    end

    @testset "SEMIX aerodynamic resistance follows atmospheric stability" begin
        c = Params()
        @test c.turbulent_flux_scheme === SemixTurbulentFluxes()
        @test c.surface_energy_balance_scheme === SemixSurfaceEnergyBalance()
        z0m_snow = c.turbulent_flux_scheme.z0m_snow
        neutral = Chion._semix_aerodynamic_resistance(c, 268.0, 268.0, 80_000.0, 5.0, z0m_snow)
        stable = Chion._semix_aerodynamic_resistance(c, 268.0, 273.0, 80_000.0, 5.0, z0m_snow)
        unstable = Chion._semix_aerodynamic_resistance(c, 273.0, 268.0, 80_000.0, 5.0, z0m_snow)
        @test stable > neutral > unstable
        weakly_suppressed = Params(turbulent_flux_scheme=SemixTurbulentFluxes(stable_coefficient=2.0))
        @test Chion._semix_aerodynamic_resistance(
            weakly_suppressed, 268.0, 273.0, 80_000.0, 5.0, weakly_suppressed.turbulent_flux_scheme.z0m_snow,
        ) < stable

        surface_temperature = 263.0
        air_temperature = 268.0
        relative_humidity = 0.8
        air_pressure = 80_000.0
        wind_speed = 5.0
        _, sensible, latent, _ = Chion._resolved_nonshortwave_surface_flux_components(
            c, air_temperature, 0.0, 86_400.0, surface_temperature,
            true, 0.0, false, 0.0, false, 0.0, true, relative_humidity,
            air_pressure, wind_speed, z0m_snow, c.ϵ_snow,
        )
        sensible_constant, sensible_linear, latent_constant, latent_linear =
            Chion._semix_turbulent_flux_linearized(
                surface_temperature, c, air_temperature, relative_humidity,
                air_pressure, wind_speed, z0m_snow,
            )
        @test sensible ≈ sensible_constant - sensible_linear * surface_temperature
        @test latent ≈ latent_constant - latent_linear * surface_temperature

        bessi = Params(turbulent_flux_scheme=BESSITurbulentFluxes())
        @test bessi.turbulent_flux_scheme === BESSITurbulentFluxes()
        # BESSI fluxes ignore the roughness length; pass the SEMIX snow value as before.
        _, bessi_sensible, bessi_latent, _ = Chion._resolved_nonshortwave_surface_flux_components(
            bessi, air_temperature, 0.0, 86_400.0, surface_temperature,
            true, 0.0, false, 0.0, false, 0.0, true, relative_humidity,
            air_pressure, wind_speed, z0m_snow, bessi.ϵ_snow,
        )
        @test bessi_sensible ≈ bessi.turbulent_flux_scheme.D_sh * (air_temperature - surface_temperature)
        @test bessi_latent ≈ Chion._bessi_latent_vapor_flux(
            surface_temperature, bessi, air_temperature, relative_humidity, air_pressure,
        )
        @test sensible != bessi_sensible
        @test latent != bessi_latent

        _, prescribed_sensible, prescribed_latent, _ =
            Chion._resolved_nonshortwave_surface_flux_components(
                c, air_temperature, 0.0, 86_400.0, surface_temperature,
                true, 0.0, true, 12.0, true, -8.0, true, relative_humidity,
                air_pressure, wind_speed, z0m_snow, c.ϵ_snow,
            )
        @test prescribed_sensible == 12.0
        @test prescribed_latent == -8.0
    end

    # The old models set `D_sh = 0`, which only acts in the BESSI turbulent-flux scheme; the
    # sensible and latent fluxes are prescribed (zero) in these solves, so it is not needed.
    @testset "Energy-flux two-layer diffusion uses physical interface conductance" begin
        model = _bessi_model(Grid(1; layering=Layering(Ntot=2, ice_substrate_layers=0));
            params = Params(ϵ_snow=0.0))
        state = model.state
        state.N[1] = 2
        state.mass[:, 1] .= 300.0
        state.density[:, 1] .= 300.0
        state.temperature[:, 1] .= (270.0, 260.0)
        state.Tsrf[1] = 270.0
        dt_seconds = 3600.0

        Chion.go_energy_flux!(
            model, 1, 270.0, 0.0, 0.0, 0.0, dt_seconds;
            q_sw_net=0.0, q_lw_down=0.0, q_sh=0.0, q_lh=0.0,
        )

        k1 = Chion._snow_thermal_conductivity(300.0, 270.0, model.params.rho_i)
        k2 = Chion._snow_thermal_conductivity(300.0, 260.0, model.params.rho_i)
        conductance = Chion.interface_conductance(k1, 1.0, k2, 1.0)
        a = dt_seconds * conductance / (model.params.ci * 300.0)
        temperature_difference = 270.0 - 260.0
        expected = (
            270.0 - a * temperature_difference / (1 + 2a),
            260.0 + a * temperature_difference / (1 + 2a),
        )
        @test state.temperature[:, 1] ≈ collect(expected) atol=1e-12
    end

    @testset "Robin surface boundary melts over a subfreezing top-cell centre" begin
        model = _bessi_model(Grid(1; layering=Layering(Ntot=1, ice_substrate_layers=0));
            params = Params(ϵ_snow=0.0))
        state = model.state
        state.N[1] = 1
        state.mass[1, 1] = 100.0
        state.density[1, 1] = 300.0
        state.temperature[1, 1] = 270.0
        state.Tsrf[1] = 270.0

        energy = Chion.go_energy_flux!(
            model, 1, 280.0, 400.0, 0.0, 0.0, 86_400.0;
            q_lw_down=0.0, q_sh=0.0, q_lh=0.0,
        )

        @test state.Tsrf[1] == model.params.T0
        @test state.temperature[1, 1] < model.params.T0
        @test energy.melt_energy_available > 0.0
    end

    @testset "Robin melting boundary retains subsurface conduction" begin
        model = _bessi_model(Grid(1; layering=Layering(Ntot=2, ice_substrate_layers=0));
            params = Params(ϵ_snow=0.0))
        state = model.state
        state.N[1] = 2
        state.mass[:, 1] .= (100.0, 300.0)
        state.density[:, 1] .= (300.0, 500.0)
        state.temperature[:, 1] .= (270.0, 260.0)
        state.Tsrf[1] = 270.0
        before = copy(state.temperature[:, 1])
        input_flux = 100.0
        dt_seconds = 86_400.0
        energy = Chion.go_energy_flux!(
            model, 1, 280.0, 0.0, 0.0, 0.0, dt_seconds;
            q_sw_net=input_flux, q_lw_down=0.0, q_sh=0.0, q_lh=0.0,
        )

        sensible_change = sum(state.mass[:, 1] .* model.params.ci .* (state.temperature[:, 1] .- before))
        @test state.Tsrf[1] == model.params.T0
        @test state.temperature[1, 1] < model.params.T0
        @test sensible_change + energy.melt_energy_available ≈ input_flux * dt_seconds atol=1e-6
    end

    @testset "Fresh snow on bare ice takes the air temperature in every new layer" begin
        model = _bessi_model(Grid(1); params = Params(diurnal_cycle=NoDiurnalCycle()))
        state = model.state
        state.ice_temperature .= 250.0
        state.Tsrf[1] = 250.0
        # 30 mm w.e. of snowfall is split over several fine near-surface layers.
        forcing = Chion.SnowpackStepForcing(250.0, 1.0, 30.0 / 86_400, 0.0, 0.0, 5.0;
            q_sw_net=0.0, has_q_sw_net=true, q_lw_down=0.0, has_q_lw_down=true,
            q_sh=0.0, has_q_sh=true, q_lh=0.0, has_q_lh=true)
        Chion.column_step!(Chion.get_fields(state), model.params, model.grid.layering, model.boundary, 1, forcing,
            model.cache.workspace)
        @test state.N[1] >= 2
        # Only radiative cooling acts, so no layer may end up warmer than the air.
        @test maximum(state.temperature[1:state.N[1], 1]) <= 250.0 + 1e-9
    end

    @testset "Trace snowfall on bare ice cannot store rain unseen" begin
        model = _bessi_model(Grid(1); params = Params(diurnal_cycle=NoDiurnalCycle()))
        state = model.state
        fields = Chion.get_fields(state)
        workspace = model.cache.workspace
        storage() = sum(state.mass[k, 1] + state.mass_w[k, 1] for k in 1:state.N[1]; init=0.0)
        fluxes(sf, rf) = (; q_sw_net=0.0, has_q_sw_net=true, q_lw_down=300.0, has_q_lw_down=true,
            q_sh=0.0, has_q_sh=true, q_lh=0.0, has_q_lh=true)
        # A trace of snowfall (far below the snow threshold), then a rainy day.
        for (snow, rain) in ((1e-24, 0.0), (0.0, 20.0 / 86_400))
            before = (S=storage(), smb=state.smb_ice[1], R=state.runoff[1], V=state.vapor_mass[1])
            forcing = Chion.SnowpackStepForcing(277.15, 1.0, snow, rain, 0.0, 5.0; fluxes(snow, rain)...)
            Chion.column_step!(fields, model.params, model.grid.layering, model.boundary, 1, forcing, workspace)
            precipitation = (snow + rain) * 86_400
            budget = (storage() - before.S) + (state.smb_ice[1] - before.smb) -
                     (precipitation - (state.runoff[1] - before.R) + (state.vapor_mass[1] - before.V))
            @test abs(budget) < 1e-9
            @test state.N[1] == 0
        end
    end

    @testset "Public energy solve couples the ice substrate on bare ice" begin
        model = _bessi_model(Grid(1))
        state = model.state
        state.ice_temperature .= model.params.T0
        state.Tsrf[1] = model.params.T0
        energy = Chion.go_energy_flux!(model, 1, 250.0, 0.0, 0.0, 0.0, 86_400.0;
            q_sw_net=0.0, q_lw_down=0.0, q_sh=0.0, q_lh=0.0)
        @test !energy.needs_melt
        @test state.Tsrf[1] < model.params.T0
        @test state.ice_temperature[1, 1] < model.params.T0
        @test_throws ErrorException Chion.go_energy_flux!(model, 1, 250.0, 0.0, 0.0, 0.0, 86_400.0;
            scratch=Chion.EnergyWorkspace(state.mass, Float64, 1, 1), q_sw_net=0.0, q_lw_down=0.0, q_sh=0.0, q_lh=0.0)
    end

    @testset "BESSIModel scheme types" begin
        grid = Grid(1)
        m1 = _bessi_model(grid; params = Params(albedo_scheme=DynamicAlbedo()))
        m2 = _bessi_model(grid; params = Params(albedo_scheme=ConstantAlbedo()))
        m3 = _bessi_model(grid; params = Params(albedo_scheme=AgingAlbedo()))
        m4 = _bessi_model(grid; params = Params(densification_scheme=HTESSELDensification()))
        m5 = _bessi_model(grid; params = Params(fresh_snow_density_scheme=ParameterizedFreshSnowDensity()))
        m6 = _bessi_model(grid)
        m7 = _bessi_model(grid; params = Params(turbulent_flux_scheme=BESSITurbulentFluxes()))
        @test m1.params.albedo_scheme isa DynamicAlbedo
        @test m2.params.albedo_scheme isa ConstantAlbedo
        @test m3.params.albedo_scheme isa AgingAlbedo
        @test m3.params.alpha_dry == 0.81
        @test m3.params.alpha_wet == 0.70
        @test m3.params.albedo_scheme.cold_timescale_days == 20.0
        @test m3.params.albedo_scheme.melting_timescale_days == 2.0
        @test m4.params.densification_scheme isa HTESSELDensification
        @test m1.params.densification_scheme isa BESSIDensification
        @test m5.params.fresh_snow_density_scheme isa ParameterizedFreshSnowDensity
        @test m1.params.fresh_snow_density_scheme isa ConstantFreshSnowDensity
        # Scheme names are accepted only in TOML configurations.
        @test Chion._scheme_from_config(Chion.AbstractAlbedoScheme, "semix") isa SemixAlbedo
        @test Chion._scheme_from_config(Chion.AbstractDensificationScheme, "htessel") isa HTESSELDensification
        @test typeof(m2.params) == typeof(Params(albedo_scheme=ConstantAlbedo()))
        @test m6.params.turbulent_flux_scheme === SemixTurbulentFluxes()
        @test m6.params.surface_energy_balance_scheme === SemixSurfaceEnergyBalance()
        @test m7.params.turbulent_flux_scheme === BESSITurbulentFluxes()

        for (kind, name) in (
            (Chion.AbstractFreshSnowDensityScheme, "bessi"),
            (Chion.AbstractFreshSnowDensityScheme, "htessel"),
            (Chion.AbstractAlbedoScheme, "legacy"),
            (Chion.AbstractAlbedoScheme, "bessi"),
        )
            alias_error = _captured_exception(() -> Chion._scheme_from_config(kind, name))
            @test alias_error isa Exception
            @test occursin("Use \"", sprint(showerror, alias_error))
        end
        @test Params(diurnal_cycle=DiurnalSubsteps()).diurnal_cycle isa DiurnalSubsteps
        @test_throws ArgumentError Params(diurnal_shortwave=true)
        @test_throws MethodError BESSIModel(grid; forcing=_placeholder_forcing(), diurnal_shortwave=true)
        @test PDDParams(method=PISMPDD()).method isa PISMPDD
        @test PDDParams().method isa SimplePDD
        @test PDDModel(grid; forcing=_placeholder_forcing(), params=PDDParams(method=PISMPDD())).params.method isa PISMPDD
        @test Chion._scheme_from_config(Chion.AbstractPDDMethod, "pism") isa PISMPDD
        @test_throws Exception PDDParams(method=:calov_greve)
        @test_throws Exception Chion._scheme_from_config(Chion.AbstractPDDMethod, "calov_greve")
        @test Params(
            albedo_scheme=SemixAlbedo(snow_bands=WarrenWiscombeSnowAlbedo()),
        ).albedo_scheme.snow_bands === WarrenWiscombeSnowAlbedo()

        accumulation_state = m6.state
        @test_throws MethodError Chion._apply_accumulation!(
            accumulation_state.N,
            accumulation_state.mass,
            accumulation_state.mass_w,
            accumulation_state.density,
            accumulation_state.temperature,
            accumulation_state.mass_base,
            accumulation_state.smb_ice,
            accumulation_state.runoff,
            accumulation_state.Tsrf,
            accumulation_state.albedo,
            1,
            m6.params,
            m6.grid.layering.Ntot,
            m6.grid.layering.mass_max,
            m6.grid.layering.mass_split,
            m6.grid.layering.mass_min,
            0.0,
            0.0,
            m6.params.seconds_per_day;
            T_air=260.0,
        )

        err = _captured_exception() do
            Params(albedo_scheme=AgingAlbedo(), alpha_dry=0.5, alpha_wet=0.6)
        end
        @test err isa Exception
        @test occursin("alpha_wet", sprint(showerror, err))
    end

    @testset "aging albedo follows snowfall age" begin
        model = _bessi_model(Grid(1); params = Params(albedo_scheme=AgingAlbedo()))
        state = model.state
        params = model.params
        state.N[1] = 1
        state.mass[1, 1] = 300.0
        state.temperature[1, 1] = params.T0 - 1.0

        @test state.albedo[1] == params.alpha_dry
        @test state.snow_age_days[1] == 0.0

        Chion.update_surface_albedo!(model, 1, 1.0)
        @test state.snow_age_days[1] == 1.0
        @test state.albedo[1] ≈ params.alpha_wet +
                                 (params.alpha_dry - params.alpha_wet) *
                                 exp(-1 / params.albedo_scheme.cold_timescale_days) atol=1e-12

        state.temperature[1, 1] = params.T0
        Chion.update_surface_albedo!(model, 1, 1.0)
        @test state.snow_age_days[1] == 2.0
        @test state.albedo[1] ≈ params.alpha_wet +
                                 (params.alpha_dry - params.alpha_wet) *
                                 exp(
                                     -1 / params.albedo_scheme.cold_timescale_days -
                                     1 / params.albedo_scheme.melting_timescale_days,
                                 ) atol=1e-12

        Chion._update_aging_surface_albedo_arrays!(
            state.N,
            state.mass,
            state.temperature,
            state.albedo,
            state.snow_age_days,
            1,
            params,
            1.0 / params.seconds_per_day,
            1.0,
        )
        @test state.snow_age_days[1] == 0.0
        @test state.albedo[1] == params.alpha_dry

        state.N[1] = 0
        Chion.update_surface_albedo!(model, 1, 1.0)
        @test state.snow_age_days[1] == 0.0
        @test state.albedo[1] == params.alpha_ice
    end

    @testset "SnowpackForcing conversions and validation" begin
        noleap_times = Chion._synthesized_time_values(fill(1.0, 365))
        @test first(noleap_times) == DateTime(2001, 1, 1, 12)
        @test last(noleap_times) == DateTime(2001, 12, 31, 12)
        @test all(value -> !(month(value) == 2 && day(value) == 29), noleap_times)

        forcing = SnowpackForcing(
            dt_days           = [1.0, 2.0],
            ncol              = 2,
            air_temperature_c = [-10.0, -5.0],
            snowfall_mm_day   = 1.5,
            rainfall_mm_day   = [0.0, 2.0],
            shortwave_down    = [100.0 120.0; 140.0 160.0],
            wind_speed        = 3.5,
            q_lw_down         = [250.0, 255.0],
            has_q_lw_down     = [true, false],
        )
        @test size(forcing.air_temperature) == (2, 2)
        @test forcing.air_temperature[1, 1] ≈ 263.15 atol=1e-8
        @test forcing.air_temperature[2, 2] ≈ 268.15 atol=1e-8
        @test forcing.snowfall_rate[2, 1]   ≈ 1.5 / 86_400.0 atol=1e-12
        @test forcing.rainfall_rate[1, 2]   ≈ 2.0 / 86_400.0 atol=1e-12
        @test forcing.shortwave_down == [100.0 120.0; 140.0 160.0]
        @test forcing.wind_speed == fill(3.5, 2, 2)
        @test forcing.has_q_lw_down == [true false; true false]
        @test length(forcing.time_values) == 2

        err = _captured_exception() do
            SnowpackForcing(
                dt_days=[1.0, 1.0], ncol=2,
                air_temperature_c=[-10.0],
                snowfall_mm_day=0.0, rainfall_mm_day=0.0, shortwave_down=0.0,
            )
        end
        @test err isa Exception
        @test occursin("air_temperature_c", sprint(showerror, err))
    end

    @testset "SnowpackForcing keeps invariant fields compact" begin
        forcing = SnowpackForcing(
            dt_days=[1.0, 1.0],
            ncol=3,
            air_temperature_c=[-10.0, -9.0],
            snowfall_mm_day=1.0,
            rainfall_mm_day=0.0,
            shortwave_down=100.0,
            latitude_deg=[60.0, 70.0, 80.0],
            surface_height=[0.0, 500.0, 1000.0],
        )

        @test forcing.q_lw_down isa Chion.ConstantForcingMatrix
        @test forcing.has_q_lw_down isa Chion.ConstantForcingMatrix
        @test forcing.air_temperature isa Chion.TimeForcingMatrix
        @test forcing.shortwave_down isa Chion.ConstantForcingMatrix
        @test forcing.latitude_deg isa Chion.ColumnForcingMatrix
        @test forcing.surface_height isa Chion.ColumnForcingMatrix
        @test Matrix(forcing.q_lw_down) == zeros(3, 2)
        @test Matrix(forcing.latitude_deg) == [60.0 60.0; 70.0 70.0; 80.0 80.0]
    end

    @testset "missing column fields and annual-mean pressure are stored per column" begin
        dates = DateTime(2001, 1, 1) .+ Day.(0:3)
        T = [250.0 252.0 254.0 256.0; 260.0 262.0 264.0 266.0]
        bare = SnowpackForcing(dt_days=ones(4), air_temperature=T, snowfall_rate=0.0, rainfall_rate=0.0,
            shortwave_down=0.0, time_values=dates)
        for name in (:latitude_deg, :surface_height, :ice_thickness, :annual_pdd)
            field = getproperty(bare, name)
            @test field isa Chion.ColumnForcingMatrix
            @test size(field) == (2, 4) && all(isnan, field)
        end
        # written per column, as Kryonomos and update_air_pressure! do
        bare.surface_height .= [100.0, 200.0]
        @test Matrix(bare.surface_height) == [fill(100.0, 1, 4); fill(200.0, 1, 4)]

        height = [0.0, 1500.0]
        annual = SnowpackForcing(dt_days=ones(4), air_temperature=T, snowfall_rate=0.0, rainfall_rate=0.0,
            shortwave_down=0.0, time_values=dates, surface_height=height, air_pressure_temperature_mode=:annual_mean)
        dense = Chion.air_pressure_from_surface_height(height, T; dt_days=ones(4), time_values=dates,
            temperature_mode=:annual_mean)
        @test annual.air_pressure isa Chion.ColumnForcingMatrix
        @test Matrix(annual.air_pressure) == dense                       # bitwise
        update_air_pressure!(annual; temperature_mode=:annual_mean)
        @test Matrix(annual.air_pressure) == dense
        # a time-varying write into a per-column field is an error, not a silent truncation
        err = _captured_exception(() -> update_air_pressure!(annual; temperature_mode=:instantaneous))
        @test err isa ErrorException && occursin("one value per column", sprint(showerror, err))
        @test Matrix(annual.air_pressure) == dense

        # instantaneous pressure and forcing spanning two years keep the full matrix
        inst = SnowpackForcing(dt_days=ones(4), air_temperature=T, snowfall_rate=0.0, rainfall_rate=0.0,
            shortwave_down=0.0, time_values=dates, surface_height=height)
        @test inst.air_pressure isa Matrix
        two_years = [DateTime(2001, 12, 30), DateTime(2001, 12, 31), DateTime(2002, 1, 1), DateTime(2002, 1, 2)]
        split = SnowpackForcing(dt_days=ones(4), air_temperature=T, snowfall_rate=0.0, rainfall_rate=0.0,
            shortwave_down=0.0, time_values=two_years, surface_height=height, air_pressure_temperature_mode=:annual_mean)
        @test split.air_pressure isa Matrix
        @test split.air_pressure == Chion.air_pressure_from_surface_height(height, T; dt_days=ones(4),
            time_values=two_years, temperature_mode=:annual_mean)

        # rebuilding and slicing keep the compact fields compact
        rebuilt = SnowpackForcing(dt_days=annual.dt_days, air_temperature=annual.air_temperature, snowfall_rate=0.0,
            rainfall_rate=0.0, shortwave_down=0.0, time_values=annual.time_values,
            surface_height=annual.surface_height, air_pressure=annual.air_pressure)
        @test rebuilt.air_pressure isa Chion.ColumnForcingMatrix && rebuilt.surface_height isa Chion.ColumnForcingMatrix
        sliced = Chion.slice_forcing(annual, 2:3)
        @test sliced.air_pressure isa Chion.ColumnForcingMatrix && size(sliced.air_pressure) == (2, 2)
        @test Matrix(sliced.air_pressure) == dense[:, 2:3]
    end

    @testset "transposed layer storage preserves logical indexing" begin
        logical = reshape(collect(1.0:12.0), 3, 4)
        storage = Chion.TransposedLayerMatrix(permutedims(logical, (2, 1)))
        @test size(storage) == size(logical)
        @test Array(storage) == logical
        storage[2, 3] = -1.0
        @test storage.parent[3, 2] == -1.0

        scratch = similar(storage, Float64, 3, 4)
        @test scratch isa Chion.TransposedLayerMatrix
        @test size(scratch.parent) == (4, 3)
    end

    # State output is an `Output` with one snapshot record per forcing step.
    @testset "Simulation saves an exact state field" begin
        mktempdir() do dir
            model, forcing, grid = _sample_model_forcing_grid()
            path = joinpath(dir, "save_exact.nc")
            sim = Simulation(model; verbose=false,
                output=[Output(path; variables=[:thickness], frequency=:step, kind=:snapshot)])
            run!(sim)
            @test finished(sim)
            @test isfile(path)
            ds = NCDataset(path)
            @test ds.attrib["status"] == "complete"
            @test haskey(ds, "thickness")
            @test haskey(ds, "time")
            @test !haskey(ds, "cycle")
            @test !haskey(ds, "month_cycle")
            @test !haskey(ds, "step_cycle")
            @test !haskey(ds, "history_mean_thickness")
            @test size(ds["thickness"]) == (2, 2, 3)
            @test ds["time"].var[:] ≈ Chion._netcdf_time_days.(forcing.time_values)
            @test ds["time"].attrib["standard_name"] == "time"
            @test ds["time"].attrib["calendar"] == "proleptic_gregorian"
            @test ds["time"].attrib["axis"] == "T"
            @test "time" in NCDatasets.unlimited(ds)
            @test ds.dim["time"] == 3
            @test !haskey(ds.attrib, "cycles_completed")
            # the last snapshot is the final state
            update_diagnostics!(model)
            @test [ds["thickness"][grid.is[c], grid.js[c], 3] for c in 1:grid.ncol] ≈
                  Float32.(model.cache.thickness)
            close(ds)
        end
    end

    @testset "Simulation saves all default state fields" begin
        mktempdir() do dir
            model, forcing, grid = _sample_model_forcing_grid()
            path = joinpath(dir, "save_group.nc")
            sim = Simulation(model; verbose=false,
                output=[Output(path; variables=available_output_variables(model), frequency=:step, kind=:snapshot)])
            run!(sim)
            @test finished(sim)
            ds = NCDataset(path)
            @test ds.attrib["status"] == "complete"
            @test haskey(ds, "thickness")
            @test haskey(ds, "mass_base")
            @test !haskey(ds, "history_mean_thickness")
            close(ds)
        end
    end

    @testset "Simulation saves mixed fields" begin
        mktempdir() do dir
            model, forcing, grid = _sample_model_forcing_grid()
            path = joinpath(dir, "save_mixed.nc")
            sim = Simulation(model; verbose=false,
                output=[Output(path; variables=[:thickness, :mass_base], frequency=:step, kind=:snapshot)])
            run!(sim)
            @test finished(sim)
            ds = NCDataset(path)
            @test ds.attrib["status"] == "complete"
            @test haskey(ds, "thickness")
            @test haskey(ds, "mass_base")
            @test !haskey(ds, "history_mean_thickness")
            close(ds)
        end
    end

    @testset "Direct monthly variable selection records values" begin
        mktempdir() do dir
            model, forcing, _ = _sample_model_forcing_grid()
            path = joinpath(dir, "monthly_lhf.nc")
            sim = run!(Simulation(model; verbose=false,
                output=[Output(path; variables=[:latent_heat_flux], frequency=:monthly)]))
            @test finished(sim)
            ds = NCDataset(path)
            @test ds.dim["time"] == 3
            @test ds["time"].var[:] ≈ Chion._netcdf_time_days.(forcing.time_values)
            @test all(isfinite, Float64.(ds["latent_heat_flux"][:, :, 1:3]))
            close(ds)
        end
    end

    @testset "Simulation can skip NetCDF entirely" begin
        model, forcing, _ = _sample_model_forcing_grid()
        sim = Simulation(model; verbose=false)
        @test run!(sim) === sim
        @test finished(sim)
        @test isempty(sim.output.outputs)
        @test isempty(sim.io.recorders)
    end

    @testset "Bare-ice rainfall is routed directly to runoff" begin
        params = Params(diurnal_cycle=NoDiurnalCycle())
        rainfall_mm = 7.0
        neutral_longwave_down = params.σ * params.ϵ_snow * params.T0^4
        forcing = SnowpackForcing(
            dt_days=[1.0],
            air_temperature=[params.T0],
            snowfall_rate=[0.0],
            rainfall_rate=[rainfall_mm / params.seconds_per_day],
            shortwave_down=[0.0],
            q_lw_down=[neutral_longwave_down],
            q_sh=[0.0],
            q_lh=[0.0],
        )
        model = BESSIModel(Grid(1; layering=Layering(Ntot=4)); forcing, params)
        sim = Simulation(model; verbose=false)

        run!(sim)

        @test finished(sim)
        state = model.state
        @test state.runoff[1] ≈ rainfall_mm atol=1e-12
        @test state.smb_ice[1] ≈ 0.0 atol=1e-12
        @test state.melt[1] ≈ 0.0 atol=1e-12
        @test sum(state.mass) ≈ 0.0 atol=1e-12
        @test sum(state.mass_w) ≈ 0.0 atol=1e-12
    end

    @testset "Prescribed latent heat controls snow vapor mass" begin
        params = Params(albedo_scheme=ConstantAlbedo(), diurnal_cycle=NoDiurnalCycle())
        q_lh = 20.0
        forcing = SnowpackForcing(
            dt_days=[1.0],
            air_temperature=[263.15],
            snowfall_rate=[0.0],
            rainfall_rate=[0.0],
            shortwave_down=[0.0],
            q_sw_net=[0.0],
            q_lw_down=[0.0],
            q_sh=[0.0],
            q_lh=[q_lh],
            relative_humidity=[0.0],
            air_pressure=[101325.0],
        )
        model = BESSIModel(Grid(1; layering=Layering(Ntot=1)); forcing, params)
        state = model.state
        state.N[1] = 1
        state.mass[1, 1] = 100.0
        state.density[1, 1] = 400.0
        state.temperature[1, 1] = 263.15
        state.Tsrf[1] = 263.15
        sim = Simulation(model; verbose=false)
        @test finished(run!(sim))
        expected = q_lh * params.seconds_per_day / (params.Lv + params.Lm)

        @test state.vapor_mass[1] ≈ expected rtol=1e-12
        @test state.mass[1, 1] ≈ 100.0 + expected rtol=1e-12
    end

    @testset "BESSI KA interval step matches serial column stepping" begin
        serial, forcing, grid = _sample_model_forcing_grid()
        stepped = BESSIModel(grid; forcing)
        serial_fields = Chion.get_fields(serial.state)

        for time_index in eachindex(forcing.time_values)
            @inbounds for idx in 1:grid.ncol
                Chion.column_step!(
                    serial_fields,
                    serial.params,
                    serial.grid.layering,
                    serial.boundary,
                    idx,
                    Chion._step_forcing_at(serial.forcing.fields, idx, Int(time_index)),
                    serial.cache.workspace,
                )
            end
        end

        integrate!(stepped, 1:length(forcing.time_values))

        @test stepped.state.N == serial.state.N
        @test stepped.state.mass ≈ serial.state.mass
        @test stepped.state.mass_w ≈ serial.state.mass_w
        @test stepped.state.density ≈ serial.state.density
        @test stepped.state.temperature ≈ serial.state.temperature
        @test stepped.state.smb_ice ≈ serial.state.smb_ice
        @test stepped.state.runoff ≈ serial.state.runoff
        @test stepped.state.Tsrf ≈ serial.state.Tsrf
        @test stepped.state.albedo ≈ serial.state.albedo
        @test stepped.state.snow_age_days ≈ serial.state.snow_age_days
    end

    # The simulation steps the model's state in place; `cpu_state` gives an independent copy.
    @testset "Simulation owns reference and current state" begin
        model, forcing, _ = _sample_model_forcing_grid()
        sim = Simulation(model; verbose=false)
        ref = cpu_state(model)
        @test ref !== model.state
        @test ref.mass !== model.state.mass
        @test all(ref.mass .== 0.0)
        run!(sim)
        @test all(ref.mass .== 0.0)
        @test sum(model.state.mass) > 0.0
    end

    @testset "time_step! loop matches run wrapper" begin
        model_a, forcing, grid = _sample_model_forcing_grid()
        sim_run = Simulation(model_a; verbose=false)
        run!(sim_run)

        model_b = BESSIModel(grid; forcing)
        sim_manual = Simulation(model_b; verbose=false)
        steps = 0
        while !finished(sim_manual)
            time_step!(sim_manual)
            steps += 1
        end

        @test finished(sim_run) && finished(sim_manual)
        @test steps == length(forcing.time_values)
        @test sim_run.clock.iteration == sim_manual.clock.iteration == length(forcing.time_values)
        @test model_a.state.mass ≈ model_b.state.mass
        @test model_a.state.smb_ice ≈ model_b.state.smb_ice
    end

    @testset "yearly_step! returns weighted annual coupling fields" begin
        expected_model, forcing, grid = _sample_model_forcing_grid()
        expected_sim = Simulation(expected_model; verbose=false)
        expected_smb_before = copy(expected_model.state.smb_ice)
        expected_temperature_sum = zeros(Float64, grid.ncol)
        for time_index in eachindex(forcing.dt_days)
            time_step!(expected_sim)
            @. expected_temperature_sum += expected_model.state.Tsrf * forcing.dt_days[time_index]
        end

        coupled_sim = Simulation(BESSIModel(grid; forcing); verbose=false)
        coupling = yearly_step!(coupled_sim)

        @test coupling.year == 2001
        @test coupling.mean_T_srf_K ≈ expected_temperature_sum ./ sum(forcing.dt_days)
        @test coupling.smb_ice_yearly ≈ expected_model.state.smb_ice .- expected_smb_before
        expected_grid = fill(NaN, size(grid.mask))
        for column in 1:grid.ncol
            expected_grid[grid.js[column], grid.is[column]] = coupling.mean_T_srf_K[column]
        end
        @test coupling.mean_T_srf_K_grid ≈ expected_grid
    end

    @testset "external forcing step matches scheduled forcing" begin
        grid = Grid(1; layering=Layering(Ntot=4))
        params = Params(diurnal_cycle=NoDiurnalCycle())
        forcing = SnowpackForcing(
            dt_days=[1.0],
            air_temperature_c=[-8.0],
            snowfall_mm_day=[2.0],
            rainfall_mm_day=[0.0],
            shortwave_down=[100.0],
        )
        scheduled = BESSIModel(grid; forcing, params)
        run!(Simulation(scheduled; verbose=false))

        external = BESSIModel(grid; forcing, params)
        sim = Simulation(external; verbose=false)
        time_step!(sim)

        @test finished(sim)
        @test external.state.mass ≈ scheduled.state.mass
        @test external.state.smb_ice ≈ scheduled.state.smb_ice
    end

    @testset "model exposes mutable scheduled forcing" begin
        grid = Grid(2; layering=Layering(Ntot=4))
        forcing = SnowpackForcing(
            dt_days=[1.0, 1.0],
            air_temperature_c=[-8.0 -7.0; -10.0 -9.0],
            snowfall_mm_day=1.0,
            rainfall_mm_day=0.0,
            shortwave_down=100.0,
        )
        model = BESSIModel(grid; forcing, params=Params(diurnal_cycle=NoDiurnalCycle()))
        surface_height = [0.0, 1000.0]

        set_input!(model, :surface_height, surface_height)
        expected = Chion.air_pressure_from_surface_height(
            surface_height,
            forcing.air_temperature;
            dt_days=forcing.dt_days,
            time_values=forcing.time_values,
            temperature_mode=:instantaneous,
        )
        # `set_input!` does not recompute the air pressure: pass it as well.
        set_input!(model, :air_pressure, expected)
        @test Matrix(model.forcing.data.surface_height) == [surface_height surface_height]
        @test model.forcing.data.air_pressure ≈ expected
        @test model.forcing.fields.air_pressure ≈ expected
    end

    # Grids without coordinates write the columns along a `column` dimension.
    @testset "non-spatial grids write NetCDF output per column" begin
        grid = Grid(1; layering=Layering(Ntot=4))
        params = Params(diurnal_cycle=NoDiurnalCycle())
        forcing = SnowpackForcing(
            dt_days=[1.0, 1.0],
            air_temperature_c=[-10.0, -9.0],
            snowfall_mm_day=0.5,
            rainfall_mm_day=0.0,
            shortwave_down=[100.0, 120.0],
        )
        sim = run!(Simulation(BESSIModel(grid; forcing, params); verbose=false))
        @test finished(sim)
        @test isempty(sim.output.outputs)

        mktempdir() do dir
            path = joinpath(dir, "column.nc")
            run!(Simulation(BESSIModel(grid; forcing, params); verbose=false,
                output=[Output(path; variables=[:thickness], frequency=:step, kind=:snapshot)]))
            ds = NCDataset(path)
            @test !haskey(ds.dim, "x")
            @test ds.dim["column"] == 1
            @test size(ds["thickness"]) == (1, 2)
            close(ds)
        end
    end

    @testset "forcing-file loading and run!" begin
        mktempdir() do dir
            forcing_path = _write_sample_forcing_file(joinpath(dir, "forcing_sample.nc"))
            loaded = load_forcing_file(forcing_path)
            @test size(loaded.forcing.air_temperature) == (4, 2)
            @test size(loaded.forcing.snowfall_rate) == (4, 2)
            @test loaded.forcing.time_values ==
                  [DateTime(2001, 1, 1, 12), DateTime(2001, 1, 2, 12)]
            @test loaded.grid isa Grid
            @test sort(collect(zip(loaded.grid.js, loaded.grid.is))) == [(1, 1), (1, 2), (2, 1), (2, 2)]

            masked = load_forcing_file(forcing_path; mask_name="MSK", mask_threshold=50.0)
            @test size(masked.forcing.air_temperature) == (2, 2)
            @test collect(zip(masked.grid.js, masked.grid.is)) == [(1, 1), (1, 2)]
            @test masked.grid.mask == [100.0 75.0; 0.0 25.0]

            model = BESSIModel(Grid(loaded.grid; layering=Layering(Ntot=4));
                forcing=Forcing(loaded.forcing; grid=loaded.grid),
                params=Params(diurnal_cycle=NoDiurnalCycle()))
            sim = run!(Simulation(model; verbose=false))
            @test finished(sim)
            @test isempty(sim.output.outputs)
        end
    end

    @testset "forcing-file loader accepts singleton-layer 4D variables" begin
        mktempdir() do dir
            forcing_path = _write_sample_forcing_file_4d(joinpath(dir, "forcing_sample_4d.nc"))
            loaded = load_forcing_file(forcing_path)
            @test size(loaded.forcing.air_temperature) == (4, 2)
            @test loaded.forcing.air_temperature[1, 1] ≈ 258.15 atol=1e-5
            @test loaded.forcing.air_temperature[1, 2] ≈ 259.15 atol=1e-5
            @test all(loaded.forcing.has_q_lw_down)
        end
    end

    @testset "PDDModel runs bulk positive degree days" begin
        grid = Grid(1)
        f = SnowpackForcing(
            dt_days=[1.0, 1.0],
            air_temperature_c=[-5.0, 1.0],
            snowfall_mm_day=[2.0, 0.0],
            rainfall_mm_day=0.0,
            shortwave_down=150.0,
        )
        params = PDDParams(ddf_snow=3.0, ddf_ice=8.0, refreezing_fraction=0.0)
        pdd = PDDModel(grid; forcing=f, params)
        pdd_sim = run!(Simulation(pdd; verbose=false))
        @test finished(pdd_sim)
        @test pdd.state.snowpack_swe[1] ≈ 0.0 atol=1e-12
        @test pdd.state.smb_ice[1] ≈ -(8.0 / 3.0) atol=1e-12
        @test pdd.state.runoff[1] ≈ 2.0 + 8.0 / 3.0 atol=1e-12
        @test pdd.state.pdd_sum[1] ≈ 1.0 atol=1e-12

        mktempdir() do dir
            spatial_grid = Grid(1;
                x=[0.0],
                y=[0.0],
                js=[1],
                is=[1],
                mask=ones(1, 1),
            )
            spatial_sim = run!(Simulation(PDDModel(spatial_grid; forcing=f, params); verbose=false))
            @test finished(spatial_sim)
            @test isempty(spatial_sim.output.outputs)
            @test spatial_sim.clock.year_index == 2

            # Two calendar years: the forcing year run again as 2002.
            two_years = Forcing(y -> relabel_year(f, y), 2001:2002; grid=spatial_grid)
            output_path = joinpath(dir, "pdd.nc")
            output_model = PDDModel(spatial_grid; forcing=two_years, params)
            output_sim = Simulation(output_model; verbose=false,
                output=[Output(output_path; variables=collect(Chion.PDD_OUTPUT_VARS), frequency=:step,
                    kind=:snapshot)])
            run!(output_sim)
            @test finished(output_sim)

            ds = NCDataset(output_path)
            @test ds.attrib["status"] == "complete"
            @test all(haskey(ds, name) for name in ("snowpack_swe", "smb_ice", "runoff", "pdd_sum"))
            @test !haskey(ds, "thickness")
            @test size(ds["smb_ice"]) == (1, 1, 4)
            expected_times = vcat(f.time_values, f.time_values .+ Year(1))
            @test ds["time"].var[:] ≈ Chion._netcdf_time_days.(expected_times)
            @test "time" in NCDatasets.unlimited(ds)
            @test ds.dim["time"] == 4
            # cumulative counters are written as totals per record
            @test sum(ds["smb_ice"][1, 1, :]) ≈ output_model.state.smb_ice[1] atol=1e-5
            @test sum(ds["pdd_sum"][1, 1, :]) ≈ output_model.state.pdd_sum[1] atol=1e-5
            @test ds["snowpack_swe"][1, 1, 4] ≈ output_model.state.snowpack_swe[1] atol=1e-5
            close(ds)

            monthly_path = joinpath(dir, "pdd_monthly.nc")
            monthly_model = PDDModel(spatial_grid; forcing=f, params)
            monthly_sim = Simulation(monthly_model; verbose=false,
                output=[Output(monthly_path; variables=[:smb_ice, :runoff, :pdd_sum], frequency=:monthly)])
            @test finished(run!(monthly_sim))
            ds = NCDataset(monthly_path)
            @test ds.dim["time"] == 1
            @test size(ds["smb_ice"]) == (1, 1, 1)
            # the record is stamped with the start of the period, its bounds span the steps
            @test ds["time"].var[:] ≈ Chion._netcdf_time_days.([f.time_values[1]])
            @test ds["time_bnds"].var[2, 1] ≈ Chion._netcdf_time_days(f.time_values[end])
            @test ds["smb_ice"][1, 1, 1] ≈ monthly_model.state.smb_ice[1] atol=1e-5
            @test ds["runoff"][1, 1, 1] ≈ monthly_model.state.runoff[1] atol=1e-5
            @test ds["pdd_sum"][1, 1, 1] ≈ monthly_model.state.pdd_sum[1] atol=1e-5
            close(ds)

            err = _captured_exception() do
                run!(Simulation(PDDModel(spatial_grid; forcing=f); verbose=false,
                    output=[Output(joinpath(dir, "unsupported.nc"); variables=[:thickness], frequency=:step)]))
            end
            @test err isa Exception
            @test occursin("is not an output variable of PDDModel", sprint(showerror, err))
        end
    end

    @testset "PDD uses the capped ice-facing Fortran budget" begin
        params = PDDParams(
            ddf_snow=3.0,
            ddf_ice=8.0,
            refreezing_fraction=0.6,
            H_snow_max=100.0,
            method=SimplePDD(),
        )

        cold_accumulation = SnowpackForcing(
            dt_days=[1.0],
            air_temperature=[params.T0 - 10.0],
            snowfall_rate=[150.0 / params.seconds_per_day],
            rainfall_rate=[20.0 / params.seconds_per_day],
            shortwave_down=[0.0],
        )
        accumulated = _pdd_after(params, cold_accumulation).state

        @test accumulated.snowpack_swe[1] ≈ 100.0 atol=1e-12
        @test accumulated.smb_ice[1] ≈ 50.0 atol=1e-12
        @test accumulated.runoff[1] ≈ 20.0 atol=1e-12
        @test accumulated.snowpack_swe[1] +
              accumulated.smb_ice[1] +
              accumulated.runoff[1] ≈ 170.0 atol=1e-12

        initial = initial_state(Grid(1), params)
        initial.snowpack_swe[1] = 10.0
        melt_forcing = SnowpackForcing(
            dt_days=[1.0],
            air_temperature=[params.T0 + 3.0],
            snowfall_rate=[0.0],
            rainfall_rate=[0.0],
            shortwave_down=[0.0],
        )
        capacity_limited_model = _pdd_after(params, melt_forcing; state=initial)
        capacity_limited = capacity_limited_model.state

        @test capacity_limited.snowpack_swe[1] ≈ 1.0 atol=1e-12
        @test capacity_limited.smb_ice[1] ≈ 0.6 atol=1e-12
        @test capacity_limited.runoff[1] ≈ 8.4 atol=1e-12
        @test capacity_limited.snowpack_swe[1] +
              capacity_limited.smb_ice[1] +
              capacity_limited.runoff[1] ≈ 10.0 atol=1e-12

        extreme_melt = SnowpackForcing(
            dt_days=[1.0],
            air_temperature=[params.T0 + 1000.0],
            snowfall_rate=[0.0],
            rainfall_rate=[0.0],
            shortwave_down=[0.0],
        )
        extreme = _pdd_after(params, extreme_melt; state=cpu_state(capacity_limited_model)).state
        @test extreme.snowpack_swe[1] == 0.0
        @test extreme.smb_ice[1] < 0.0
    end

    @testset "PDD method and physical constants are explicit" begin
        cold = SnowpackForcing(
            dt_days=[1.0],
            air_temperature_c=[-5.0],
            snowfall_mm_day=[0.0],
            rainfall_mm_day=[0.0],
            shortwave_down=[0.0],
        )
        simple_state = _pdd_after(PDDParams(method=SimplePDD()), cold).state
        pism_state = _pdd_after(PDDParams(method=PISMPDD()), cold).state

        @test simple_state.pdd_sum[1] == 0.0
        @test pism_state.pdd_sum[1] ≈ 0.4165773529384319 atol=1e-12

        custom_constants = PDDParams(
            method=SimplePDD(),
            T0=270.0,
            seconds_per_day=100.0,
        )
        custom_forcing = SnowpackForcing(
            dt_days=[1.0],
            air_temperature=[271.0],
            snowfall_rate=[1.0],
            rainfall_rate=[0.0],
            shortwave_down=[0.0],
        )
        custom_state = _pdd_after(custom_constants, custom_forcing).state
        @test custom_state.pdd_sum[1] ≈ 1.0 atol=1e-12
        @test custom_state.snowpack_swe[1] ≈ 97.0 atol=1e-12
        @test custom_state.smb_ice[1] ≈ 3.0 atol=1e-12

        @test _captured_exception(() -> PDDParams(method=:unknown)) isa Exception
        @test _captured_exception(() -> PDDParams(H_snow_max=0.0)) isa Exception
    end

    @testset "PDD honors active columns through the runtime path" begin
        grid = Grid(2)
        forcing = SnowpackForcing(
            dt_days=[1.0],
            ncol=2,
            air_temperature_c=-10.0,
            snowfall_mm_day=5.0,
            rainfall_mm_day=0.0,
            shortwave_down=0.0,
        )
        model = PDDModel(grid; forcing)
        sim = Simulation(model; verbose=false)
        set_active_mask!(sim, [true, false])
        time_step!(sim)

        @test model.state.snowpack_swe ≈ [5.0, 0.0] atol=1e-12
        @test model.state.smb_ice == [0.0, 0.0]
        @test model.state.runoff == [0.0, 0.0]
        @test model.state.pdd_sum == [0.0, 0.0]
    end

    if Chion.cuda_available()
        @testset "PDD GPU path matches CPU" begin
            grid = Grid(2)
            params = PDDParams(H_snow_max=8.0, method=PISMPDD())
            forcing = SnowpackForcing(
                dt_days=[1.0, 1.0],
                ncol=2,
                air_temperature_c=[-5.0, 2.0],
                snowfall_mm_day=5.0,
                rainfall_mm_day=1.0,
                shortwave_down=0.0,
            )
            cpu_model = PDDModel(grid; forcing, params)
            gpu_model = to_backend(PDDModel(grid; forcing, params), CUDABackend())
            run!(Simulation(cpu_model; verbose=false))
            run!(Simulation(gpu_model; verbose=false))
            cpu = cpu_model.state
            gpu = cpu_state(gpu_model)

            @test gpu.snowpack_swe ≈ cpu.snowpack_swe rtol=1e-6
            @test gpu.smb_ice ≈ cpu.smb_ice rtol=1e-6
            @test gpu.runoff ≈ cpu.runoff rtol=1e-6
            @test gpu.pdd_sum ≈ cpu.pdd_sum rtol=1e-6
        end
    end

    @testset "ITMModel matches the Fortran bulk budget" begin
        grid = Grid(2)
        forcing = SnowpackForcing(
            dt_days=[2.0],
            ncol=2,
            air_temperature=reshape([274.0, 274.0], 2, 1),
            snowfall_rate=reshape([0.0, 0.0], 2, 1),
            rainfall_rate=reshape([0.0, 0.0], 2, 1),
            shortwave_down=reshape([0.0, 0.0], 2, 1),
            q_sw_net=reshape([400.0, 400.0], 2, 1),
            latitude_deg=reshape([72.0, 65.0], 2, 1),
            surface_height=reshape([1500.0, 500.0], 2, 1),
            ice_thickness=reshape([1000.0, 0.0], 2, 1),
            annual_pdd=reshape([200.0, 800.0], 2, 1),
        )
        model = ITMModel(grid; forcing)
        params = model.params
        sim = Simulation(model; verbose=false)
        @test sim.model isa ITMModel
        @test model.state.H_snow == fill(params.H_snow_max, 2)

        @test finished(run!(sim))
        state = model.state
        @test all(state.melt .> 0.0)
        @test all(state.H_snow .<= params.H_snow_max)
        @test state.alb_s[1] == params.alb_snow_wet
        @test state.smb_ice[1] ≈ state.smbi[1] * 2.0 atol=1e-12

        no_qsw = SnowpackForcing(
            dt_days=[1.0], air_temperature=[274.0], snowfall_rate=[0.0], rainfall_rate=[0.0],
            shortwave_down=[0.0], latitude_deg=[72.0], surface_height=[1500.0],
            ice_thickness=[1000.0], annual_pdd=[200.0],
        )
        no_qsw_model = ITMModel(Grid(1); forcing=no_qsw)
        @test finished(run!(Simulation(no_qsw_model; verbose=false)))
        @test no_qsw_model.state.melt[1] == 0.0

        incomplete = SnowpackForcing(
            dt_days=[1.0], air_temperature=[274.0], snowfall_rate=[0.0], rainfall_rate=[0.0],
            shortwave_down=[0.0], latitude_deg=[72.0], surface_height=[1500.0],
        )
        @test _captured_exception(() -> ITMModel(Grid(1); forcing=incomplete)) isa Exception

        albedo = Chion._itm_surface_albedo(params, 1500.0, 1000.0, 0.0, 200.0)
        @test albedo == params.alb_ice
    end

    @testset "PDD and ITM run in Float32" begin
        forcing() = SnowpackForcing(
            dt_days=fill(1.0, 4),
            ncol=2,
            air_temperature_c=[-3.0 1.0 4.0 2.0; 2.0 5.0 6.0 3.0],
            snowfall_mm_day=[4.0 1.0 0.0 2.0; 1.0 0.0 0.0 0.5],
            rainfall_mm_day=[0.0 0.5 1.0 0.0; 0.0 1.0 2.0 0.0],
            shortwave_down=fill(250.0, 2, 4),
            latitude_deg=[72.0, 65.0],
            surface_height=[1500.0, 500.0],
            ice_thickness=[1000.0, 0.0],
            annual_pdd=[200.0, 800.0],
        )
        for make in (
                (NF) -> PDDModel(Grid(2; NF); forcing=forcing(), params=PDDParams(; NF, method=PISMPDD())),
                (NF) -> ITMModel(Grid(2; NF); forcing=forcing(), params=ITMParams(; NF)),
            )
            model64 = make(Float64)
            model32 = make(Float32)
            @test Chion.number_type(model32.params) == Float32
            run!(Simulation(model64; verbose=false))
            @test finished(run!(Simulation(model32; verbose=false)))
            for name in fieldnames(typeof(model32.state))
                values32 = getfield(model32.state, name)
                @test eltype(values32) == Float32
                @test values32 ≈ getfield(model64.state, name) rtol=1e-4 atol=1e-3
            end
        end
    end
end
