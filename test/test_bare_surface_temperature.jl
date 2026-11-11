using Test, Chion

# A cold, dry column (East Antarctic plateau analogue) started empty. Each day
# lays ~0.08 kg m⁻² of snow on the bare surface. If the surface temperature of
# an emptied column is reset to T0, the next step linearises the surface fluxes
# about 273 K, the thin snow comes out too warm and sublimates within the step,
# and the column stays bare for good with all its snowfall sublimated.
function _cold_plateau_forcing(doy; t_mean=228.15, snowfall_year=30.0, seconds_per_day=86_400.0)
    toa = max(230.0 + 300.0 * cos(2π * (doy - 355) / 365), 0.0)
    return Chion.SnowpackStepForcing(
        t_mean + 15.0 * cos(2π * (doy - 15) / 365), 1.0,
        snowfall_year / (365 * seconds_per_day), 0.0, 0.6 * toa, 5.0;
        relative_humidity=0.8, has_relative_humidity=true, air_pressure=70_000.0,
        latitude_deg=-80.0, surface_height=3000.0, day_of_year=Float64(doy),
        solar_longitude_deg=mod(360.0 * (doy - 80) / 365, 360.0),
    )
end

_column_storage(s) = sum(s.mass[1:s.N[1], 1]) + sum(s.mass_w[1:s.N[1], 1])

@testset "Bare-surface temperature" begin
    @testset "Snow sublimating away over the substrate keeps the solved surface temperature" begin
        model = BESSIModel(SnowpackGrid(1))
        state = BESSIState(model)
        state.ice_temperature .= 240.0
        state.Tsrf[1] = 240.0
        forcing = Chion.SnowpackStepForcing(240.0, 1.0, 1e-7, 0.0, 0.0, 8.0;
            relative_humidity=0.2, has_relative_humidity=true)
        config = Chion._step_config_from_keywords(; Chion._bessi_step_kwargs(model)...)
        Chion.column_step!(state, 1, forcing, config, Chion.ColumnarStepWorkspace(state))
        @test state.N[1] == 0
        @test state.sublimation[1] > 0
        @test abs(state.Tsrf[1] - 240.0) < 10.0
    end

    @testset "Cold, dry column builds snow from empty (substrate layers $n_ice)" for n_ice in (5, 0)
        model = BESSIModel(SnowpackGrid(1); ice_substrate_layers=n_ice)
        state = BESSIState(model)
        workspace = Chion.ColumnarStepWorkspace(state)
        config = Chion._step_config_from_keywords(; Chion._bessi_step_kwargs(model)...)
        snowfall = 0.0
        tsrf_mean = albedo_mean = air_mean = 0.0
        for year in 1:100, doy in 1:365
            forcing = _cold_plateau_forcing(doy; seconds_per_day=model.c.seconds_per_day)
            snowfall += forcing.snowfall_rate * model.c.seconds_per_day
            Chion.column_step!(state, 1, forcing, config, workspace)
            if year == 100
                tsrf_mean += state.Tsrf[1] / 365
                albedo_mean += state.albedo[1] / 365
                air_mean += forcing.air_temperature / 365
            end
        end
        storage = _column_storage(state)
        # Snowfall = storage + runoff + ice transfer - net vapour exchange.
        @test storage + state.runoff[1] + state.smb_ice[1] - state.vapor_mass[1] ≈ snowfall rtol=1e-10
        @test storage > 2 / 3 * snowfall
        @test state.sublimation[1] < snowfall / 3
        @test state.runoff[1] == 0.0
        @test abs(tsrf_mean - air_mean) < 5.0
        @test albedo_mean > 0.75
    end
end
