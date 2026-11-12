using Test
using Dates
using Chion

# A warm, snow-free column ablating at the melting point (the Fortran chion's
# test_seb_fluxes, "ablate"): 40 days at 280 K, light rain, 320 W m⁻².
function _ablating_column(; nday=40, kwargs...)
    model = BESSIModel(SnowpackGrid(1); kwargs...)
    shortwave = 320.2
    forcing = SnowpackForcing(
        dt_days=fill(1.0, nday),
        air_temperature=fill(280.02, 1, nday),
        snowfall_rate=fill(0.0, 1, nday),
        rainfall_rate=fill(1.0e-6, 1, nday),
        shortwave_down=fill(shortwave, 1, nday),
        wind_speed=5.0,
        relative_humidity=0.8,
        latitude_deg=69.95,
        surface_height=1005.0,
        time_values=[DateTime(2001, 6, 1, 12) + Day(k - 1) for k in 1:nday],
    )
    return model, forcing, shortwave
end

const _SEB_SUMS = (:shortwave_down_sum, :shortwave_net_sum, :longwave_down_sum, :longwave_net_sum,
    :sensible_heat_flux_sum, :latent_heat_flux_sum, :rain_heat_flux_sum)

_seb_sums(state) = NamedTuple{_SEB_SUMS}(map(name -> getfield(state, name)[1], _SEB_SUMS))

@testset "Surface energy balance diagnostics" begin
    @testset "melt energy closes on the sum of the terms (ice substrate=$(nice))" for nice in (5, 0)
        model, forcing, shortwave = _ablating_column(; ice_substrate_layers=nice)
        sim = Simulation(model; forcing=forcing, years=1)
        integrator = init_integrator(sim; io=devnull)
        state = integrator.model_runtime.data.state
        nday = length(forcing.dt_days)
        ok_shortwave = true
        sums_half = nothing
        melt_half = 0.0
        for k in 1:nday
            before = _seb_sums(state)
            step!(integrator)
            after = _seb_sums(state)
            dt = forcing.dt_days[k]
            sw_down = (after.shortwave_down_sum - before.shortwave_down_sum) / dt
            sw_net = (after.shortwave_net_sum - before.shortwave_net_sum) / dt
            ok_shortwave &= isapprox(sw_down, shortwave; rtol=1e-3) && 0 < sw_net <= sw_down
            if k == nday ÷ 2
                sums_half = after
                melt_half = state.melt[1]
            end
        end
        finalize!(integrator)
        @test all(iszero, state.mass)
        @test ok_shortwave

        # Over the second half, once the ice substrate has warmed (1e-4 with
        # the substrate, exact without).
        sums = _seb_sums(state)
        energy_in = sum(sums[name] - sums_half[name] for name in
            (:shortwave_net_sum, :longwave_net_sum, :sensible_heat_flux_sum, :latent_heat_flux_sum, :rain_heat_flux_sum))
        melt_energy = (state.melt[1] - melt_half) * model.c.Lm
        @test melt_energy > 0
        @test melt_energy ≈ energy_in * model.c.seconds_per_day rtol=1e-3
        @test sums.longwave_down_sum > 0
    end

    @testset "sums start at zero" begin
        model, forcing, _ = _ablating_column(; nday=2)
        state = BESSIState(model)
        @test all(name -> iszero(getfield(state, name)[1]), _SEB_SUMS)
        sim = Simulation(model; forcing=forcing, years=1)
        run!(sim; io=devnull)
        @test sim.now.shortwave_down_sum[1] > 0
        @test sim.ref.shortwave_down_sum[1] == 0
    end
end
