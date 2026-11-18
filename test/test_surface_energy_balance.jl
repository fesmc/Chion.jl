using Test
using Dates
using NCDatasets
using Chion

# A warm, snow-free column ablating at the melting point (the Fortran chion's
# test_seb_fluxes, "ablate"): 40 days at 280 K, light rain, 320 W m⁻².
function _ablating_column(; nday = 40, ice_substrate_layers = 5)
    grid = Grid(1; layering = Layering(; ice_substrate_layers), x = [0.0], y = [0.0], js = [1], is = [1], mask = ones(1, 1))
    shortwave = 320.2
    forcing = SnowpackForcing(
        dt_days = fill(1.0, nday),
        air_temperature = fill(280.02, 1, nday),
        snowfall_rate = fill(0.0, 1, nday),
        rainfall_rate = fill(1.0e-6, 1, nday),
        shortwave_down = fill(shortwave, 1, nday),
        wind_speed = 5.0,
        relative_humidity = 0.8,
        latitude_deg = 69.95,
        surface_height = 1005.0,
        time_values = [DateTime(2001, 6, 1, 12) + Day(k - 1) for k in 1:nday],
    )
    return BESSIModel(grid; forcing), forcing, shortwave
end

const _SEB_SUMS = (:shortwave_down_sum, :shortwave_net_sum, :longwave_down_sum, :longwave_net_sum,
    :sensible_heat_flux_sum, :latent_heat_flux_sum, :rain_heat_flux_sum)

_seb_sums(state) = NamedTuple{_SEB_SUMS}(map(name -> getfield(state, name)[1], _SEB_SUMS))

@testset "Surface energy balance diagnostics" begin
    @testset "melt energy closes on the sum of the terms (ice substrate=$(nice))" for nice in (5, 0)
        model, forcing, shortwave = _ablating_column(; ice_substrate_layers = nice)
        sim = Simulation(model; verbose = false)
        state = model.state
        nday = length(forcing.dt_days)
        ok_shortwave = true
        sums_half = nothing
        melt_half = 0.0
        for k in 1:nday
            before = _seb_sums(state)
            time_step!(sim)
            after = _seb_sums(state)
            dt = forcing.dt_days[k]
            sw_down = (after.shortwave_down_sum - before.shortwave_down_sum) / dt
            sw_net = (after.shortwave_net_sum - before.shortwave_net_sum) / dt
            ok_shortwave &= isapprox(sw_down, shortwave; rtol = 1e-3) && 0 < sw_net <= sw_down
            if k == nday ÷ 2
                sums_half = after
                melt_half = state.melt[1]
            end
        end
        @test finished(sim)
        @test all(iszero, state.mass)
        @test ok_shortwave

        # Over the second half, once the ice substrate has warmed (1e-4 with
        # the substrate, exact without).
        sums = _seb_sums(state)
        energy_in = sum(sums[name] - sums_half[name] for name in
            (:shortwave_net_sum, :longwave_net_sum, :sensible_heat_flux_sum, :latent_heat_flux_sum, :rain_heat_flux_sum))
        melt_energy = (state.melt[1] - melt_half) * model.Lm
        @test melt_energy > 0
        @test melt_energy ≈ energy_in * model.seconds_per_day rtol = 1e-3
        @test sums.longwave_down_sum > 0
    end

    @testset "monthly means in the monthly output" begin
        model, forcing, shortwave = _ablating_column()   # 30 days of June, 10 of July
        mktempdir() do dir
            path = joinpath(dir, "monthly_seb.nc")
            sim = Simulation(model; outputs = Output(path; variables = :monthly, frequency = :monthly), verbose = false)
            run!(sim)
            @test finished(sim)
            NCDataset(path) do ds
                days = [30.0, 10.0]
                for (flux, flux_sum) in Chion.MONTHLY_FLUX_SUMS
                    @test ds[String(flux)].attrib["units"] == "W m-2"
                    means = Float64.(ds[String(flux)][1, 1, 1:2])
                    @test sum(means .* days) ≈ getfield(model.state, flux_sum)[1] rtol = 1e-6
                end
                @test Float64.(ds["shortwave_down"][1, 1, 1:2]) ≈ [shortwave, shortwave] rtol = 1e-6
            end
        end
    end

    @testset "sums start at zero" begin
        model, _, _ = _ablating_column(; nday = 2)
        @test all(name -> iszero(getfield(model.state, name)[1]), _SEB_SUMS)
        run!(Simulation(model; verbose = false))
        @test model.state.shortwave_down_sum[1] > 0
        @test initial_state(model.grid, model.params).shortwave_down_sum[1] == 0
    end
end
