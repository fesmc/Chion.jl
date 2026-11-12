using Test
using Chion

# Mirrors the Fortran chion's tests/test_forcing_host.f90.
@testset "Host forcing helpers" begin
    @testset "shortwave from a direct transmissivity" begin
        greenland = (trans0=0.587, trans_per_km=0.052)
        @test shortwave_down_from_transmissivity(500.0, 1500.0; greenland...) ≈ 500.0 * (0.587 + 0.052 * 1.5) rtol=1e-12
        @test shortwave_down_from_transmissivity(500.0, -50.0; greenland...) ≈ 500.0 * 0.587 rtol=1e-12
        @test shortwave_down_from_transmissivity(500.0, NaN; greenland...) ≈ 500.0 * 0.587 rtol=1e-12
        @test shortwave_down_from_transmissivity(500.0, 1.0e4; trans0=0.68, trans_per_km=0.041) == 500.0
        @test shortwave_down_from_transmissivity(500.0, 1000.0; trans0=-0.2, trans_per_km=0.1) == 0.0
    end

    @testset "wind from elevation and slope" begin
        u = wind_speed_from_geometry(fill(2000.0, 40, 30), 16.0e3; u0=2.5, u_z=1.27)
        @test u[10, 10] ≈ 2.5 + 1.27 * 2.0 rtol=1e-12
        # 20 m/km along x, on a 16 km grid and on one twice as coarse.
        z = [20.0 * 16.0 * (i - 1) for i in 1:40, j in 1:30]
        z2 = [20.0 * 32.0 * (i - 1) for i in 1:20, j in 1:15]
        u = wind_speed_from_geometry(z, 16.0e3; u0=3.0, u_s=0.2)
        u2 = wind_speed_from_geometry(z2, 32.0e3; u0=3.0, u_s=0.2)
        @test u[20, 15] ≈ 3.0 + 0.2 * 20.0 rtol=1e-12
        @test u2[10, 8] ≈ u[20, 15] rtol=1e-12
        @test u[1, 1] ≈ 3.0 + 0.2 * 20.0 rtol=1e-12   # one-sided at the edge
        @test all(iszero, wind_speed_from_geometry(zero(z), 16.0e3; u0=-1.0))
    end

    @testset "sigma of the daily temperature" begin
        era5 = (sigma0=2.0, slope=-0.14, smin=2.0, smax=6.0)
        @test temperature_variability_sigma(263.15, 273.15; era5...) ≈ 3.4 rtol=1e-12
        @test temperature_variability_sigma(278.15, 273.15; era5...) == 2.0
        @test temperature_variability_sigma(223.15, 273.15; era5...) == 6.0
    end

    @testset "red noise, zero monthly mean" begin
        xi = red_noise_monthly(360, 1, 30, 4.0, 1)
        @test length(xi) == 360
        @test all(abs(sum(xi[30m+1:30m+30])) / 30 < 1e-12 for m in 0:11)
        @test sum(abs2, xi) / 360 ≈ 1.0 rtol=1e-12
        r1 = sum(xi[1:359] .* xi[2:360]) / sum(abs2, xi)
        @test abs(r1 - exp(-0.25)) < 0.15
        @test red_noise_monthly(360, 1, 30, 4.0, 1) == xi
        @test red_noise_monthly(360, 1, 30, 4.0, 2) != xi
        xb = red_noise_monthly(360, 5, 30, 4.0, 1)
        @test length(xb) == 72
        @test all(abs(sum(xb[6m+1:6m+6])) / 6 < 1e-12 for m in 0:11)
        # The Fortran chion's sequence (single precision) for the same seed.
        @test xi[1:3] ≈ [-0.7955155, 0.0267786, -1.6122738] atol=1e-6
    end
end
