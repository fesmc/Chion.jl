using Test
using Chion

@testset "Dynamic albedo wetness elapsed time" begin
    for NF in (Float32, Float64)
        c = Chion.SnowpackPhysicalConstants(NF; albedo_scheme=:dynamic)
        mass = fill(NF(100), 1, 1)
        density = fill(NF(300), 1, 1)
        # Disable temperature aging to isolate wetness subdivision invariance.
        temperature = fill(c.T0 - NF(30), 1, 1)
        pore_volume = mass[1] / density[1] - mass[1] / c.rho_i
        function evolve(r, steps)
            water = fill(NF(r) * c.max_lwc_albedo * pore_volume * c.rho_w, 1, 1)
            albedo = [c.alpha_dry]
            for dt in steps
                Chion._update_surface_albedo_arrays!(
                    [1], mass, water, density, temperature, albedo, 1, c, NF(dt))
            end
            only(albedo)
        end
        for r in (0, 0.1, 0.5, 0.99, 1, 2)
            daily = evolve(r, [1])
            expected = c.alpha_dry - (c.alpha_dry - c.alpha_wet) * NF(clamp(r, 0, 1))
            @test daily ≈ expected
            @test evolve(r, fill(1/24, 24)) ≈ daily
            @test evolve(r, [0.1, 0.2, 0.7]) ≈ daily
            @test evolve(r, [2]) ≈ evolve(r, [1, 1])
            @test evolve(r, [0]) == c.alpha_dry
            @test c.alpha_wet <= daily <= c.alpha_dry
        end
    end
end
