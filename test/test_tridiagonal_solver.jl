import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test
using Chion

const SM = Chion.SnowpackModel

function random_diagonally_dominant_tridiagonal(n::Int)
    lower_diagonal = rand(n - 1)
    upper_diagonal = rand(n - 1)
    main_diagonal = rand(n) .+ 2.0
    right_hand_side = rand(n)
    return lower_diagonal, main_diagonal, upper_diagonal, right_hand_side
end

@testset "Custom Tridiagonal Solver" begin
    @testset "Thomas Solver Matches LinearAlgebra" begin
        for n in (3, 5, 10, 15, 30)
            lower_diagonal, main_diagonal, upper_diagonal, right_hand_side =
                random_diagonally_dominant_tridiagonal(n)

            linear_solution = SM._solve_tridiagonal_system(
                lower_diagonal,
                main_diagonal,
                upper_diagonal,
                right_hand_side;
                solver=:linear_algebra,
            )
            thomas_solution = SM._solve_tridiagonal_system(
                lower_diagonal,
                main_diagonal,
                upper_diagonal,
                right_hand_side;
                solver=:thomas,
            )

            @test thomas_solution ≈ linear_solution atol = 1e-12 rtol = 1e-12
        end
    end

    @testset "Energy Flux Supports Both Solver Backends" begin
        column_linear = SM.SnowpackColumn(Ntot = 10, N = 5, temperature_init = 260.0, density_init = 320.0)
        column_thomas = SM.SnowpackColumn(Ntot = 10, N = 5, temperature_init = 260.0, density_init = 320.0)

        column_linear.mass[1:5] .= [200.0, 220.0, 240.0, 260.0, 280.0]
        column_linear.density[1:5] .= [300.0, 330.0, 360.0, 390.0, 420.0]
        column_linear.temperature[1:5] .= [270.0, 265.0, 260.0, 255.0, 250.0]

        column_thomas.mass .= column_linear.mass
        column_thomas.density .= column_linear.density
        column_thomas.temperature .= column_linear.temperature

        out_linear = SM.go_energy_flux!(
            column_linear,
            260.0,
            0.0,
            0.0,
            0.0,
            3600.0;
            q_sw_net=0.0,
            q_lw_down=0.0,
            q_sh=0.0,
            q_lh=0.0,
            tridiagonal_solver=:linear_algebra,
        )

        out_thomas = SM.go_energy_flux!(
            column_thomas,
            260.0,
            0.0,
            0.0,
            0.0,
            3600.0;
            q_sw_net=0.0,
            q_lw_down=0.0,
            q_sh=0.0,
            q_lh=0.0,
            tridiagonal_solver=:thomas,
        )

        @test column_thomas.temperature[1:5] ≈ column_linear.temperature[1:5] atol = 1e-12 rtol = 1e-12
        @test out_thomas.heating ≈ out_linear.heating atol = 1e-12 rtol = 1e-12
        @test out_thomas.needs_melt == out_linear.needs_melt
    end
end
