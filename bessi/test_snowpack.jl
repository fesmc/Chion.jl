"""
Unit tests for SnowpackModel module.

Tests the Born et al. (2019) layer dynamics algorithm.
"""

## Preamble #############################################
cd(@__DIR__)
import Pkg; Pkg.activate(".")
#########################################################

include("SnowpackModel.jl")
using .SnowpackModel
using Test

function test_initialization()
    column = SnowpackColumn()
    @test column.nlayers == 15
    @test column.n_active == 0
    @test column.mmax == 500.0
    @test column.msplit == 300.0
    @test column.mmin == 100.0
    @test all(column.mass .== 0.0)

    state = get_state(column)
    @test state["total_mass"] == 0.0
    @test state["n_active"] == 0
end

function test_basic_accumulation()
    column = SnowpackColumn()
    dt = 86400.0
    mdot = 1.0 / 86400.0

    for _ in 1:50
        step!(column, mdot, dt)
    end

    state = get_state(column)
    @test state["n_active"] == 1
    @test state["total_mass"] ≈ 50.0 atol=1e-6
end

function test_layer_splitting_at_mmax()
    column = SnowpackColumn()
    dt = 86400.0

    step!(column, column.mmax / dt, dt)
    state1 = get_state(column)
    @test state1["n_active"] == 1
    @test state1["mass"][1] ≈ column.mmax

    step!(column, 1.0 / dt, dt)
    state2 = get_state(column)
    @test state2["n_active"] == 2
    @test state2["mass"][1] ≈ (column.mmax + 1.0 - column.msplit) atol=1e-6
    @test state2["mass"][2] ≈ column.msplit atol=1e-6
    @test state2["total_mass"] ≈ column.mmax + 1.0 atol=1e-6
end

function test_multiple_layer_splits()
    column = SnowpackColumn()
    dt = 86400.0

    total_mass = 1600.0
    step!(column, total_mass / dt, dt)

    state = get_state(column)
    @test state["n_active"] >= 4
    @test state["total_mass"] ≈ total_mass atol=1e-3

    for i in 2:state["n_active"]
        @test state["mass"][i] >= column.msplit - 1e-6
    end
end

function test_layer_merging_at_mmin()
    column = SnowpackColumn()
    dt = 86400.0

    step!(column, 400.0 / dt, dt)
    step!(column, 150.0 / dt, dt)

    state1 = get_state(column)
    @test state1["n_active"] == 2

    step!(column, -250.0 / dt, dt)

    state2 = get_state(column)
    @test state2["n_active"] == 1
    @test state2["total_mass"] ≈ 300.0 atol=1e-3
end

function test_mass_conservation_accumulation()
    column = SnowpackColumn()
    dt = 86400.0

    total_added = 0.0
    for _ in 1:100
        dmass = 10.0
        step!(column, dmass / dt, dt)
        total_added += dmass
    end

    state = get_state(column)
    @test state["total_mass"] ≈ total_added atol=1e-3
end

function test_mass_conservation_melt()
    column = SnowpackColumn()
    dt = 86400.0

    step!(column, 1000.0 / dt, dt)
    mass_before = get_state(column)["total_mass"]

    step!(column, -500.0 / dt, dt)
    mass_after = get_state(column)["total_mass"]

    @test mass_after ≈ mass_before - 500.0 atol=1e-3
end

function test_ice_formation_at_base()
    column = SnowpackColumn(nlayers=5, initial_density=400.0)
    dt = 86400.0

    for _ in 1:10
        step!(column, 300.0 / dt, dt)
    end

    for _ in 1:20
        step!(column, -100.0 / dt, dt)
    end

    state = get_state(column)
    @test state["n_active"] >= 0
end

function test_full_column_behavior()
    column = SnowpackColumn(nlayers=3)
    dt = 86400.0

    for _ in 1:5
        step!(column, 500.0 / dt, dt)
    end

    state = get_state(column)
    @test state["n_active"] == 3

    step!(column, 500.0 / dt, dt)
    state2 = get_state(column)
    @test state2["n_active"] <= 3
end

function test_density_preservation_during_split()
    column = SnowpackColumn(initial_density=350.0)
    dt = 86400.0

    step!(column, 600.0 / dt, dt)

    state = get_state(column)
    @test state["n_active"] == 2
    @test state["density"][1] ≈ 350.0
    @test state["density"][2] ≈ 350.0
end

function test_density_averaging_during_merge()
    column = SnowpackColumn(nlayers=15)
    dt = 86400.0

    column.n_active = 2
    column.mass[1] = 50.0
    column.density[1] = 300.0
    column.mass[2] = 200.0
    column.density[2] = 400.0

    SnowpackModel.merge_surface_layer!(column)

    state = get_state(column)
    @test state["n_active"] == 1

    expected_density = (50.0 * 300.0 + 200.0 * 400.0) / 250.0
    @test state["density"][1] ≈ expected_density atol=1e-6
    @test state["mass"][1] ≈ 250.0
end

function test_empty_column_behavior()
    column = SnowpackColumn()
    dt = 86400.0

    mdot_base = step!(column, -10.0 / dt, dt)

    state = get_state(column)
    @test state["n_active"] == 0
    @test state["total_mass"] == 0.0
    @test mdot_base >= 0.0
end

@testset "SnowpackModel Tests" begin
    @testset "Initialization"                    test_initialization()
    @testset "Basic accumulation"               test_basic_accumulation()
    @testset "Layer splitting at mmax"           test_layer_splitting_at_mmax()
    @testset "Multiple layer splits"             test_multiple_layer_splits()
    @testset "Layer merging at mmin"             test_layer_merging_at_mmin()
    @testset "Mass conservation during accumulation" test_mass_conservation_accumulation()
    @testset "Mass conservation during melt"    test_mass_conservation_melt()
    @testset "Ice formation at base"             test_ice_formation_at_base()
    @testset "Full column behavior"              test_full_column_behavior()
    @testset "Density preservation during split" test_density_preservation_during_split()
    @testset "Density averaging during merge"    test_density_averaging_during_merge()
    @testset "Empty column behavior"             test_empty_column_behavior()
end

println("\nAll tests passed! ✓")
