using Test
using Chion

# One cold, snow-free column on an ice substrate (or one snow layer), without diurnal cycle,
# shortwave, rain or humidity: only the boundary conditions move heat.
function _bc_model(boundary; ndays = 30, snowfall_mm_day = 0.0, layering = Layering(), temperature_init = 250.0)
    forcing = SnowpackForcing(; dt_days = ones(ndays), ncol = 1, air_temperature_c = -20.0,
        snowfall_mm_day, rainfall_mm_day = 0.0, shortwave_down = 0.0, latitude_deg = 70.0)
    params = Params(; temperature_init, diurnal_cycle = NoDiurnalCycle())
    return BESSIModel(Grid(1; layering); forcing, params, boundary)
end

_enthalpy(model) = column_enthalpy(model, 1)
_boundary_enthalpy(model) = model.cache.workspace.boundary_enthalpy[1]
_temperatures(model) = vcat(model.temperature[1:model.N[1], 1], model.ice_temperature[:, 1])

@testset "Boundary conditions" begin
    seconds = 30 * 86_400.0

    @testset "prescribed fluxes add exactly their energy" begin
        # Neumann surface: ΔH = F t.
        model = _bc_model(BoundaryConditions(surface = PrescribedSurfaceFlux(5.0)))
        H0 = _enthalpy(model)
        integrate!(model, 1:30)
        @test _enthalpy(model) - H0 ≈ 5.0 * seconds rtol = 1e-9
        @test _boundary_enthalpy(model) ≈ 5.0 * seconds rtol = 1e-9
        @test all(<(model.T0), _temperatures(model))
        # Neumann bottom under an insulated surface: ΔH = q t, warming from below.
        model = _bc_model(BoundaryConditions(surface = PrescribedSurfaceFlux(0.0), bottom = GeothermalFlux(0.5)))
        H0 = _enthalpy(model)
        integrate!(model, 1:30)
        @test _enthalpy(model) - H0 ≈ 0.5 * seconds rtol = 1e-9
        T = _temperatures(model)
        @test issorted(T) && T[end] > 250.0
    end

    @testset "budget closes with the energy balance and snowfall" begin
        for bottom in (GeothermalFlux(0.05), PrescribedBottomTemperature(265.0))
            model = _bc_model(BoundaryConditions(; bottom); snowfall_mm_day = 4.0)
            H0, B0 = _enthalpy(model), _boundary_enthalpy(model)
            for k in 1:30
                integrate!(model, k:k)
                @test _enthalpy(model) - H0 ≈ _boundary_enthalpy(model) - B0 rtol = 1e-9
            end
            @test model.N[1] > 0
        end
    end

    @testset "Dirichlet at both ends relaxes to their temperature" begin
        # Ice substrate (5 rows) and a single snow row (the surface row is also the bottom row).
        for layering in (Layering(), Layering(Ntot = 1, ice_substrate_layers = 0))
            boundary = BoundaryConditions(surface = PrescribedSurfaceTemperature(260.0),
                bottom = PrescribedBottomTemperature(260.0))
            model = _bc_model(boundary; ndays = 365, layering)
            if layering.ice_substrate_layers == 0
                model.N[1] = 1; model.mass[1, 1] = 200.0; model.density[1, 1] = 400.0; model.temperature[1, 1] = 250.0
            end
            H0, B0 = _enthalpy(model), _boundary_enthalpy(model)
            integrate!(model, 1:365)
            @test all(T -> isapprox(T, 260.0; atol = 1e-6), _temperatures(model))
            @test model.Tsrf[1] == 260.0
            @test _enthalpy(model) - H0 ≈ _boundary_enthalpy(model) - B0 rtol = 1e-9
        end
        # Different ends: a steady profile between them, warming with depth.
        model = _bc_model(BoundaryConditions(surface = PrescribedSurfaceTemperature(250.0),
            bottom = PrescribedBottomTemperature(270.0)); ndays = 365)
        integrate!(model, 1:365)
        T = _temperatures(model)
        @test issorted(T) && 250.0 < T[1] && T[end] < 270.0
    end

    @testset "one value per column" begin
        # Three columns with their own geothermal flux: each column gains its own q t.
        q = [0.2, 0.5, 1.0]
        forcing = SnowpackForcing(; dt_days = ones(30), ncol = 3, air_temperature_c = -20.0, snowfall_mm_day = 0.0,
            rainfall_mm_day = 0.0, shortwave_down = 0.0, latitude_deg = 70.0)
        params = Params(; temperature_init = 250.0, diurnal_cycle = NoDiurnalCycle())
        model = BESSIModel(Grid(3); forcing, params,
            boundary = BoundaryConditions(surface = PrescribedSurfaceFlux(0.0), bottom = GeothermalFlux(q)))
        H0 = [column_enthalpy(model, c) for c in 1:3]
        integrate!(model, 1:30)
        @test [column_enthalpy(model, c) for c in 1:3] .- H0 ≈ q .* seconds rtol = 1e-9
        # Per-column Dirichlet ends, updated with set_boundary! (e.g. from an ice-sheet model).
        model = BESSIModel(Grid(3); forcing = SnowpackForcing(; dt_days = ones(365), ncol = 3, air_temperature_c = -20.0,
                snowfall_mm_day = 0.0, rainfall_mm_day = 0.0, shortwave_down = 0.0, latitude_deg = 70.0), params,
            boundary = BoundaryConditions(surface = PrescribedSurfaceTemperature([255.0, 260.0, 265.0]),
                bottom = PrescribedBottomTemperature(fill(255.0, 3))))
        set_boundary!(model, :bottom, [255.0, 260.0, 265.0])
        integrate!(model, 1:365)
        @test all(c -> all(T -> isapprox(T, (255.0, 260.0, 265.0)[c]; atol = 1e-6), model.ice_temperature[:, c]), 1:3)
        @test_throws ArgumentError BESSIModel(Grid(3); forcing, params, boundary = BoundaryConditions(bottom = GeothermalFlux([1.0, 2.0])))
        @test_throws ArgumentError set_boundary!(model, :bottom, [1.0, 2.0])
        @test_throws ArgumentError set_boundary!(model, :bottom, [255.0, -1.0, 265.0])
    end

    @testset "construction, TOML and metadata" begin
        @test BoundaryConditions() == BoundaryConditions(SurfaceEnergyBalance(), Insulated())
        model = _bc_model(BoundaryConditions(bottom = GeothermalFlux(0.06)))
        @test model.boundary.bottom === GeothermalFlux(0.06)
        model32 = BESSIModel(Grid(1; NF = Float32); forcing = model.forcing.source, params = Params(NF = Float32),
            boundary = BoundaryConditions(bottom = GeothermalFlux(0.06)))
        @test model32.boundary.bottom === GeothermalFlux(0.06f0)
        @test_throws ArgumentError _bc_model(BoundaryConditions(surface = PrescribedSurfaceTemperature(-1.0)))
        @test_throws ArgumentError _bc_model(BoundaryConditions(bottom = GeothermalFlux(Inf)))
        table = Dict{String, Any}("surface" => "energy_balance",
            "bottom" => Dict{String, Any}("type" => "geothermal_flux", "flux" => 0.07))
        @test Chion._config_boundary(table) == BoundaryConditions(bottom = GeothermalFlux(0.07))
        @test Chion._config_boundary(Dict{String, Any}()) == BoundaryConditions()
        sim = Simulation(model; verbose = false)
        meta = Chion._run_metadata(sim)
        @test Chion._config_boundary(meta["boundary"]) == model.boundary
        @test sprint(show, PrescribedBottomTemperature(265.0)) == "PrescribedBottomTemperature(temperature = 265.0)"
    end
end
