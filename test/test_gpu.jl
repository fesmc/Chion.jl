using Test
using CUDA
using Chion

if Chion.cuda_available()
    @testset "GPU execution" begin
        gpu = CUDABackend()

        @testset "BESSI transposed storage matches CPU" begin
            ncol = 4
            forcing = SnowpackForcing(dt_days = [1.0, 1.0], ncol = ncol, air_temperature_c = [-8.0, -6.0],
                snowfall_mm_day = [2.0, 1.0], rainfall_mm_day = 0.0, shortwave_down = [100.0, 120.0],
                latitude_deg = fill(70.0, ncol))
            grid = Grid(ncol; layering = Layering(Ntot = 4))
            params = Params(albedo_scheme = AgingAlbedo())
            cpu = BESSIModel(grid; forcing, params)
            device = BESSIModel(grid; forcing, params, backend = gpu)
            @test device.state.mass isa Chion.TransposedLayerMatrix
            @test size(device.state.mass.parent) == (ncol, 4)
            @test device.forcing.fields.q_lw_down isa Chion.ConstantForcingMatrix
            @test device.forcing.fields.latitude_deg isa Chion.ColumnForcingMatrix
            @test device.forcing.fields.air_temperature.values isa CuArray   # compact: varies in time only
            # A changed input reaches the device.
            set_input!(device, :latitude_deg, fill(72.0, ncol))
            set_input!(cpu, :latitude_deg, fill(72.0, ncol))
            @test Array(device.forcing.fields.latitude_deg.values) == fill(72.0, ncol)

            run!(Simulation(cpu; verbose = false))
            run!(Simulation(device; verbose = false))
            host = cpu_state(device)
            for name in (:mass, :mass_w, :density, :temperature, :smb_ice, :albedo, :snow_age_days)
                @test getfield(host, name) ≈ getfield(cpu.state, name)
            end
            @test Array(device.thickness) ≈ cpu.thickness
            # Back to the CPU.
            back = to_backend(device, CPU())
            @test back.state.mass == host.mass && back.state.mass isa Matrix
        end

        @testset "BESSI GPU handles a partial final workgroup and bare ice" begin
            # The production GrIS mask does not generally contain an exact
            # multiple of the CUDA workgroup size (256) active columns.
            ncol = 257
            for (temperature, snowfall, shortwave) in ((-5.0, 1.0, 100.0), (2.0, 0.0, 250.0))
                forcing = SnowpackForcing(dt_days = [1.0], ncol = ncol, air_temperature_c = temperature,
                    snowfall_mm_day = snowfall, rainfall_mm_day = 0.0, shortwave_down = shortwave, latitude_deg = 70.0)
                model = BESSIModel(Grid(ncol; layering = Layering(Ntot = 4)); forcing,
                    params = Params(albedo_scheme = AgingAlbedo()), backend = gpu)
                sim = run!(Simulation(model; verbose = false))
                @test finished(sim)
                @test all(isfinite, Array(model.albedo))
            end
        end

        @testset "Per-column boundary conditions on the GPU" begin
            forcing = SnowpackForcing(dt_days = ones(20), ncol = 3, air_temperature_c = -20.0, snowfall_mm_day = 1.0,
                rainfall_mm_day = 0.0, shortwave_down = 50.0, latitude_deg = 70.0)
            params = Params(temperature_init = 250.0)
            boundary = BoundaryConditions(bottom = GeothermalFlux([0.2, 0.5, 1.0]))
            cpu = BESSIModel(Grid(3); forcing, params, boundary)
            device = BESSIModel(Grid(3); forcing, params, boundary, backend = gpu)
            @test device.boundary.bottom.flux isa CuArray
            integrate!(cpu, 1:20); integrate!(device, 1:20)
            @test cpu_state(device).ice_temperature ≈ cpu.state.ice_temperature rtol = 1e-10
            set_boundary!(device, :bottom, [1.0, 1.0, 1.0])
            @test Array(device.boundary.bottom.flux) == [1.0, 1.0, 1.0]
            sim = Simulation(device; verbose = false, output = OutputConfig(dir = mktempdir()))
            @test haskey(Chion._run_metadata(sim)["boundary"]["bottom"], "flux")
        end

        @testset "PDD active indices" begin
            ncol = 3
            forcing = SnowpackForcing(dt_days = [1.0, 30.0], ncol = ncol, air_temperature_c = [-2.0, 2.0],
                snowfall_mm_day = [1.0, 1.0], rainfall_mm_day = [0.0, 1.0], shortwave_down = 100.0)
            model = PDDModel(Grid(ncol); forcing, backend = gpu)
            fields = model.forcing.fields
            @test fields.air_temperature.values isa CuArray
            @test !hasproperty(fields, :shortwave_down) && !hasproperty(fields, :air_pressure)
            sim = Simulation(model; verbose = false)
            set_active_mask!(sim, [true, false, true])
            run!(sim)
            pdd = Array(model.pdd_sum)
            @test pdd[1] > 0.0 && pdd[2] == 0.0 && pdd[3] > 0.0
        end

        @testset "ITM GPU matches CPU" begin
            ncol = 3
            forcing = SnowpackForcing(dt_days = [1.0, 1.0, 1.0], ncol = ncol,
                air_temperature = [270.0 272.0 275.0; 268.0 271.0 274.0; 265.0 270.0 276.0],
                snowfall_rate = fill(2.0 / 86_400, ncol, 3), rainfall_rate = zeros(ncol, 3),
                shortwave_down = fill(450.0, ncol, 3), latitude_deg = [65.0, 70.0, 75.0],
                surface_height = [300.0, 1200.0, 2500.0], ice_thickness = [1000.0, 1500.0, 2500.0],
                annual_pdd = [400.0, 150.0, 20.0])
            models = Dict(backend => ITMModel(Grid(ncol); forcing, backend) for backend in (CPU(), gpu))
            for model in values(models)
                sim = Simulation(model; years = [2001, 2001], verbose = false)
                run!(sim)
                @test finished(sim)
            end
            cpu, device = models[CPU()].state, cpu_state(models[gpu])
            @test any(cpu.melt_cum .> 0.0)
            for name in (:H_snow, :smb_cum, :smb_ice, :melt_cum, :runoff_cum, :refreezing_cum, :Tsrf)
                @test getfield(device, name) ≈ getfield(cpu, name) rtol = 1e-10
            end
        end
    end
else
    @info "Skipping GPU tests because CUDA is not functional on this node."
end
