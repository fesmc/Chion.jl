using Test
using Chion

if Chion.cuda_available()
    @testset "GPU execution" begin
        @testset "BESSI transposed storage matches CPU" begin
            ncol = 4
            grid = SnowpackGrid(ncol)
            forcing = SnowpackForcing(
                dt_days=[1.0, 1.0],
                ncol=ncol,
                air_temperature_c=[-8.0, -6.0],
                snowfall_mm_day=[2.0, 1.0],
                rainfall_mm_day=0.0,
                shortwave_down=[100.0, 120.0],
                latitude_deg=fill(70.0, ncol),
            )
            model = BESSIModel(grid; Ntot=4, albedo=:aging)
            cpu_sim = Simulation(model; forcing=forcing, backend=:threads, years=1)
            gpu_sim = Simulation(model; forcing=forcing, backend=:gpu, years=1)

            integrator = init_integrator(gpu_sim; io=devnull)
            gpu_state = integrator.model_runtime.data.state
            @test gpu_state.mass isa Chion.TransposedLayerMatrix
            @test size(gpu_state.mass.parent) == (ncol, model.Ntot)
            @test integrator.model_runtime.data.step_fields.q_lw_down isa Chion.ConstantForcingMatrix
            @test integrator.model_runtime.data.step_fields.latitude_deg isa Chion.ColumnForcingMatrix
            forcing.air_temperature[:, 1] .= 264.0
            forcing.air_pressure .*= 0.99
            sync_forcing!(integrator)
            @test Array(integrator.model_runtime.data.step_fields.air_temperature.values) ==
                forcing.air_temperature.values
            @test Array(integrator.model_runtime.data.step_fields.air_pressure) ==
                forcing.air_pressure

            cpu_result = run!(cpu_sim; io=devnull)
            run!(integrator)
            gpu_result = finalize!(integrator)

            @test cpu_result.status == gpu_result.status == :complete
            @test gpu_sim.now.mass ≈ cpu_sim.now.mass
            @test gpu_sim.now.mass_w ≈ cpu_sim.now.mass_w
            @test gpu_sim.now.density ≈ cpu_sim.now.density
            @test gpu_sim.now.temperature ≈ cpu_sim.now.temperature
            @test gpu_sim.now.smb_ice ≈ cpu_sim.now.smb_ice
            @test gpu_sim.now.albedo ≈ cpu_sim.now.albedo
            @test gpu_sim.now.snow_age_days ≈ cpu_sim.now.snow_age_days
        end

        @testset "BESSI GPU handles a partial final workgroup" begin
            # The production GrIS mask does not generally contain an exact
            # multiple of the CUDA workgroup size (256) active columns.
            ncol = 257
            forcing = SnowpackForcing(
                dt_days=[1.0],
                ncol=ncol,
                air_temperature_c=-5.0,
                snowfall_mm_day=1.0,
                rainfall_mm_day=0.0,
                shortwave_down=100.0,
                latitude_deg=70.0,
            )
            sim = Simulation(BESSIModel(SnowpackGrid(ncol); Ntot=4, albedo=:aging);
                forcing=forcing, backend=:gpu, years=1)
            @test run!(sim; io=devnull).status == :complete
        end

        @testset "BESSI GPU handles bare-ice columns" begin
            ncol = 257
            forcing = SnowpackForcing(
                dt_days=[1.0],
                ncol=ncol,
                air_temperature_c=2.0,
                snowfall_mm_day=0.0,
                rainfall_mm_day=0.0,
                shortwave_down=250.0,
                latitude_deg=70.0,
            )
            model = BESSIModel(SnowpackGrid(ncol); Ntot=4, albedo=:aging)
            sim = Simulation(model;
                forcing=forcing, backend=:gpu, years=1)
            @test run!(sim; io=devnull).status == :complete
            @test all(isfinite, sim.now.albedo)
        end

        @testset "PDD active indices" begin
            ncol = 3
            forcing = SnowpackForcing(
                dt_days=[1.0, 30.0],
                ncol=ncol,
                air_temperature_c=[-2.0, 2.0],
                snowfall_mm_day=[1.0, 1.0],
                rainfall_mm_day=[0.0, 1.0],
                shortwave_down=100.0,
            )
            sim = Simulation(PDDModel(SnowpackGrid(ncol));
                forcing=forcing,
                backend=:gpu,
                years=1,
            )
            integrator = init_integrator(sim; io=devnull)
            runtime_forcing = integrator.model_runtime.data.step_fields
            @test runtime_forcing.air_temperature !== forcing.air_temperature
            @test runtime_forcing.snowfall_rate !== forcing.snowfall_rate
            @test runtime_forcing.rainfall_rate !== forcing.rainfall_rate
            @test !hasproperty(runtime_forcing, :shortwave_down)
            @test !hasproperty(runtime_forcing, :air_pressure)
            set_active_mask!(integrator, [true, false, true])
            run!(integrator)
            result = finalize!(integrator)

            @test result.status == :complete
            @test sim.now.pdd_sum[1] > 0.0
            @test sim.now.pdd_sum[2] == 0.0
            @test sim.now.pdd_sum[3] > 0.0
        end

        @testset "ITM GPU matches CPU" begin
            ncol = 3
            forcing = SnowpackForcing(
                dt_days=[1.0, 1.0, 1.0],
                ncol=ncol,
                air_temperature=[270.0 272.0 275.0; 268.0 271.0 274.0; 265.0 270.0 276.0],
                snowfall_rate=fill(2.0 / 86_400, ncol, 3),
                rainfall_rate=zeros(ncol, 3),
                shortwave_down=fill(450.0, ncol, 3),
                latitude_deg=[65.0, 70.0, 75.0],
                surface_height=[300.0, 1200.0, 2500.0],
                ice_thickness=[1000.0, 1500.0, 2500.0],
                annual_pdd=[400.0, 150.0, 20.0],
            )
            runs = Dict(backend => Simulation(ITMModel(SnowpackGrid(ncol)); forcing, backend, years=2)
                        for backend in (:threads, :gpu))
            for sim in values(runs)
                @test run!(sim; io=devnull).status == :complete
            end
            cpu, gpu = runs[:threads].now, runs[:gpu].now
            @test any(cpu.melt_cum .> 0.0)
            for name in (:H_snow, :smb_cum, :smb_ice, :melt_cum, :runoff_cum, :refreezing_cum, :Tsrf)
                @test Array(getproperty(gpu, name)) ≈ Array(getproperty(cpu, name)) rtol=1e-10
            end
        end
    end
else
    @info "Skipping GPU tests because CUDA is not functional on this node."
end
