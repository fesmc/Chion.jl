using Test
using Dates
using NCDatasets
using Chion

function _captured_exception(f::Function)
    try
        f()
        return nothing
    catch err
        return err
    end
end

function _sample_model_forcing_grid()
    grid = SnowpackGrid(4;
        x    = [0.0, 10_000.0],
        y    = [0.0, 10_000.0],
        js   = [1, 1, 2, 2],
        is   = [1, 2, 1, 2],
        mask = ones(2, 2),
    )
    model = BESSIModel(grid; Ntot=4)
    forcing = SnowpackForcing(
        dt_days           = [1.0, 1.0, 1.0],
        ncol              = 4,
        air_temperature_c = [-15.0, -14.0, -13.0],
        snowfall_mm_day   = 1.0,
        rainfall_mm_day   = 0.0,
        shortwave_down    = [120.0, 140.0, 160.0],
        wind_speed        = 3.0,
        time_values       = [DateTime(2001, 1, 15, 12), DateTime(2001, 2, 15, 12), DateTime(2001, 3, 15, 12)],
    )
    return model, forcing, grid
end

function _write_sample_forcing_file(path::AbstractString)
    nx = 2; ny = 2; ntime = 2
    ds = NCDataset(path, "c")
    try
        defDim(ds, "x", nx)
        defDim(ds, "y", ny)
        defDim(ds, "time", ntime)

        defVar(ds, "x", Float64, ("x",))[:] = [0.0, 10_000.0]
        defVar(ds, "y", Float64, ("y",))[:] = [0.0, 10_000.0]
        defVar(ds, "MSK", Float64, ("x", "y"))[:, :] = [100.0 0.0; 75.0 25.0]

        time = defVar(ds, "time", Float64, ("time",))
        time.attrib["units"] = "days since 2001-01-01 12:00:00"
        time[:] = [0.0, 1.0]

        tt = fill(-15.0, ntime, ny, nx); tt[2, :, :] .= -14.0
        defVar(ds, "TT", Float64, ("time", "y", "x"))[:, :, :] = tt
        defVar(ds, "SF", Float64, ("time", "y", "x"))[:, :, :] = fill(1.0, ntime, ny, nx)
        defVar(ds, "RF", Float64, ("time", "y", "x"))[:, :, :] = fill(0.0, ntime, ny, nx)
        defVar(ds, "SWD", Float64, ("time", "y", "x"))[:, :, :] = fill(100.0, ntime, ny, nx)
        defVar(ds, "LWD", Float64, ("time", "y", "x"))[:, :, :] = fill(250.0, ntime, ny, nx)
        defVar(ds, "SHF", Float64, ("time", "y", "x"))[:, :, :] = fill(0.0, ntime, ny, nx)
        defVar(ds, "LHF", Float64, ("time", "y", "x"))[:, :, :] = fill(0.0, ntime, ny, nx)
    finally
        close(ds)
    end
    return path
end

function _write_sample_forcing_file_4d(path::AbstractString)
    nx = 2; ny = 2; ntime = 2; nlayer = 1
    ds = NCDataset(path, "c")
    try
        defDim(ds, "TIME", ntime)
        defDim(ds, "layer", nlayer)
        defDim(ds, "y", ny)
        defDim(ds, "x", nx)

        defVar(ds, "x", Float64, ("x",))[:] = [0.0, 10_000.0]
        defVar(ds, "y", Float64, ("y",))[:] = [0.0, 10_000.0]
        defVar(ds, "TIME", Float64, ("TIME",))[:] = [0.0, 1.0]
        defVar(ds, "YYYY", Float64, ("TIME",))[:] = [2001.0, 2001.0]
        defVar(ds, "MM", Float64, ("TIME",))[:] = [1.0, 1.0]
        defVar(ds, "DD", Float64, ("TIME",))[:] = [1.0, 2.0]
        defVar(ds, "HH", Float64, ("TIME",))[:] = [12.0, 12.0]

        tt = fill(-15.0, nx, ny, nlayer, ntime); tt[:, :, :, 2] .= -14.0
        defVar(ds, "TT", Float32, ("x", "y", "layer", "TIME"))[:, :, :, :] = tt
        defVar(ds, "SF", Float32, ("x", "y", "TIME"))[:, :, :] = fill(1.0, nx, ny, ntime)
        defVar(ds, "RF", Float32, ("x", "y", "TIME"))[:, :, :] = fill(0.0, nx, ny, ntime)
        defVar(ds, "SWD", Float32, ("x", "y", "TIME"))[:, :, :] = fill(100.0, nx, ny, ntime)
        defVar(ds, "LWD", Float32, ("x", "y", "TIME"))[:, :, :] = fill(250.0, nx, ny, ntime)
    finally
        close(ds)
    end
    return path
end

@testset "Run API" begin
    @testset "BESSIModel scheme types" begin
        grid = SnowpackGrid(1)
        m1 = BESSIModel(grid; albedo=:dynamic)
        m2 = BESSIModel(grid; albedo=:constant)
        m3 = BESSIModel(grid; albedo=:aging)
        m4 = BESSIModel(grid; densification=:htessel)
        m5 = BESSIModel(grid; fresh_snow_density=:parameterized)
        @test m1.c.albedo_scheme == Chion.ALBEDO_DYNAMIC
        @test m2.c.albedo_scheme == Chion.ALBEDO_CONSTANT
        @test m3.c.albedo_scheme == Chion.ALBEDO_AGING
        @test m3.c.alpha_dry == 0.81
        @test m3.c.alpha_wet == 0.70
        @test m3.c.aging_cold_timescale_days == 20.0
        @test m3.c.aging_melting_timescale_days == 5.0
        @test m4.c.low_density_densification == Chion.LOW_DENSIFICATION_HTESSEL
        @test m5.c.fresh_snow_density_scheme == Chion.FRESH_SNOW_DENSITY_PARAMETERIZED

        err = _captured_exception() do
            BESSIModel(grid; albedo=:aging, alpha_dry=0.5, alpha_wet=0.6)
        end
        @test err isa Exception
        @test occursin("alpha_wet", sprint(showerror, err))
    end

    @testset "aging albedo follows snowfall age" begin
        model = BESSIModel(SnowpackGrid(1); albedo=:aging)
        state = BESSIState(model)
        state.N[1] = 1
        state.mass[1, 1] = 300.0
        state.temperature[1, 1] = model.c.T0 - 1.0

        @test state.albedo[1] == model.c.alpha_dry
        @test state.snow_age_days[1] == 0.0

        Chion.update_surface_albedo!(state, 1, 1.0)
        @test state.snow_age_days[1] == 1.0
        @test state.albedo[1] ≈ model.c.alpha_wet +
                                 (model.c.alpha_dry - model.c.alpha_wet) * exp(-1 / 20) atol=1e-12

        state.temperature[1, 1] = model.c.T0
        Chion.update_surface_albedo!(state, 1, 1.0)
        @test state.snow_age_days[1] == 2.0
        @test state.albedo[1] ≈ model.c.alpha_wet +
                                 (model.c.alpha_dry - model.c.alpha_wet) *
                                 exp(-1 / 20 - 1 / 5) atol=1e-12

        Chion._update_aging_surface_albedo_arrays!(
            state.N,
            state.mass,
            state.temperature,
            state.albedo,
            state.snow_age_days,
            1,
            state.c,
            1.0 / state.c.seconds_per_day,
            1.0,
        )
        @test state.snow_age_days[1] == 0.0
        @test state.albedo[1] == model.c.alpha_dry

        state.N[1] = 0
        Chion.update_surface_albedo!(state, 1, 1.0)
        @test state.snow_age_days[1] == 0.0
        @test state.albedo[1] == model.c.alpha_ice
    end

    @testset "SnowpackForcing conversions and validation" begin
        forcing = SnowpackForcing(
            dt_days           = [1.0, 2.0],
            ncol              = 2,
            air_temperature_c = [-10.0, -5.0],
            snowfall_mm_day   = 1.5,
            rainfall_mm_day   = [0.0, 2.0],
            shortwave_down    = [100.0 120.0; 140.0 160.0],
            wind_speed        = 3.5,
            q_lw_down         = [250.0, 255.0],
            has_q_lw_down     = [true, false],
        )
        @test size(forcing.air_temperature) == (2, 2)
        @test forcing.air_temperature[1, 1] ≈ 263.15 atol=1e-8
        @test forcing.air_temperature[2, 2] ≈ 268.15 atol=1e-8
        @test forcing.snowfall_rate[2, 1]   ≈ 1.5 / 86_400.0 atol=1e-12
        @test forcing.rainfall_rate[1, 2]   ≈ 2.0 / 86_400.0 atol=1e-12
        @test forcing.shortwave_down == [100.0 120.0; 140.0 160.0]
        @test forcing.wind_speed == fill(3.5, 2, 2)
        @test forcing.has_q_lw_down == [true false; true false]
        @test length(forcing.time_values) == 2

        err = _captured_exception() do
            SnowpackForcing(
                dt_days=[1.0, 1.0], ncol=2,
                air_temperature_c=[-10.0],
                snowfall_mm_day=0.0, rainfall_mm_day=0.0, shortwave_down=0.0,
            )
        end
        @test err isa Exception
        @test occursin("air_temperature_c", sprint(showerror, err))
    end

    @testset "SnowpackForcing keeps invariant fields compact" begin
        forcing = SnowpackForcing(
            dt_days=[1.0, 1.0],
            ncol=3,
            air_temperature_c=[-10.0, -9.0],
            snowfall_mm_day=1.0,
            rainfall_mm_day=0.0,
            shortwave_down=100.0,
            latitude_deg=[60.0, 70.0, 80.0],
            surface_height=[0.0, 500.0, 1000.0],
        )

        @test forcing.q_lw_down isa Chion.ConstantForcingMatrix
        @test forcing.has_q_lw_down isa Chion.ConstantForcingMatrix
        @test forcing.air_temperature isa Chion.TimeForcingMatrix
        @test forcing.shortwave_down isa Chion.ConstantForcingMatrix
        @test forcing.latitude_deg isa Chion.ColumnForcingMatrix
        @test forcing.surface_height isa Chion.ColumnForcingMatrix
        @test Matrix(forcing.q_lw_down) == zeros(3, 2)
        @test Matrix(forcing.latitude_deg) == [60.0 60.0; 70.0 70.0; 80.0 80.0]
    end

    @testset "transposed layer storage preserves logical indexing" begin
        logical = reshape(collect(1.0:12.0), 3, 4)
        storage = Chion.TransposedLayerMatrix(permutedims(logical, (2, 1)))
        @test size(storage) == size(logical)
        @test Array(storage) == logical
        storage[2, 3] = -1.0
        @test storage.parent[3, 2] == -1.0

        scratch = similar(storage, Float64, 3, 4)
        @test scratch isa Chion.TransposedLayerMatrix
        @test size(scratch.parent) == (4, 3)
    end

    @testset "Simulation saves an exact state field" begin
        mktempdir() do dir
            model, forcing, grid = _sample_model_forcing_grid()
            sim = Simulation(model; forcing=forcing,
                netcdf_variables=:thickness,
                write_netcdf=true,
                output_dir=dir,
                netcdf_path=joinpath(dir, "save_exact.nc"),
                backend=:cpu,
                years=1,
                history_year_stride=1,
            )
            result = run!(sim)
            @test result.status == :complete
            @test isfile(result.netcdf_path)
            ds = NCDataset(result.netcdf_path)
            @test haskey(ds, "thickness")
            @test !haskey(ds, "cycle")
            @test !haskey(ds, "month_cycle")
            @test !haskey(ds, "step_cycle")
            @test !haskey(ds, "history_mean_thickness")
            @test size(ds["thickness"]) == (3, 2, 2)
            @test ds["t"].attrib["calendar"] == "proleptic_gregorian"
            @test ds["t"][:] == forcing.time_values
            @test ds.attrib["records_written"] == "3"
            @test !haskey(ds.attrib, "cycles_completed")
            close(ds)
        end
    end

    @testset "Simulation saves all default state fields" begin
        mktempdir() do dir
            model, forcing, grid = _sample_model_forcing_grid()
            sim = Simulation(model; forcing=forcing,
                netcdf_variables=:all,
                write_netcdf=true,
                output_dir=dir,
                netcdf_path=joinpath(dir, "save_group.nc"),
                backend=:threads,
                years=1,
                history_year_stride=1,
            )
            result = run!(sim)
            @test result.status == :complete
            ds = NCDataset(result.netcdf_path)
            @test haskey(ds, "thickness")
            @test haskey(ds, "mass_base")
            @test !haskey(ds, "history_mean_thickness")
            close(ds)
        end
    end

    @testset "Simulation saves mixed fields" begin
        mktempdir() do dir
            model, forcing, grid = _sample_model_forcing_grid()
            sim = Simulation(model; forcing=forcing,
                netcdf_variables=[:thickness, :mass_base],
                write_netcdf=true,
                output_dir=dir,
                netcdf_path=joinpath(dir, "save_mixed.nc"),
                backend=:threads,
                years=1,
                history_year_stride=1,
            )
            result = run!(sim)
            @test result.status == :complete
            ds = NCDataset(result.netcdf_path)
            @test haskey(ds, "thickness")
            @test haskey(ds, "mass_base")
            @test !haskey(ds, "history_mean_thickness")
            close(ds)
        end
    end

    @testset "Direct monthly variable selection records values" begin
        mktempdir() do dir
            model, forcing, _ = _sample_model_forcing_grid()
            path = joinpath(dir, "monthly_lhf.nc")
            result = run!(Simulation(model;
                forcing=forcing,
                netcdf_variables=:latent_heat_flux,
                write_netcdf=true,
                netcdf_path=path,
                backend=:threads,
                years=1,
            ); io=devnull)
            @test result.status == :complete
            ds = NCDataset(path)
            @test ds.attrib["records_written"] == "3"
            @test all(isfinite, Float64.(ds["latent_heat_flux"][1:3, :, :]))
            close(ds)
        end
    end

    @testset "Simulation can skip NetCDF entirely" begin
        model, forcing, _ = _sample_model_forcing_grid()
        sim = Simulation(model; forcing=forcing,
            backend=:threads,
            years=1,
        )
        result = run!(sim)
        @test result.status == :complete
        @test result.netcdf_path == ""
    end

    @testset "BESSI simulation can skip annual metrics" begin
        model, forcing, _ = _sample_model_forcing_grid()
        result = run!(Simulation(model; forcing=forcing,
            backend=:threads,
            years=2,
            write_netcdf=false,
            compute_year_metrics=false,
        ); io=devnull)
        timing_keys = [row.key for row in first(timing_rows(result.timings))]
        @test result.status == :complete
        @test isempty(result.history)
        @test !(:year_metrics in timing_keys)
        @test !(:summarize_columns_year in timing_keys)
    end

    @testset "BESSI KA interval step matches serial column stepping" begin
        model, forcing, _ = _sample_model_forcing_grid()
        serial = BESSIState(model)
        stepped = BESSIState(model)
        serial_workspace = Chion.ColumnarStepWorkspace(serial)
        stepped_workspace = Chion.ColumnarStepWorkspace(stepped)
        kwargs = Chion._bessi_step_kwargs(model)
        config = Chion._step_config_from_keywords(; kwargs...)

        for time_index in eachindex(forcing.time_values)
            @inbounds for idx in 1:serial.ncol
                Chion.column_step!(
                    serial,
                    idx,
                    Chion._step_forcing_at(forcing, idx, Int(time_index)),
                    config,
                    serial_workspace,
                )
            end
        end

        step!(
            stepped,
            forcing,
            stepped_workspace,
            1:stepped.ncol;
            kwargs...,
        )

        @test stepped.N == serial.N
        @test stepped.mass ≈ serial.mass
        @test stepped.mass_w ≈ serial.mass_w
        @test stepped.density ≈ serial.density
        @test stepped.temperature ≈ serial.temperature
        @test stepped.smb_ice ≈ serial.smb_ice
        @test stepped.runoff ≈ serial.runoff
        @test stepped.Tsrf ≈ serial.Tsrf
        @test stepped.albedo ≈ serial.albedo
        @test stepped.snow_age_days ≈ serial.snow_age_days
    end

    @testset "Simulation owns reference and current state" begin
        model, forcing, _ = _sample_model_forcing_grid()
        sim = Simulation(model; forcing=forcing, years=1)
        @test sim.ref !== sim.now
        @test all(sim.ref.mass .== 0.0)
        run!(sim; io=devnull)
        @test all(sim.ref.mass .== 0.0)
        @test sum(sim.now.mass) > 0.0
    end

    @testset "initialized integrator matches run wrapper" begin
        model_a, forcing, _ = _sample_model_forcing_grid()
        sim_run = Simulation(model_a; forcing=forcing, years=1)
        result_run = run!(sim_run; io=devnull)

        model_b = BESSIModel(model_a.grid; Ntot=4)
        sim_manual = Simulation(model_b; forcing=forcing, years=1)
        integrator = init_integrator(sim_manual; io=devnull)
        run!(integrator)
        result_manual = finalize!(integrator)

        @test result_run.status == result_manual.status == :complete
        @test result_run.years_completed == result_manual.years_completed == 1
        @test sim_run.now.mass ≈ sim_manual.now.mass
        @test sim_run.now.smb_ice ≈ sim_manual.now.smb_ice
    end

    @testset "yearly_step! returns weighted annual coupling fields" begin
        model, forcing, grid = _sample_model_forcing_grid()

        expected_sim = Simulation(model; forcing=forcing, years=1)
        expected_integrator = init_integrator(expected_sim; io=devnull)
        expected_smb_before = copy(expected_sim.now.smb_ice)
        expected_temperature_sum = zeros(Float64, grid.ncol)
        for time_index in eachindex(forcing.dt_days)
            step!(expected_integrator)
            @. expected_temperature_sum +=
                expected_integrator.model_runtime.data.state.Tsrf * forcing.dt_days[time_index]
        end

        coupled_sim = Simulation(BESSIModel(grid; Ntot=4); forcing=forcing, years=1)
        coupled_integrator = init_integrator(coupled_sim; io=devnull)
        coupling = yearly_step!(coupled_integrator)

        @test coupling.mean_T_srf_K ≈ expected_temperature_sum ./ sum(forcing.dt_days)
        @test coupling.ice_sheet_net_forcing_yearly ≈
            expected_sim.now.smb_ice .- expected_smb_before
        @test coupling.mean_T_srf_K_grid ≈
            Chion.scatter_to_grid(coupling.mean_T_srf_K, grid.js, grid.is, size(grid.mask))
    end

    @testset "external forcing step matches scheduled forcing" begin
        grid = SnowpackGrid(1)
        forcing = SnowpackForcing(
            dt_days=[1.0],
            air_temperature_c=[-8.0],
            snowfall_mm_day=[2.0],
            rainfall_mm_day=[0.0],
            shortwave_down=[100.0],
        )
        scheduled = Simulation(BESSIModel(grid; Ntot=4); forcing=forcing, years=1)
        run!(scheduled; io=devnull)

        external = Simulation(BESSIModel(grid; Ntot=4); forcing=forcing, years=1)
        integrator = init_integrator(external; io=devnull)
        step!(integrator, 1.0, true)
        result = finalize!(integrator)

        @test result.status == :complete
        @test external.now.mass ≈ scheduled.now.mass
        @test external.now.smb_ice ≈ scheduled.now.smb_ice
    end

    @testset "integrator exposes mutable scheduled forcing" begin
        grid = SnowpackGrid(2)
        forcing = SnowpackForcing(
            dt_days=[1.0, 1.0],
            air_temperature_c=[-8.0 -7.0; -10.0 -9.0],
            snowfall_mm_day=1.0,
            rainfall_mm_day=0.0,
            shortwave_down=100.0,
        )
        sim = Simulation(BESSIModel(grid; Ntot=4); forcing=forcing, years=1)
        integrator = init_integrator(sim; io=devnull)
        surface_height = [0.0, 1000.0]

        integrator.sim.forcing.surface_height .= surface_height
        update_air_pressure!(integrator.sim.forcing)
        sync_forcing!(integrator)

        expected = Chion.air_pressure_from_surface_height(
            surface_height,
            forcing.air_temperature;
            dt_days=forcing.dt_days,
            time_values=forcing.time_values,
            temperature_mode=:instantaneous,
        )
        @test integrator.sim.forcing.air_pressure ≈ expected
        @test integrator.model_runtime.data.step_fields.air_pressure ≈ expected
    end

    @testset "non-spatial grids only require coordinates for NetCDF" begin
        grid = SnowpackGrid(1)
        model = BESSIModel(grid; Ntot=4)
        forcing = SnowpackForcing(
            dt_days=[1.0, 1.0],
            air_temperature_c=[-10.0, -9.0],
            snowfall_mm_day=0.5,
            rainfall_mm_day=0.0,
            shortwave_down=[100.0, 120.0],
        )
        result = run!(Simulation(model; forcing=forcing, years=1))
        @test result.status == :complete
        @test result.netcdf_path == ""

        err = _captured_exception() do
            run!(Simulation(model; forcing=forcing, netcdf_variables=:thickness, write_netcdf=true, years=1))
        end
        @test err isa Exception
        @test occursin("spatial coordinates", sprint(showerror, err))
    end

    @testset "forcing-file loading and run!" begin
        mktempdir() do dir
            forcing_path = _write_sample_forcing_file(joinpath(dir, "forcing_sample.nc"))
            loaded = load_forcing_file(forcing_path)
            @test size(loaded.forcing.air_temperature) == (4, 2)
            @test size(loaded.forcing.snowfall_rate) == (4, 2)
            @test loaded.forcing.time_values ==
                  [DateTime(2001, 1, 1, 12), DateTime(2001, 1, 2, 12)]
            @test sort(collect(zip(loaded.grid.js, loaded.grid.is))) == [(1, 1), (1, 2), (2, 1), (2, 2)]

            masked = load_forcing_file(forcing_path; mask_name="MSK", mask_threshold=50.0)
            @test size(masked.forcing.air_temperature) == (2, 2)
            @test collect(zip(masked.grid.js, masked.grid.is)) == [(1, 1), (1, 2)]
            @test masked.grid.mask == [100.0 75.0; 0.0 25.0]

            sim = Simulation(BESSIModel(loaded.grid; Ntot=4);
                forcing=loaded.forcing,
                years=1,
                backend=:threads,
            )
            result = run!(sim)
            @test result.status == :complete
            @test result.netcdf_path == ""
        end
    end

    @testset "forcing-file loader accepts singleton-layer 4D variables" begin
        mktempdir() do dir
            forcing_path = _write_sample_forcing_file_4d(joinpath(dir, "forcing_sample_4d.nc"))
            loaded = load_forcing_file(forcing_path)
            @test size(loaded.forcing.air_temperature) == (4, 2)
            @test loaded.forcing.air_temperature[1, 1] ≈ 258.15 atol=1e-5
            @test loaded.forcing.air_temperature[1, 2] ≈ 259.15 atol=1e-5
            @test all(loaded.forcing.has_q_lw_down)
        end
    end

    @testset "PDDModel runs bulk positive degree days" begin
        grid = SnowpackGrid(1)
        f = SnowpackForcing(
            dt_days=[1.0, 1.0],
            air_temperature_c=[-5.0, 1.0],
            snowfall_mm_day=[2.0, 0.0],
            rainfall_mm_day=0.0,
            shortwave_down=150.0,
        )
        pdd = build_model(:pdd, grid; ddf_snow=3.0, ddf_ice=8.0, refreezing_fraction=0.0)
        pdd_sim = Simulation(pdd; forcing=f, years=1)
        result = run!(pdd_sim)
        @test result.status == :complete
        @test pdd_sim.now.snowpack_swe[1] ≈ 0.0 atol=1e-12
        @test pdd_sim.now.smb_ice[1] ≈ -(8.0 / 3.0) atol=1e-12
        @test pdd_sim.now.runoff[1] ≈ 2.0 + 8.0 / 3.0 atol=1e-12
        @test pdd_sim.now.pdd_sum[1] ≈ 1.0 atol=1e-12

        named = Simulation(:pdd, grid;
            model_kwargs=(ddf_snow=3.0, ddf_ice=8.0, refreezing_fraction=0.0),
            forcing=f,
            years=1,
        )
        @test run!(named).status == :complete

        mktempdir() do dir
            spatial_grid = SnowpackGrid(1;
                x=[0.0],
                y=[0.0],
                js=[1],
                is=[1],
                mask=ones(1, 1),
            )
            spatial_pdd = build_model(:pdd, spatial_grid; ddf_snow=3.0, ddf_ice=8.0, refreezing_fraction=0.0)
            result = run!(Simulation(spatial_pdd;
                forcing=f,
                years=1,
            ))
            @test result.status == :complete
            @test result.netcdf_path == ""
            @test result.years_completed == 1

            output_path = joinpath(dir, "pdd.nc")
            output_sim = Simulation(spatial_pdd;
                forcing=f,
                write_netcdf=true,
                netcdf_variables=:all,
                netcdf_path=output_path,
                years=1,
            )
            output_result = run!(output_sim)
            @test output_result.status == :complete
            @test output_result.netcdf_path == output_path

            ds = NCDataset(output_path)
            @test all(haskey(ds, name) for name in ("snowpack_swe", "smb_ice", "runoff", "pdd_sum"))
            @test !haskey(ds, "thickness")
            @test size(ds["smb_ice"]) == (2, 1, 1)
            @test ds.attrib["records_written"] == "2"
            @test ds["smb_ice"][2, 1, 1] ≈ output_sim.now.smb_ice[1] atol=1e-5
            @test ds["pdd_sum"][2, 1, 1] ≈ output_sim.now.pdd_sum[1] atol=1e-5
            close(ds)

            monthly_path = joinpath(dir, "pdd_monthly.nc")
            monthly_sim = Simulation(PDDModel(spatial_grid;
                    ddf_snow=3.0,
                    ddf_ice=8.0,
                    refreezing_fraction=0.0,
                );
                forcing=f,
                write_netcdf=true,
                netcdf_variables=:monthly,
                netcdf_path=monthly_path,
                years=1,
            )
            @test run!(monthly_sim).status == :complete
            ds = NCDataset(monthly_path)
            @test ds.attrib["records_written"] == "1"
            @test ds["smb_ice"][1, 1, 1] ≈ monthly_sim.now.smb_ice[1] atol=1e-5
            @test ds["runoff"][1, 1, 1] ≈ monthly_sim.now.runoff[1] atol=1e-5
            @test ds["pdd_sum"][1, 1, 1] ≈ monthly_sim.now.pdd_sum[1] atol=1e-5
            close(ds)

            err = _captured_exception() do
                init_integrator(Simulation(PDDModel(spatial_grid);
                    forcing=f,
                    write_netcdf=true,
                    netcdf_variables=:thickness,
                    netcdf_path=joinpath(dir, "unsupported.nc"),
                    years=1,
                ))
            end
            @test err isa Exception
            @test occursin("Unsupported NetCDF variables for PDDModel", sprint(showerror, err))
        end
    end

    @testset "ITMModel is not public API yet" begin
        grid = SnowpackGrid(1)
        @test !(Symbol("ITMModel") in Set(names(Chion)))
        @test _captured_exception(() -> build_model(:itm, grid)) isa Exception
    end

    @testset "legacy names are not exported or documented" begin
        exported_names = Set(names(Chion))
        for name in (:ForcingData, :SnowpackStepFields, :GridLayout, :LoadedProblem,
                :RunConfig, :SimulationOptions, :OutputOptions,
                :run_case, :prescribed_case, :synthetic_case,
                :CurrentState, :StochasticMonthlyPDD,
                :DynamicAlbedo, :ConstantAlbedo, :PrescribedAlbedo,
                :BESSIDensification, :HTESSELDensification,
                :ConstantFreshSnowDensity, :ParameterizedFreshSnowDensity)
            @test !(name in exported_names)
        end

        root = dirname(dirname(@__FILE__))
        docs_text = join(read.(filter(endswith(".md"), readdir(joinpath(root, "docs", "src"); join=true)), String), "\n")
        for token in ("ForcingData", "SnowpackStepFields", "GridLayout", "LoadedProblem",
                "RunConfig", "SimulationOptions", "OutputOptions",
                "run_case", "prescribed_case", "synthetic_case")
            @test !occursin(token, docs_text)
        end
    end
end
