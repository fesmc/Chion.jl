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
        m3 = BESSIModel(grid; densification=:htessel)
        m4 = BESSIModel(grid; fresh_snow_density=:parameterized)
        @test m1.c.albedo_scheme == Chion.ALBEDO_DYNAMIC
        @test m2.c.albedo_scheme == Chion.ALBEDO_CONSTANT
        @test m3.c.low_density_densification == Chion.LOW_DENSIFICATION_HTESSEL
        @test m4.c.fresh_snow_density_scheme == Chion.FRESH_SNOW_DENSITY_PARAMETERIZED
    end

    @testset "SnowpackForcing conversions and validation" begin
        noleap_times = Chion._synthesized_time_values(fill(1.0, 365))
        @test first(noleap_times) == DateTime(2001, 1, 1, 12)
        @test last(noleap_times) == DateTime(2001, 12, 31, 12)
        @test all(value -> !(month(value) == 2 && day(value) == 29), noleap_times)

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
            @test haskey(ds, "t")
            @test !haskey(ds, "cycle")
            @test !haskey(ds, "month_cycle")
            @test !haskey(ds, "step_cycle")
            @test !haskey(ds, "history_mean_thickness")
            @test size(ds["thickness"]) == (2, 2, 3)
            @test ds["t"].var[:] ≈ Chion._netcdf_time_days.(forcing.time_values)
            @test ds["t"].attrib["standard_name"] == "time"
            @test ds["t"].attrib["calendar"] == "365_day"
            @test ds["t"].attrib["axis"] == "T"
            @test "t" in NCDatasets.unlimited(ds)
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
            @test ds["t"].var[:] ≈ Chion._netcdf_time_days.(forcing.time_values)
            @test all(isfinite, Float64.(ds["latent_heat_flux"][:, :, 1:3]))
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

    @testset "Bare-ice rainfall is routed directly to runoff" begin
        grid = SnowpackGrid(1)
        model = BESSIModel(grid; Ntot=4)
        rainfall_mm = 7.0
        neutral_longwave_down = model.c.σ * model.c.ϵ_snow * model.c.T0^4
        forcing = SnowpackForcing(
            dt_days=[1.0],
            air_temperature=[model.c.T0],
            snowfall_rate=[0.0],
            rainfall_rate=[rainfall_mm / model.c.seconds_per_day],
            shortwave_down=[0.0],
            q_lw_down=[neutral_longwave_down],
            q_sh=[0.0],
            q_lh=[0.0],
        )
        sim = Simulation(model; forcing=forcing, backend=:cpu, years=1)

        result = run!(sim; io=devnull)

        @test result.status == :complete
        @test sim.now.runoff[1] ≈ rainfall_mm atol=1e-12
        @test sim.now.smb_ice[1] ≈ 0.0 atol=1e-12
        @test sim.now.melt[1] ≈ 0.0 atol=1e-12
        @test sum(sim.now.mass) ≈ 0.0 atol=1e-12
        @test sum(sim.now.mass_w) ≈ 0.0 atol=1e-12
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
                years=2,
            )
            output_result = run!(output_sim)
            @test output_result.status == :complete
            @test output_result.netcdf_path == output_path

            ds = NCDataset(output_path)
            @test all(haskey(ds, name) for name in ("snowpack_swe", "smb_ice", "runoff", "pdd_sum"))
            @test !haskey(ds, "thickness")
            @test size(ds["smb_ice"]) == (1, 1, 4)
            expected_times = vcat(f.time_values, f.time_values .+ Year(1))
            @test ds["t"].var[:] ≈ Chion._netcdf_time_days.(expected_times)
            @test "t" in NCDatasets.unlimited(ds)
            @test ds.attrib["records_written"] == "4"
            @test ds["smb_ice"][1, 1, 4] ≈ output_sim.now.smb_ice[1] atol=1e-5
            @test ds["pdd_sum"][1, 1, 4] ≈ output_sim.now.pdd_sum[1] atol=1e-5
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
            @test size(ds["smb_ice"]) == (1, 1, 1)
            @test ds["t"].var[:] ≈ Chion._netcdf_time_days.([f.time_values[end]])
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

    @testset "PDD uses the capped ice-facing Fortran budget" begin
        grid = SnowpackGrid(1)
        model = PDDModel(
            grid;
            ddf_snow=3.0,
            ddf_ice=8.0,
            refreezing_fraction=0.6,
            H_snow_max=100.0,
            pdd_method=:simple,
        )

        cold_accumulation = SnowpackForcing(
            dt_days=[1.0],
            air_temperature=[model.c.T0 - 10.0],
            snowfall_rate=[150.0 / model.c.seconds_per_day],
            rainfall_rate=[20.0 / model.c.seconds_per_day],
            shortwave_down=[0.0],
        )
        accumulated = PDDState(model)
        Chion.pdd_step!(model, accumulated, cold_accumulation)

        @test accumulated.snowpack_swe[1] ≈ 100.0 atol=1e-12
        @test accumulated.smb_ice[1] ≈ 50.0 atol=1e-12
        @test accumulated.runoff[1] ≈ 20.0 atol=1e-12
        @test accumulated.snowpack_swe[1] +
              accumulated.smb_ice[1] +
              accumulated.runoff[1] ≈ 170.0 atol=1e-12

        capacity_limited = PDDState(model)
        capacity_limited.snowpack_swe[1] = 10.0
        melt_forcing = SnowpackForcing(
            dt_days=[1.0],
            air_temperature=[model.c.T0 + 3.0],
            snowfall_rate=[0.0],
            rainfall_rate=[0.0],
            shortwave_down=[0.0],
        )
        Chion.pdd_step!(model, capacity_limited, melt_forcing)

        @test capacity_limited.snowpack_swe[1] ≈ 1.0 atol=1e-12
        @test capacity_limited.smb_ice[1] ≈ 0.6 atol=1e-12
        @test capacity_limited.runoff[1] ≈ 8.4 atol=1e-12
        @test capacity_limited.snowpack_swe[1] +
              capacity_limited.smb_ice[1] +
              capacity_limited.runoff[1] ≈ 10.0 atol=1e-12

        extreme_melt = SnowpackForcing(
            dt_days=[1.0],
            air_temperature=[model.c.T0 + 1000.0],
            snowfall_rate=[0.0],
            rainfall_rate=[0.0],
            shortwave_down=[0.0],
        )
        Chion.pdd_step!(model, capacity_limited, extreme_melt)
        @test capacity_limited.snowpack_swe[1] == 0.0
        @test capacity_limited.smb_ice[1] < 0.0
    end

    @testset "PDD method and physical constants are explicit" begin
        grid = SnowpackGrid(1)
        simple = PDDModel(grid; pdd_method=:simple)
        pism = PDDModel(grid; pdd_method=:pism)
        cold = SnowpackForcing(
            dt_days=[1.0],
            air_temperature_c=[-5.0],
            snowfall_mm_day=[0.0],
            rainfall_mm_day=[0.0],
            shortwave_down=[0.0],
        )
        simple_state = PDDState(simple)
        pism_state = PDDState(pism)
        Chion.pdd_step!(simple, simple_state, cold)
        Chion.pdd_step!(pism, pism_state, cold)

        @test simple_state.pdd_sum[1] == 0.0
        @test pism_state.pdd_sum[1] ≈ 0.4165773529384319 atol=1e-12

        custom_constants = PDDModel(
            grid;
            pdd_method=:simple,
            T0=270.0,
            seconds_per_day=100.0,
        )
        custom_state = PDDState(custom_constants)
        custom_forcing = SnowpackForcing(
            dt_days=[1.0],
            air_temperature=[271.0],
            snowfall_rate=[1.0],
            rainfall_rate=[0.0],
            shortwave_down=[0.0],
        )
        Chion.pdd_step!(custom_constants, custom_state, custom_forcing)
        @test custom_state.pdd_sum[1] ≈ 1.0 atol=1e-12
        @test custom_state.snowpack_swe[1] ≈ 97.0 atol=1e-12
        @test custom_state.smb_ice[1] ≈ 3.0 atol=1e-12

        @test _captured_exception(() -> PDDModel(grid; pdd_method=:unknown)) isa Exception
        @test _captured_exception(() -> PDDModel(grid; H_snow_max=0.0)) isa Exception
    end

    @testset "PDD honors active columns through the runtime path" begin
        grid = SnowpackGrid(2)
        model = PDDModel(grid)
        forcing = SnowpackForcing(
            dt_days=[1.0],
            ncol=2,
            air_temperature_c=-10.0,
            snowfall_mm_day=5.0,
            rainfall_mm_day=0.0,
            shortwave_down=0.0,
        )
        sim = Simulation(model; forcing=forcing, backend=:cpu, years=1)
        integrator = init_integrator(sim; io=devnull)
        set_active_mask!(integrator, [true, false])
        step!(integrator)

        @test sim.now.snowpack_swe ≈ [5.0, 0.0] atol=1e-12
        @test sim.now.smb_ice == [0.0, 0.0]
        @test sim.now.runoff == [0.0, 0.0]
        @test sim.now.pdd_sum == [0.0, 0.0]
    end

    if Chion.cuda_available()
        @testset "PDD GPU path matches CPU" begin
            grid = SnowpackGrid(2)
            model = PDDModel(grid; H_snow_max=8.0, pdd_method=:pism)
            forcing = SnowpackForcing(
                dt_days=[1.0, 1.0],
                ncol=2,
                air_temperature_c=[-5.0, 2.0],
                snowfall_mm_day=5.0,
                rainfall_mm_day=1.0,
                shortwave_down=0.0,
            )
            cpu = Simulation(model; forcing=forcing, backend=:cpu, years=1)
            gpu = Simulation(model; forcing=forcing, backend=:gpu, years=1)

            run!(cpu; io=devnull)
            run!(gpu; io=devnull)

            @test gpu.now.snowpack_swe ≈ cpu.now.snowpack_swe rtol=1e-6
            @test gpu.now.smb_ice ≈ cpu.now.smb_ice rtol=1e-6
            @test gpu.now.runoff ≈ cpu.now.runoff rtol=1e-6
            @test gpu.now.pdd_sum ≈ cpu.now.pdd_sum rtol=1e-6
        end
    end

    @testset "ITMModel matches the Fortran bulk budget" begin
        grid = SnowpackGrid(2)
        model = ITMModel(grid)
        forcing = SnowpackForcing(
            dt_days=[2.0],
            ncol=2,
            air_temperature=reshape([274.0, 274.0], 2, 1),
            snowfall_rate=reshape([0.0, 0.0], 2, 1),
            rainfall_rate=reshape([0.0, 0.0], 2, 1),
            shortwave_down=reshape([0.0, 0.0], 2, 1),
            q_sw_net=reshape([400.0, 400.0], 2, 1),
            latitude_deg=reshape([72.0, 65.0], 2, 1),
            surface_height=reshape([1500.0, 500.0], 2, 1),
            ice_thickness=reshape([1000.0, 0.0], 2, 1),
            annual_pdd=reshape([200.0, 800.0], 2, 1),
        )
        sim = Simulation(:itm, grid; forcing=forcing, years=1)
        @test sim.model isa ITMModel
        @test sim.now.H_snow == fill(model.H_snow_max, 2)

        result = run!(sim; io=devnull)
        @test result.status == :complete
        @test all(sim.now.melt .> 0.0)
        @test all(sim.now.H_snow .<= model.H_snow_max)
        @test sim.now.alb_s[1] == model.alb_snow_wet
        @test sim.now.smb_ice[1] ≈ sim.now.smbi[1] * 2.0 atol=1e-12

        no_qsw = SnowpackForcing(
            dt_days=[1.0], air_temperature=[274.0], snowfall_rate=[0.0], rainfall_rate=[0.0],
            shortwave_down=[0.0], latitude_deg=[72.0], surface_height=[1500.0],
            ice_thickness=[1000.0], annual_pdd=[200.0],
        )
        no_qsw_sim = Simulation(ITMModel(SnowpackGrid(1)); forcing=no_qsw, years=1)
        @test run!(no_qsw_sim; io=devnull).status == :complete
        @test no_qsw_sim.now.melt[1] == 0.0

        incomplete = SnowpackForcing(
            dt_days=[1.0], air_temperature=[274.0], snowfall_rate=[0.0], rainfall_rate=[0.0],
            shortwave_down=[0.0], latitude_deg=[72.0], surface_height=[1500.0],
        )
        @test _captured_exception(() -> init_integrator(Simulation(ITMModel(SnowpackGrid(1)); forcing=incomplete))) isa Exception

        albedo = Chion._itm_surface_albedo(model, 1500.0, 1000.0, 0.0, 200.0)
        @test albedo == model.alb_ice
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
