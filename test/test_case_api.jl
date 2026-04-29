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
    grid = SnowpackGrid(CPU(), 4;
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
        grid = SnowpackGrid(CPU(), 1)
        m1 = BESSIModel(grid; albedo=DynamicAlbedo())
        m2 = BESSIModel(grid; albedo=ConstantAlbedo())
        m3 = BESSIModel(grid; densification=HTESSELDensification())
        m4 = BESSIModel(grid; fresh_snow_density=ParameterizedFreshSnowDensity())
        @test m1.c.albedo_scheme == Chion.ALBEDO_DYNAMIC
        @test m2.c.albedo_scheme == Chion.ALBEDO_CONSTANT
        @test m3.c.low_density_densification == Chion.LOW_DENSIFICATION_HTESSEL
        @test m4.c.fresh_snow_density_scheme == Chion.FRESH_SNOW_DENSITY_PARAMETERIZED
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

    @testset "Simulation saves an exact symbol" begin
        mktempdir() do dir
            model, forcing, grid = _sample_model_forcing_grid()
            sim = Simulation(model; forcing=forcing,
                save=:final_thickness,
                output_dir=dir,
                netcdf_path=joinpath(dir, "save_exact.nc"),
                backend=:cpu,
                write_outputs=false,
                years=1,
                history_year_stride=1,
            )
            result = run!(sim)
            @test result.status == :complete
            @test isfile(result.netcdf_path)
            ds = NCDataset(result.netcdf_path)
            @test haskey(ds, "final_thickness")
            @test haskey(ds, "year")
            @test haskey(ds, "month")
            @test haskey(ds, "step_year")
            @test !haskey(ds, "cycle")
            @test !haskey(ds, "month_cycle")
            @test !haskey(ds, "step_cycle")
            @test !haskey(ds, "history_mean_thickness")
            @test size(ds["final_thickness"]) == (2, 2)
            @test haskey(ds.attrib, "years_completed")
            @test !haskey(ds.attrib, "cycles_completed")
            close(ds)
        end
    end

    @testset "Simulation saves a group" begin
        mktempdir() do dir
            model, forcing, grid = _sample_model_forcing_grid()
            sim = Simulation(model; forcing=forcing,
                save=:history,
                output_dir=dir,
                netcdf_path=joinpath(dir, "save_group.nc"),
                backend=:threads,
                write_outputs=false,
                years=1,
                history_year_stride=1,
            )
            result = run!(sim)
            @test result.status == :complete
            ds = NCDataset(result.netcdf_path)
            @test haskey(ds, "history_mean_thickness")
            @test haskey(ds, "history_mean_base_mass")
            @test !haskey(ds, "final_thickness")
            close(ds)
        end
    end

    @testset "Simulation saves mixed fields" begin
        mktempdir() do dir
            model, forcing, grid = _sample_model_forcing_grid()
            sim = Simulation(model; forcing=forcing,
                save=[:final_thickness, :history],
                output_dir=dir,
                netcdf_path=joinpath(dir, "save_mixed.nc"),
                backend=:threads,
                write_outputs=false,
                years=1,
                history_year_stride=1,
            )
            result = run!(sim)
            @test result.status == :complete
            ds = NCDataset(result.netcdf_path)
            @test haskey(ds, "final_thickness")
            @test haskey(ds, "history_mean_thickness")
            close(ds)
        end
    end

    @testset "Simulation can skip NetCDF entirely" begin
        model, forcing, _ = _sample_model_forcing_grid()
        sim = Simulation(model; forcing=forcing,
            save=:none,
            backend=:threads,
            write_outputs=false,
            years=1,
        )
        result = run!(sim)
        @test result.status == :complete
        @test result.netcdf_path == ""
        @test result.summary_path == ""
        @test result.history_csv_path == ""
    end

    @testset "Simulation owns reference and current state" begin
        model, forcing, _ = _sample_model_forcing_grid()
        sim = Simulation(model; forcing=forcing, save=:none, years=1, write_outputs=false)
        @test sim.ref !== sim.now
        @test all(sim.ref.domain.mass .== 0.0)
        run!(sim; io=devnull)
        @test all(sim.ref.domain.mass .== 0.0)
        @test sum(sim.now.domain.mass) > 0.0
    end

    @testset "initialized integrator matches run wrapper" begin
        model_a, forcing, _ = _sample_model_forcing_grid()
        sim_run = Simulation(model_a; forcing=forcing, save=:none, years=1, write_outputs=false)
        result_run = run!(sim_run; io=devnull)

        model_b = BESSIModel(model_a.grid; Ntot=4)
        sim_manual = Simulation(model_b; forcing=forcing, save=:none, years=1, write_outputs=false)
        integrator = init_integrator(sim_manual; io=devnull)
        run!(integrator)
        result_manual = finalize!(integrator)

        @test result_run.status == result_manual.status == :complete
        @test result_run.years_completed == result_manual.years_completed == 1
        @test sim_run.now.domain.mass ≈ sim_manual.now.domain.mass
        @test sim_run.now.domain.smb_ice ≈ sim_manual.now.domain.smb_ice
    end

    @testset "external forcing step matches scheduled forcing" begin
        grid = SnowpackGrid(CPU(), 1)
        forcing = SnowpackForcing(
            dt_days=[1.0],
            air_temperature_c=[-8.0],
            snowfall_mm_day=[2.0],
            rainfall_mm_day=[0.0],
            shortwave_down=[100.0],
        )
        scheduled = Simulation(BESSIModel(grid; Ntot=4); forcing=forcing, save=:none, years=1, write_outputs=false)
        run!(scheduled; io=devnull)

        external = Simulation(BESSIModel(grid; Ntot=4); forcing=forcing, save=:none, years=1, write_outputs=false)
        integrator = init_integrator(external; io=devnull)
        set_forcing!(
            integrator;
            air_temperature_c=-8.0,
            snowfall_mm_day=2.0,
            rainfall_mm_day=0.0,
            shortwave_down=100.0,
        )
        step!(integrator, 1.0, true)
        result = finalize!(integrator)

        @test result.status == :complete
        @test external.now.domain.mass ≈ scheduled.now.domain.mass
        @test external.now.domain.smb_ice ≈ scheduled.now.domain.smb_ice
    end

    @testset "non-spatial grids only require coordinates for NetCDF" begin
        grid = SnowpackGrid(CPU(), 1)
        model = BESSIModel(grid; Ntot=4)
        forcing = SnowpackForcing(
            dt_days=[1.0, 1.0],
            air_temperature_c=[-10.0, -9.0],
            snowfall_mm_day=0.5,
            rainfall_mm_day=0.0,
            shortwave_down=[100.0, 120.0],
        )
        result = run!(Simulation(model; forcing=forcing, save=:none, years=1, write_outputs=false))
        @test result.status == :complete
        @test result.netcdf_path == ""

        err = _captured_exception() do
            run!(Simulation(model; forcing=forcing, save=:final_thickness, years=1, write_outputs=false))
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

            sim = Simulation(BESSIModel(loaded.grid; Ntot=4);
                forcing=loaded.forcing,
                output=OutputOptions(save=:none, write_outputs=false),
                options=SimulationOptions(years=1, backend=:threads),
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
        grid = SnowpackGrid(CPU(), 1)
        f = SnowpackForcing(
            dt_days=[1.0, 1.0],
            air_temperature_c=[-5.0, 1.0],
            snowfall_mm_day=[2.0, 0.0],
            rainfall_mm_day=0.0,
            shortwave_down=150.0,
        )
        pdd = build_model(:pdd, grid; ddf_snow=3.0, ddf_ice=8.0, refreezing_fraction=0.0)
        pdd_sim = Simulation(pdd; forcing=f, save=:none, years=1, write_outputs=false)
        result = run!(pdd_sim)
        @test result.status == :complete
        @test pdd_sim.now.snowpack_swe[1] ≈ 0.0 atol=1e-12
        @test pdd_sim.now.smb_ice[1] ≈ -(8.0 / 3.0) atol=1e-12
        @test pdd_sim.now.runoff[1] ≈ 2.0 + 8.0 / 3.0 atol=1e-12
        @test pdd_sim.now.pdd_sum[1] ≈ 1.0 atol=1e-12

        named = Simulation(:pdd, grid;
            model_kwargs=(ddf_snow=3.0, ddf_ice=8.0, refreezing_fraction=0.0),
            forcing=f,
            save=:none,
            years=1,
            write_outputs=false,
        )
        @test run!(named).status == :complete

        mktempdir() do dir
            spatial_grid = SnowpackGrid(CPU(), 1;
                x=[0.0],
                y=[0.0],
                js=[1],
                is=[1],
                mask=ones(1, 1),
            )
            spatial_pdd = build_model(:pdd, spatial_grid; ddf_snow=3.0, ddf_ice=8.0, refreezing_fraction=0.0)
            result = run!(Simulation(spatial_pdd;
                forcing=f,
                save=[:final_wet_mass, :final_ice_sheet_smb, :history, :step_pdd],
                netcdf_path=joinpath(dir, "pdd.nc"),
                years=1,
                write_outputs=false,
            ))
            @test isfile(result.netcdf_path)
            ds = NCDataset(result.netcdf_path)
            @test haskey(ds, "final_wet_mass")
            @test haskey(ds, "final_ice_sheet_smb")
            @test haskey(ds, "history_mean_wet_mass")
            @test haskey(ds, "step_pdd")
            @test ds["final_wet_mass"][1, 1] ≈ 0.0 atol=1e-12
            @test ds["final_ice_sheet_smb"][1, 1] ≈ -(8.0 / 3.0) atol=1e-6
            @test ds["step_pdd"][1, 1, 1] ≈ 1.0 atol=1e-6
            close(ds)
        end
    end

    @testset "ITMModel stub errors on run!" begin
        grid = SnowpackGrid(CPU(), 1)
        f = SnowpackForcing(
            dt_days=fill(1.0, 3),
            air_temperature_c=fill(-10.0, 3),
            snowfall_mm_day=1.0,
            rainfall_mm_day=0.0,
            shortwave_down=150.0,
        )
        itm_sim = Simulation(ITMModel(grid); forcing=f)
        @test _captured_exception(() -> run!(itm_sim)) isa Exception
    end

    @testset "legacy names are not exported or documented" begin
        exported_names = Set(names(Chion))
        for name in (:ForcingData, :SnowpackStepFields, :GridLayout, :LoadedProblem,
                :RunConfig, :run_case, :prescribed_case, :synthetic_case)
            @test !(name in exported_names)
        end

        root = dirname(dirname(@__FILE__))
        docs_text = join(read.(filter(endswith(".md"), readdir(joinpath(root, "docs", "src"); join=true)), String), "\n")
        for token in ("ForcingData", "SnowpackStepFields", "GridLayout", "LoadedProblem",
                "RunConfig", "run_case", "prescribed_case", "synthetic_case")
            @test !occursin(token, docs_text)
        end
    end
end
