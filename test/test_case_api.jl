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
        @test m1.domain.c.albedo_scheme == Chion.ALBEDO_DYNAMIC
        @test m2.domain.c.albedo_scheme == Chion.ALBEDO_CONSTANT
        @test m3.domain.c.low_density_densification == Chion.LOW_DENSIFICATION_HTESSEL
        @test m4.domain.c.fresh_snow_density_scheme == Chion.FRESH_SNOW_DENSITY_PARAMETERIZED
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
                cycles=1,
                history_stride=1,
            )
            result = run!(sim)
            @test result.status == :cycles
            @test isfile(result.netcdf_path)
            ds = NCDataset(result.netcdf_path)
            @test haskey(ds, "final_thickness")
            @test !haskey(ds, "history_mean_thickness")
            @test size(ds["final_thickness"]) == (2, 2)
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
                cycles=1,
                history_stride=1,
            )
            result = run!(sim)
            @test result.status == :cycles
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
                cycles=1,
                history_stride=1,
            )
            result = run!(sim)
            @test result.status == :cycles
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
            cycles=1,
        )
        result = run!(sim)
        @test result.status == :cycles
        @test result.netcdf_path == ""
        @test result.summary_path == ""
        @test result.history_csv_path == ""
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
        result = run!(Simulation(model; forcing=forcing, save=:none, cycles=1, write_outputs=false))
        @test result.status == :cycles
        @test result.netcdf_path == ""

        err = _captured_exception() do
            run!(Simulation(model; forcing=forcing, save=:final_thickness, cycles=1, write_outputs=false))
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
                options=SimulationOptions(cycles=1, backend=:threads),
            )
            result = run!(sim)
            @test result.status == :cycles
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

    @testset "PDDModel and ITMModel stubs error on run!" begin
        grid = SnowpackGrid(CPU(), 1)
        f = SnowpackForcing(
            dt_days=fill(1.0, 3),
            air_temperature_c=fill(-10.0, 3),
            snowfall_mm_day=1.0,
            rainfall_mm_day=0.0,
            shortwave_down=150.0,
        )
        pdd_sim = Simulation(PDDModel(grid); forcing=f)
        itm_sim = Simulation(ITMModel(grid); forcing=f)
        @test _captured_exception(() -> run!(pdd_sim)) isa Exception
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
