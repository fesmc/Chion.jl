using Test
using Dates
using HDF5
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

function _sample_domain_forcing_layout()
    physics = Chion.physics()
    domain = Chion.SnowpackDomain(c=physics, Ntot=4, ncol=4)
    forcing = Chion.ForcingData(
        dt_days=[1.0, 1.0, 1.0],
        ncol=4,
        air_temperature_c=[-15.0, -14.0, -13.0],
        snowfall_mm_day=1.0,
        rainfall_mm_day=0.0,
        shortwave_down=[120.0, 140.0, 160.0],
        wind_speed=3.0,
        time_values=[DateTime(2001, 1, 15, 12), DateTime(2001, 2, 15, 12), DateTime(2001, 3, 15, 12)],
    )
    layout = Chion.GridLayout(
        [0.0, 10_000.0],
        [0.0, 10_000.0],
        [1, 1, 2, 2],
        [1, 2, 1, 2],
        ones(2, 2),
    )
    return domain, forcing, layout
end

function _write_sample_gris_forcing_file(path::AbstractString)
    nx = 2
    ny = 2
    nlayer = 2
    ntime = 2
    h5open(path, "w") do file
        file["x"] = [0.0, 10_000.0]
        file["y"] = [0.0, 10_000.0]
        file["MSK"] = fill(1.0, nx, ny)
        file["OUTLAY_bnds"] = [0.0 0.5; 0.5 1.0]

        tt = fill(-15.0, nx, ny, ntime)
        tt[:, :, 2] .= -14.0
        file["TT"] = tt
        file["SF"] = fill(0.0, nx, ny, ntime)
        file["RF"] = fill(0.0, nx, ny, ntime)
        file["SWD"] = fill(100.0, nx, ny, ntime)
        file["LWD"] = fill(250.0, nx, ny, ntime)
        file["SHF"] = fill(0.0, nx, ny, ntime)
        file["LHF"] = fill(0.0, nx, ny, ntime)

        zn3 = fill(1.0, nx, ny, ntime)
        file["ZN3"] = zn3

        ro1 = fill(0.0, nx, ny, nlayer, ntime)
        ro1[:, :, 1, :] .= 300.0
        ro1[:, :, 2, :] .= 340.0
        file["RO1"] = ro1

        ti1 = fill(0.0, nx, ny, nlayer, ntime)
        ti1[:, :, 1, :] .= -12.0
        ti1[:, :, 2, :] .= -8.0
        file["TI1"] = ti1

        wa1 = fill(0.0, nx, ny, nlayer, ntime)
        file["WA1"] = wa1

        yyyy = fill(2001, ntime)
        mm = fill(1, ntime)
        dd = fill(1, ntime)
        dd[2] = 2
        hh = fill(12, ntime)
        file["YYYY"] = yyyy
        file["MM"] = mm
        file["DD"] = dd
        file["HH"] = hh
    end
    return path
end

@testset "Run API" begin
    @testset "physics helper" begin
        actual = Chion.physics(
            albedo=:constant,
            densification=:htessel,
            fresh_snow_density=:parameterized,
        )
        expected = Chion.SnowpackPhysicalConstants(
            albedo_scheme=:constant,
            low_density_densification=:htessel,
            fresh_snow_density_scheme=:parameterized,
        )
        @test actual.albedo_scheme == expected.albedo_scheme
        @test actual.low_density_densification == expected.low_density_densification
        @test actual.fresh_snow_density_scheme == expected.fresh_snow_density_scheme
    end

    @testset "ForcingData conversions and validation" begin
        forcing = Chion.ForcingData(
            dt_days=[1.0, 2.0],
            ncol=2,
            air_temperature_c=[-10.0, -5.0],
            snowfall_mm_day=1.5,
            rainfall_mm_day=[0.0, 2.0],
            shortwave_down=[100.0 120.0; 140.0 160.0],
            wind_speed=3.5,
            q_lw_down=[250.0, 255.0],
            has_q_lw_down=[true, false],
        )

        @test size(forcing.air_temperature) == (2, 2)
        @test forcing.air_temperature[1, 1] ≈ 263.15 atol = 1e-8
        @test forcing.air_temperature[2, 2] ≈ 268.15 atol = 1e-8
        @test forcing.snowfall_rate[2, 1] ≈ 1.5 / 86_400.0 atol = 1e-12
        @test forcing.rainfall_rate[1, 2] ≈ 2.0 / 86_400.0 atol = 1e-12
        @test forcing.shortwave_down == [100.0 120.0; 140.0 160.0]
        @test forcing.wind_speed == fill(3.5, 2, 2)
        @test forcing.has_q_lw_down == [true false; true false]
        @test length(forcing.time_values) == 2

        err = _captured_exception() do
            Chion.ForcingData(
                dt_days=[1.0, 1.0],
                ncol=2,
                air_temperature_c=[-10.0],
                snowfall_mm_day=0.0,
                rainfall_mm_day=0.0,
                shortwave_down=0.0,
            )
        end
        @test err isa Exception
        @test occursin("air_temperature_c", sprint(showerror, err))
    end

    @testset "run! saves an exact symbol" begin
        mktempdir() do dir
            domain, forcing, layout = _sample_domain_forcing_layout()
            result = Chion.run!(
                domain,
                forcing;
                layout=layout,
                save=:final_thickness,
                output_dir=dir,
                netcdf_path=joinpath(dir, "save_exact.nc"),
                backend=:cpu,
                write_outputs=false,
                cycles=1,
                history_stride=1,
            )

            @test result.status == :cycles
            @test result.run.backend == :threads
            @test isfile(result.netcdf_path)

            ds = NCDataset(result.netcdf_path)
            @test haskey(ds, "final_thickness")
            @test !haskey(ds, "history_mean_thickness")
            @test size(ds["final_thickness"]) == (2, 2)
            close(ds)
        end
    end

    @testset "run! saves a group" begin
        mktempdir() do dir
            domain, forcing, layout = _sample_domain_forcing_layout()
            result = Chion.run!(
                domain,
                forcing;
                layout=layout,
                save=:history,
                output_dir=dir,
                netcdf_path=joinpath(dir, "save_group.nc"),
                backend=:threads,
                write_outputs=false,
                cycles=1,
                history_stride=1,
            )

            @test result.status == :cycles
            ds = NCDataset(result.netcdf_path)
            @test haskey(ds, "history_mean_thickness")
            @test haskey(ds, "history_mean_base_mass")
            @test !haskey(ds, "final_thickness")
            close(ds)
        end
    end

    @testset "run! saves mixed fields" begin
        mktempdir() do dir
            domain, forcing, layout = _sample_domain_forcing_layout()
            result = Chion.run!(
                domain,
                forcing;
                layout=layout,
                save=[:final_thickness, :history],
                output_dir=dir,
                netcdf_path=joinpath(dir, "save_mixed.nc"),
                backend=:threads,
                write_outputs=false,
                cycles=1,
                history_stride=1,
            )

            @test result.status == :cycles
            ds = NCDataset(result.netcdf_path)
            @test haskey(ds, "final_thickness")
            @test haskey(ds, "history_mean_thickness")
            close(ds)
        end
    end

    @testset "run! can skip NetCDF entirely" begin
        domain, forcing, _ = _sample_domain_forcing_layout()
        result = Chion.run!(
            domain,
            forcing;
            save=:none,
            backend=:threads,
            write_outputs=false,
            cycles=1,
        )

        @test result.status == :cycles
        @test result.netcdf_path == ""
        @test result.summary_path == ""
        @test result.history_csv_path == ""
    end

    @testset "load_gris_forcing_file_problem and run!(problem)" begin
        mktempdir() do dir
            forcing_path = _write_sample_gris_forcing_file(joinpath(dir, "gris_sample.nc"))
            problem = Chion.load_gris_forcing_file_problem(
                forcing_path;
                mask_threshold=0.0,
                ntot=4,
                physics=Chion.physics(),
            )

            @test problem isa Chion.LoadedProblem
            @test length(problem.notes) >= 1
            @test problem.metadata.path == forcing_path

            result_problem = Chion.run!(
                problem;
                save=:none,
                backend=:threads,
                write_outputs=false,
                cycles=1,
            )

            @test result_problem.status == :cycles
            @test result_problem.netcdf_path == ""

            result_direct = Chion.run!(
                deepcopy(problem.domain),
                problem.forcing;
                layout=problem.layout,
                save=:none,
                backend=:threads,
                write_outputs=false,
                cycles=1,
            )

            @test result_direct.status == :cycles
        end
    end
end
