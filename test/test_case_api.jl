using Test
using Dates
using HDF5
using Chion

function _captured_exception(f::Function)
    try
        f()
        return nothing
    catch err
        return err
    end
end

function _write_mar_fixture(path::AbstractString)
    ntime = 3
    nlayer = 2
    ny = 2
    nx = 2

    x = [0.0, 10_000.0]
    y = [0.0, 10_000.0]
    mask = [
        100.0 100.0
        20.0 100.0
    ]
    outlay_bounds = [
        0.0 0.4
        0.4 0.8
    ]

    tt = Array{Float64}(undef, ntime, ny, nx)
    sf = fill(0.5, ntime, ny, nx)
    rf = fill(0.0, ntime, ny, nx)
    swd = Array{Float64}(undef, ntime, ny, nx)
    lwd = fill(250.0, ntime, ny, nx)
    shf = fill(0.0, ntime, ny, nx)
    lhf = fill(0.0, ntime, ny, nx)
    zn3 = fill(0.8, ntime, ny, nx)
    ro1 = fill(0.0, ntime, nlayer, ny, nx)
    ti1 = fill(0.0, ntime, nlayer, ny, nx)
    wa1 = fill(0.0, ntime, nlayer, ny, nx)

    for t in 1:ntime, j in 1:ny, i in 1:nx
        tt[t, j, i] = -16.0 + 1.5 * (t - 1) + 0.4 * (i - 1) - 0.3 * (j - 1)
        swd[t, j, i] = 100.0 + 20.0 * (t - 1) + 5.0 * (i - 1)
    end
    rf[3, :, :] .= 0.2
    ro1[:, 1, :, :] .= 320.0
    ro1[:, 2, :, :] .= 450.0
    ti1[:, 1, :, :] .= -8.0
    ti1[:, 2, :, :] .= -12.0
    wa1[:, 2, :, :] .= 0.02

    yyyy = [2025, 2025, 2025]
    mm = [1, 1, 1]
    dd = [1, 2, 3]
    hh = [12, 12, 12]

    h5open(path, "w") do file
        write(file, "x", x)
        write(file, "y", y)
        write(file, "MSK", permutedims(mask, (2, 1)))
        write(file, "OUTLAY_bnds", permutedims(outlay_bounds, (2, 1)))
        write(file, "TT", permutedims(tt, (3, 2, 1)))
        write(file, "SF", permutedims(sf, (3, 2, 1)))
        write(file, "RF", permutedims(rf, (3, 2, 1)))
        write(file, "SWD", permutedims(swd, (3, 2, 1)))
        write(file, "LWD", permutedims(lwd, (3, 2, 1)))
        write(file, "SHF", permutedims(shf, (3, 2, 1)))
        write(file, "LHF", permutedims(lhf, (3, 2, 1)))
        write(file, "ZN3", permutedims(zn3, (3, 2, 1)))
        write(file, "RO1", permutedims(ro1, (4, 3, 2, 1)))
        write(file, "TI1", permutedims(ti1, (4, 3, 2, 1)))
        write(file, "WA1", permutedims(wa1, (4, 3, 2, 1)))
        write(file, "YYYY", yyyy)
        write(file, "MM", mm)
        write(file, "DD", dd)
        write(file, "HH", hh)
    end

    return path
end

@testset "Case API" begin
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

    @testset "prescribed_case builds case definitions" begin
        definition = Chion.prescribed_case(
            physics=Chion.physics(),
            ntot=4,
            nx=2,
            ny=2,
            dt_days=[1.0, 1.0, 1.0],
            air_temperature_c=[-15.0, -14.0, -13.0],
            snowfall_mm_day=1.0,
            rainfall_mm_day=0.0,
            shortwave_down=[120.0, 140.0, 160.0],
            initial_surface_mass=[120.0, 0.0, 180.0, 60.0],
            initial_density=330.0,
            initial_temperature_c=-11.0,
            input_label="prescribed_demo",
        )

        @test definition isa Chion.CaseDefinition
        @test definition.metadata.format == :prescribed
        @test definition.metadata.ncol == 4
        @test definition.input_label == "prescribed_demo"
        @test size(definition.forcing.air_temperature) == (4, 3)
        @test length(definition.layout.js) == 4
        @test Chion.get_state(definition.domain, 1)["n_active"] == 1
        @test Chion.get_state(definition.domain, 2)["n_active"] == 0
    end

    @testset "synthetic_case runs through build_case/run_case" begin
        definition = Chion.synthetic_case(variant=:multi_column, ntot=4, ntime=4, nx=2, ny=2)
        case = Chion.build_case(
            definition;
            run=Chion.RunConfig(
                name="synthetic_smoke",
                backend=:cpu,
                write_outputs=false,
                write_netcdf=false,
                cycles=1,
            ),
        )
        result = Chion.run_case(case; io=devnull)

        @test result.status == :cycles
        @test length(result.history) == 1
        @test result.summary_path == ""
        @test result.history_csv_path == ""
        @test result.run.backend == :threads
        @test case.definition === definition
    end

    @testset "MAR fixture loads and runs" begin
        mktempdir() do dir
            fixture_path = _write_mar_fixture(joinpath(dir, "mar_fixture.h5"))
            definition = Chion.mar_case(fixture_path; ntot=4, physics=Chion.physics(), mask_threshold=50.0)

            @test definition.metadata.format == :mar
            @test definition.metadata.ncol == 3
            @test definition.metadata.ntime == 3
            @test definition.input_label == fixture_path
            @test size(definition.forcing.air_temperature) == (3, 3)
            @test length(definition.layout.js) == 3
            @test all(definition.forcing.wind_speed .== 5.0)
            @test occursin("default 5.0 m s^-1", only(definition.notes))

            case = Chion.build_case(
                definition;
                run=Chion.RunConfig(
                    name="mar_fixture_smoke",
                    backend=:cpu,
                    write_outputs=false,
                    write_netcdf=false,
                    cycles=1,
                ),
            )
            result = Chion.run_case(case; io=devnull)

            @test result.status == :cycles
            @test length(result.history) == 1
            @test result.run.backend == :threads
        end
    end

    @testset "Optional real MAR smoke" begin
        real_mar_path = get(
            ENV,
            "CHION_REAL_MAR_PATH",
            "/p/projects/ou/labs/ai/Nils/MARv3.14.3-10km-daily-ERA5-2025.nc",
        )
        enabled = get(ENV, "CHION_RUN_REAL_MAR_SMOKE", "0") == "1"
        if enabled && isfile(real_mar_path)
            definition = Chion.mar_case(real_mar_path; ntot=10, physics=Chion.physics(), mask_threshold=50.0)
            case = Chion.build_case(
                definition;
                run=Chion.RunConfig(
                    name="real_mar_smoke",
                    backend=:cpu,
                    write_outputs=false,
                    write_netcdf=false,
                    cycles=1,
                ),
            )
            result = Chion.run_case(case; io=devnull)
            @test result.status == :cycles
        else
            @info "Skipping optional real MAR smoke test" enabled=enabled path=real_mar_path
            @test true
        end
    end
end
