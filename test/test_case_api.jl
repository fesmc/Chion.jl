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

        time = defVar(ds, "time", Float64, ("time",))
        time.attrib["units"] = "days since 2001-01-01 12:00:00"
        time[:] = [0.0, 1.0]

        tt = fill(-15.0, ntime, ny, nx); tt[2, :, :] .= -14.0
        defVar(ds, "TT", Float64, ("time", "y", "x"))[:, :, :] = tt
        defVar(ds, "SF", Float64, ("time", "y", "x"))[:, :, :] = fill(1.0, ntime, ny, nx)
        defVar(ds, "RF", Float64, ("time", "y", "x"))[:, :, :] = fill(0.0, ntime, ny, nx)
        defVar(ds, "SWD", Float64, ("time", "y", "x"))[:, :, :] = fill(100.0, ntime, ny, nx)
        al2 = fill(0.55, ntime, ny, nx); al2[2, :, :] .= 0.60
        defVar(ds, "AL2", Float64, ("time", "y", "x"))[:, :, :] = al2
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
        m1 = BESSIModel(grid; albedo=DynamicAlbedo())
        m2 = BESSIModel(grid; albedo=ConstantAlbedo())
        m_prescribed = BESSIModel(grid; albedo=PrescribedAlbedo())
        m3 = BESSIModel(grid; densification=HTESSELDensification())
        m4 = BESSIModel(grid; fresh_snow_density=ParameterizedFreshSnowDensity())
        @test m1.c.albedo_scheme == Chion.ALBEDO_DYNAMIC
        @test m2.c.albedo_scheme == Chion.ALBEDO_CONSTANT
        @test m_prescribed.c.albedo_scheme == Chion.ALBEDO_PRESCRIBED
        @test m3.c.low_density_densification == Chion.LOW_DENSIFICATION_HTESSEL
        @test m4.c.fresh_snow_density_scheme == Chion.FRESH_SNOW_DENSITY_PARAMETERIZED
    end

    @testset "BESSI depth cap uses reference snow depth, not total mass" begin
        c = Chion.SnowpackPhysicalConstants()
        cap_depth = Chion.BESSI_REFERENCE_LAYER_COUNT * 1.5 * Chion.DEFAULT_MASS_SPLIT /
                    Chion.BESSI_REFERENCE_DEPTH_DENSITY

        N = [5]
        mass = fill(2000.0, 5, 1)
        mass_w = zeros(5, 1)
        density = fill(700.0, 5, 1)
        temperature = fill(c.T0, 5, 1)
        mass_base = [0.0]
        smb_ice = [0.0]
        runoff = [0.0]
        Tsrf = [c.T0]
        albedo_dynamic = [c.alpha_dry]

        Chion._enforce_snow_depth_cap!(
            N, mass, mass_w, density, temperature,
            mass_base, smb_ice, runoff, Tsrf, albedo_dynamic,
            1, 5, Chion.DEFAULT_MASS_SPLIT, Chion.DEFAULT_SECONDS_PER_DAY, c,
        )
        @test sum(mass[:, 1]) ≈ 10_000.0
        @test sum(mass[:, 1] ./ density[:, 1]) < cap_depth
        @test mass_base[1] == 0.0

        fill!(mass, 1500.0)
        fill!(density, 350.0)
        fill!(mass_base, 0.0)
        fill!(smb_ice, 0.0)

        Chion._enforce_snow_depth_cap!(
            N, mass, mass_w, density, temperature,
            mass_base, smb_ice, runoff, Tsrf, albedo_dynamic,
            1, 5, Chion.DEFAULT_MASS_SPLIT, Chion.DEFAULT_SECONDS_PER_DAY, c,
        )
        @test sum(mass[:, 1] ./ density[:, 1]) ≈ cap_depth atol=1e-12
        @test mass_base[1] ≈ 750.0 atol=1e-10
        @test smb_ice[1] ≈ 750.0 atol=1e-10
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
            prescribed_albedo = [0.45, 0.50],
        )
        @test size(forcing.air_temperature) == (2, 2)
        @test forcing.air_temperature[1, 1] ≈ 263.15 atol=1e-8
        @test forcing.air_temperature[2, 2] ≈ 268.15 atol=1e-8
        @test forcing.snowfall_rate[2, 1]   ≈ 1.5 / 86_400.0 atol=1e-12
        @test forcing.rainfall_rate[1, 2]   ≈ 2.0 / 86_400.0 atol=1e-12
        @test forcing.shortwave_down == [100.0 120.0; 140.0 160.0]
        @test forcing.wind_speed == fill(3.5, 2, 2)
        @test forcing.has_q_lw_down == [true false; true false]
        @test forcing.prescribed_albedo == [0.45 0.50; 0.45 0.50]
        @test all(forcing.has_prescribed_albedo)
        @test length(forcing.time_values) == 2

        latitude_forcing = SnowpackForcing(
            dt_days=[1.0, 1.0],
            latitude_deg=[60.0, 65.0],
            air_temperature_c=[-10.0, -9.0],
            snowfall_mm_day=0.0,
            rainfall_mm_day=0.0,
            shortwave_down=100.0,
            time_values=[DateTime(2001, 3, 20, 12), DateTime(2001, 3, 21, 12)],
        )
        @test size(latitude_forcing.latitude_deg) == (2, 2)
        @test latitude_forcing.latitude_deg == [60.0 60.0; 65.0 65.0]
        @test latitude_forcing.day_of_year[1] != latitude_forcing.solar_longitude_deg[1]

        latitude_matrix_forcing = SnowpackForcing(
            dt_days=[1.0, 1.0],
            ncol=2,
            latitude_deg=[60.0 61.0; 65.0 66.0],
            air_temperature_c=[-10.0, -9.0],
            snowfall_mm_day=0.0,
            rainfall_mm_day=0.0,
            shortwave_down=100.0,
        )
        @test latitude_matrix_forcing.latitude_deg == [60.0 61.0; 65.0 66.0]

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

    @testset "adaptive diurnal shortwave substeps conserve forcing energy" begin
        terms = Chion._diurnal_shortwave_integral_terms(65.0, 90.0)
        sw_mean = 200.0
        scale = sw_mean * 2π / terms.daylight_integral
        @test scale * terms.daylight_integral / 2π ≈ sw_mean atol=1e-12
        temp_mean = 273.15 - 2.0
        temp_amplitude = 6.0
        temp_morning = Chion._diurnal_temperature_interval_average(temp_mean, temp_amplitude, -π, -π / 2)
        temp_midday = Chion._diurnal_temperature_interval_average(temp_mean, temp_amplitude, -π / 2, π / 2)
        temp_evening = Chion._diurnal_temperature_interval_average(temp_mean, temp_amplitude, π / 2, π)
        @test (temp_morning * (π / 2) + temp_midday * π + temp_evening * (π / 2)) / 2π ≈ temp_mean atol=1e-12
        @test temp_midday > temp_mean
        @test temp_morning < temp_mean
        morning_night = Chion._diurnal_shortwave_interval_average(sw_mean, 65.0, 90.0, -π, -terms.sunset_hour_angle)
        daylight = Chion._diurnal_shortwave_interval_average(sw_mean, 65.0, 90.0, -terms.sunset_hour_angle, terms.sunset_hour_angle)
        evening_night = Chion._diurnal_shortwave_interval_average(sw_mean, 65.0, 90.0, terms.sunset_hour_angle, π)
        weighted_substeps = (
            morning_night * (π - terms.sunset_hour_angle) +
            daylight * (2 * terms.sunset_hour_angle) +
            evening_night * (π - terms.sunset_hour_angle)
        ) / 2π
        @test weighted_substeps ≈ sw_mean atol=1e-12
        warm_air = 273.15 - 2.0
        cold_air = 273.15 - 9.0
        min_substep_air = 273.15 - 8.0
        @test Chion._diurnal_shortwave_substep_count(1.0, sw_mean, warm_air, min_substep_air, 65.0, 90.0, 0.0, 24) == 24
        @test Chion._diurnal_shortwave_substep_count(1.0, sw_mean, warm_air, min_substep_air, 65.0, 90.0, 0.0, 1) == 1
        @test Chion._diurnal_shortwave_substep_count(1.0, sw_mean, cold_air, min_substep_air, 65.0, 90.0, 0.0, 24) == 1
        @test Chion._diurnal_shortwave_substep_count(1.0, sw_mean, warm_air, min_substep_air, 65.0, 90.0, 10_000.0, 24) == 1

        grid = SnowpackGrid(1)
        missing_latitude_forcing = SnowpackForcing(
            dt_days=[1.0],
            air_temperature_c=[-5.0],
            snowfall_mm_day=0.0,
            rainfall_mm_day=0.0,
            shortwave_down=100.0,
        )
        err = _captured_exception() do
            run!(Simulation(
                BESSIModel(grid; diurnal_shortwave=true);
                forcing=missing_latitude_forcing,
                save=:none,
                years=1,
                write_outputs=false,
            ); io=devnull)
        end
        @test err isa Exception
        @test occursin("latitude_deg", sprint(showerror, err))

        err = _captured_exception() do
            run!(Simulation(
                BESSIModel(grid; diurnal_shortwave_substeps=true);
                forcing=missing_latitude_forcing,
                save=:none,
                years=1,
                write_outputs=false,
            ); io=devnull)
        end
        @test err isa Exception
        @test occursin("latitude_deg", sprint(showerror, err))

        adaptive_model = BESSIModel(grid;
            Ntot=3,
            diurnal_shortwave_substeps=true,
            diurnal_shortwave_max_substeps=4,
            diurnal_temperature_cycle=true,
            diurnal_temperature_amplitude_c=6.0,
        )
        adaptive_forcing = SnowpackForcing(
            dt_days=[1.0],
            air_temperature_c=[-2.0],
            snowfall_mm_day=5.0,
            rainfall_mm_day=0.0,
            shortwave_down=200.0,
            latitude_deg=65.0,
            time_values=[DateTime(2001, 6, 21, 12)],
        )
        adaptive_result = run!(Simulation(
            adaptive_model;
            forcing=adaptive_forcing,
            save=:none,
            years=1,
            write_outputs=false,
        ); io=devnull)
        @test adaptive_result.status == :complete

        alias_model = BESSIModel(grid; Ntot=3, diurnal_shortwave=true)
        @test alias_model.diurnal_shortwave_substeps
        @test BESSIModel(grid; Ntot=3, diurnal_shortwave_substeps=true, diurnal_shortwave_max_substeps=1).diurnal_shortwave_max_substeps == 1
        @test BESSIModel(grid; Ntot=3, diurnal_shortwave_substeps=true, diurnal_shortwave_max_substeps=24).diurnal_shortwave_max_substeps == 24
        @test adaptive_model.diurnal_temperature_cycle
        @test adaptive_model.diurnal_temperature_amplitude == 6.0
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

    @testset "monthly SMB is climatic SMB" begin
        mktempdir() do dir
            grid = SnowpackGrid(1;
                x=[0.0],
                y=[0.0],
                js=[1],
                is=[1],
                mask=ones(1, 1),
            )
            forcing = SnowpackForcing(
                dt_days=[1.0],
                air_temperature_c=[-20.0],
                snowfall_mm_day=[10.0],
                rainfall_mm_day=[0.0],
                shortwave_down=[0.0],
                time_values=[DateTime(2001, 1, 15, 12)],
            )
            sim = Simulation(BESSIModel(grid; Ntot=4); forcing=forcing,
                save=:monthly,
                output_dir=dir,
                netcdf_path=joinpath(dir, "monthly_smb.nc"),
                backend=:threads,
                write_outputs=false,
                years=1,
            )
            result = run!(sim; io=devnull)
            @test result.status == :complete
            ds = NCDataset(result.netcdf_path)
            @test haskey(ds, "monthly_smb")
            @test haskey(ds, "monthly_net_ice_sheet_forcing")
            @test !haskey(ds, "monthly_mean_ice_sheet_smb")
            @test !haskey(ds, "monthly_ice_sheet_smb")
            @test ds["monthly_smb"][1, 1, 1] ≈ 10.0 atol=1e-6
            @test ds["monthly_net_ice_sheet_forcing"][1, 1, 1] ≈ 0.0 atol=1e-6
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
        @test sim.now isa CurrentState
        @test sim.ref isa ReferenceState
        @test !hasproperty(sim.now, :domain)
        @test all(sim.ref.mass .== 0.0)
        run!(sim; io=devnull)
        @test all(sim.ref.mass .== 0.0)
        @test sum(sim.now.mass) > 0.0
        @test get_state(sim, 1)["total_mass"] > 0.0
        @test sim.now.thickness[1] > 0.0
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
        @test sim_run.now.mass ≈ sim_manual.now.mass
        @test sim_run.now.smb_ice ≈ sim_manual.now.smb_ice
    end

    @testset "checkpoint restart matches uninterrupted run" begin
        mktempdir() do dir
            checkpoint_path = joinpath(dir, "restart.chk")
            model_full, forcing, _ = _sample_model_forcing_grid()
            full = Simulation(model_full; forcing=forcing, save=:none, years=2, write_outputs=false)
            full_result = run!(full; io=devnull)

            model_restart = BESSIModel(model_full.grid; Ntot=4)
            restarted = Simulation(model_restart; forcing=forcing, save=:none, years=2, write_outputs=false)
            integrator = init_integrator(restarted; io=devnull)
            step!(integrator, length(forcing.time_values))
            checkpoint!(integrator, checkpoint_path)

            resumed = restart_integrator(checkpoint_path; io=devnull)
            @test resumed.completed_years == 1
            run!(resumed)
            restart_result = finalize!(resumed)

            @test restart_result.status == full_result.status == :complete
            @test restart_result.years_completed == full_result.years_completed == 2
            @test resumed.sim.now.mass ≈ full.now.mass
            @test resumed.sim.now.smb_ice ≈ full.now.smb_ice
            @test restart_result.history[end].mean_thickness ≈ full_result.history[end].mean_thickness
        end
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
        @test external.now.mass ≈ scheduled.now.mass
        @test external.now.smb_ice ≈ scheduled.now.smb_ice
    end

    @testset "active mask limits initialized BESSI stepping to selected columns" begin
        grid = SnowpackGrid(4;
            x=[0.0, 1.0],
            y=[0.0, 1.0],
            js=[1, 1, 2, 2],
            is=[1, 2, 1, 2],
            mask=ones(2, 2),
        )
        forcing = SnowpackForcing(
            dt_days=[1.0],
            ncol=4,
            air_temperature_c=fill(-5.0, 4, 1),
            snowfall_mm_day=fill(10.0, 4, 1),
            rainfall_mm_day=0.0,
            shortwave_down=0.0,
        )
        sim = Simulation(BESSIModel(grid; Ntot=4); forcing=forcing, save=:none, years=1, write_outputs=false)
        integrator = init_integrator(sim; io=devnull)
        set_active_mask!(integrator, Bool[true, false, true, false])
        step!(integrator)

        @test sim.now.mass[1, 1] > 0.0
        @test sim.now.mass[1, 2] == 0.0
        @test sim.now.mass[1, 3] > 0.0
        @test sim.now.mass[1, 4] == 0.0
    end

    @testset "prescribed zero latent heat keeps precipitation heat terms" begin
        function run_precip_heat_case(; q_lh=nothing)
            grid = SnowpackGrid(1)
            forcing = SnowpackForcing(
                dt_days=[1.0],
                ncol=1,
                air_temperature_c=[-5.0],
                snowfall_mm_day=[10.0],
                rainfall_mm_day=0.0,
                shortwave_down=0.0,
                q_lh=q_lh,
            )
            sim = Simulation(BESSIModel(grid; Ntot=4); forcing=forcing, save=:none, years=1, write_outputs=false)
            run!(sim; io=devnull)
            return sim.now.temperature[1, 1]
        end

        @test run_precip_heat_case() ≈ run_precip_heat_case(q_lh=[0.0]) atol=0.0
    end

    @testset "relative humidity latent heat fallback is used only without q_lh" begin
        function run_latent_vapor_case(; relative_humidity=nothing, q_lh=nothing)
            grid = SnowpackGrid(1)
            forcing = SnowpackForcing(
                dt_days=[1.0],
                ncol=1,
                air_temperature_c=[-5.0],
                relative_humidity=relative_humidity,
                snowfall_mm_day=[10.0],
                rainfall_mm_day=0.0,
                shortwave_down=0.0,
                q_lh=q_lh,
            )
            sim = Simulation(BESSIModel(grid; Ntot=4); forcing=forcing, save=:none, years=1, write_outputs=false)
            run!(sim; io=devnull)
            return sim.now.temperature[1, 1]
        end

        without_relative_humidity = run_latent_vapor_case()
        with_relative_humidity = run_latent_vapor_case(relative_humidity=[100.0])
        prescribed_zero = run_latent_vapor_case(relative_humidity=[100.0], q_lh=[0.0])

        @test with_relative_humidity > without_relative_humidity
        @test prescribed_zero ≈ without_relative_humidity atol=0.0
    end

    @testset "snow-covered cells include turbulent vapor mass flux" begin
        function run_snow_vapor_case(; relative_humidity=nothing, q_lh=nothing)
            grid = SnowpackGrid(1)
            forcing = SnowpackForcing(
                dt_days=[1.0],
                ncol=1,
                air_temperature_c=[-20.0],
                relative_humidity=relative_humidity,
                snowfall_mm_day=[100.0],
                rainfall_mm_day=0.0,
                shortwave_down=0.0,
                q_lh=q_lh,
            )
            sim = Simulation(BESSIModel(grid; Ntot=4); forcing=forcing, save=:none, years=1, write_outputs=false)
            run!(sim; io=devnull)
            return sim.now.mass[1, 1]
        end

        without_relative_humidity = run_snow_vapor_case()
        with_dry_air = run_snow_vapor_case(relative_humidity=[0.0])
        prescribed_zero = run_snow_vapor_case(relative_humidity=[0.0], q_lh=[0.0])

        @test with_dry_air < without_relative_humidity
        @test prescribed_zero ≈ without_relative_humidity atol=0.0
    end

    @testset "bare ice includes vapor mass flux bookkeeping" begin
        function run_bare_ice_case(; relative_humidity=nothing)
            grid = SnowpackGrid(1)
            forcing = SnowpackForcing(
                dt_days=[1.0],
                ncol=1,
                air_temperature_c=[0.0],
                relative_humidity=relative_humidity,
                snowfall_mm_day=[0.0],
                rainfall_mm_day=0.0,
                shortwave_down=0.0,
            )
            sim = Simulation(BESSIModel(grid; Ntot=4); forcing=forcing, save=:none, years=1, write_outputs=false)
            run!(sim; io=devnull)
            return sim.now.smb_ice[1], sim.now.melt[1], sim.now.runoff[1]
        end

        smb_none, melt_none, runoff_none = run_bare_ice_case()
        smb_rh, melt_rh, runoff_rh = run_bare_ice_case(relative_humidity=[100.0])

        @test smb_none == 0.0
        @test melt_none == 0.0
        @test runoff_none == 0.0
        @test smb_rh < 0.0
        @test melt_rh == 0.0
        @test runoff_rh == 0.0
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
        result = run!(Simulation(model; forcing=forcing, save=:none, years=1, write_outputs=false))
        @test result.status == :complete
        @test result.netcdf_path == ""

        err = _captured_exception() do
            run!(Simulation(model; forcing=forcing, save=:final_thickness, years=1, write_outputs=false))
        end
        @test err isa Exception
        @test occursin("spatial coordinates", sprint(showerror, err))
    end

    @testset "prescribed albedo forcing drives surface albedo" begin
        grid = SnowpackGrid(1)
        forcing = SnowpackForcing(
            dt_days=[1.0],
            air_temperature_c=[-10.0],
            snowfall_mm_day=1.0,
            rainfall_mm_day=0.0,
            shortwave_down=100.0,
            prescribed_albedo=[0.42],
        )
        sim = Simulation(
            BESSIModel(grid; Ntot=4, albedo=PrescribedAlbedo());
            forcing=forcing,
            save=:none,
            years=1,
            write_outputs=false,
        )
        result = run!(sim; io=devnull)
        @test result.status == :complete
        @test sim.now.albedo_dynamic[1] ≈ 0.42 atol=1e-12

        missing_albedo = SnowpackForcing(
            dt_days=[1.0],
            air_temperature_c=[-10.0],
            snowfall_mm_day=1.0,
            rainfall_mm_day=0.0,
            shortwave_down=100.0,
        )
        err = _captured_exception() do
            run!(Simulation(
                BESSIModel(grid; Ntot=4, albedo=PrescribedAlbedo());
                forcing=missing_albedo,
                save=:none,
                years=1,
                write_outputs=false,
            ); io=devnull)
        end
        @test err isa Exception
        @test occursin("prescribed_albedo", sprint(showerror, err))
    end

    @testset "forcing-file loading and run!" begin
        mktempdir() do dir
            forcing_path = _write_sample_forcing_file(joinpath(dir, "forcing_sample.nc"))
            loaded = load_forcing_file(forcing_path; prescribed_albedo_name="AL2")
            @test size(loaded.forcing.air_temperature) == (4, 2)
            @test size(loaded.forcing.snowfall_rate) == (4, 2)
            @test loaded.forcing.prescribed_albedo[1, 1] ≈ 0.55 atol=1e-12
            @test loaded.forcing.prescribed_albedo[1, 2] ≈ 0.60 atol=1e-12
            @test all(loaded.forcing.has_prescribed_albedo)
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
        grid = SnowpackGrid(1)
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

        monthly_f = SnowpackForcing(
            dt_days=[30.0],
            air_temperature_c=[0.0],
            snowfall_mm_day=0.0,
            rainfall_mm_day=0.0,
            shortwave_down=150.0,
        )
        monthly_pdd = PDDModel(grid;
            ddf_snow=3.0,
            ddf_ice=8.0,
            refreezing_fraction=0.0,
            monthly_method=StochasticMonthlyPDD(temperature_sigma=5.0),
        )
        monthly_sim = Simulation(monthly_pdd; forcing=monthly_f, save=:none, years=1, write_outputs=false)
        @test run!(monthly_sim).status == :complete
        expected_monthly_pdd = 30.0 * 5.0 / sqrt(2 * pi)
        @test monthly_sim.now.pdd_sum[1] ≈ expected_monthly_pdd rtol=1e-6
        @test monthly_sim.now.smb_ice[1] ≈ -8.0 * expected_monthly_pdd rtol=1e-6
        @test monthly_sim.now.runoff[1] ≈ 8.0 * expected_monthly_pdd rtol=1e-6

        named = Simulation(:pdd, grid;
            model_kwargs=(ddf_snow=3.0, ddf_ice=8.0, refreezing_fraction=0.0),
            forcing=f,
            save=:none,
            years=1,
            write_outputs=false,
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
        grid = SnowpackGrid(1)
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
