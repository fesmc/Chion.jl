#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Chion
using NCDatasets
include(joinpath(@__DIR__, "..", "shared", "script_helpers.jl"))
using .ChionExampleScriptHelpers

const DEFAULT_GRIS_FORCING_FILE_OUTPUT_DIR = joinpath(@__DIR__, "..", "plots", "gris_forcing_file_simulation")

function default_gris_forcing_file_path()
    for candidate in (
        get(ENV, "FORCING_PATH", ""),
        "/Users/niboch001/Downloads/MARv3.14.3-10km-daily-ERA5-2025.nc",
        "/Users/niboch001/Downloads/MARv3.14.3-10km-daily-ERA5-2026.nc",
    )
        isfile(candidate) && return candidate
    end
    return ""
end

function apply_gris_mask(loaded, forcing_file::AbstractString; threshold::Float64=50.0)
    mask = NCDataset(forcing_file) do ds
        Float64.(ds["MSK"][:])
    end

    ny = length(loaded.grid.y)
    nx = length(loaded.grid.x)
    if length(mask) == nx * ny && ndims(mask) == 1
        mask = reshape(mask, nx, ny)
    end
    if size(mask) == (nx, ny)
        mask = permutedims(mask, (2, 1))
    end
    size(mask) == (ny, nx) || error("Unexpected MSK shape $(size(mask)); expected ($ny, $nx).")

    rows = Int[]
    js = Int[]
    is = Int[]
    for j in 1:ny, i in 1:nx
        row = i + (j - 1) * nx
        if isfinite(mask[j, i]) &&
           mask[j, i] >= threshold &&
           isfinite(loaded.forcing.air_temperature[row, 1])
            push!(rows, row)
            push!(js, j)
            push!(is, i)
        end
    end

    f = loaded.forcing
    forcing = SnowpackForcing(
        time_values=f.time_values,
        dt_days=f.dt_days,
        air_temperature=f.air_temperature[rows, :],
        snowfall_rate=f.snowfall_rate[rows, :],
        rainfall_rate=f.rainfall_rate[rows, :],
        shortwave_down=f.shortwave_down[rows, :],
        wind_speed=f.wind_speed[rows, :],
        q_lw_down=f.q_lw_down[rows, :],
        has_q_lw_down=f.has_q_lw_down[rows, :],
        q_sh=f.q_sh[rows, :],
        has_q_sh=f.has_q_sh[rows, :],
        q_lh=f.q_lh[rows, :],
        has_q_lh=f.has_q_lh[rows, :],
        latitude_deg=f.latitude_deg[rows, :],
    )
    grid = SnowpackGrid(length(rows);
        x=loaded.grid.x,
        y=loaded.grid.y,
        js=js,
        is=is,
        mask=mask,
    )

    println("Applied GrIS mask: kept $(length(rows)) / $(nx * ny) columns")
    return (grid=grid, forcing=forcing)
end

function print_gris_api_help()
    println("Usage:")
    println("  julia --project=. examples/scripts/run_gris_forcing_file_case.jl [options]")
    println()
    println("Core options:")
    println("  --forcing-file=PATH          Prepared NetCDF forcing file")
    println("  --mask-threshold=VALUE       Apply MAR-style MSK threshold before running (default: 50)")
    println("  --output-dir=PATH            Output directory (default: examples/plots/gris_forcing_file_simulation)")
    println("  --netcdf-path=PATH           Output NetCDF path")
    println("  --no-output                  Skip summary/CSV file output")
    println("  --no-nc                      Skip NetCDF output")
    println("  --netcdf-vars=SPEC           NetCDF variables to write")
    println("  --model=NAME                 bessi|pdd (default: bessi)")
    println("  --ntot=N                     Maximum active layers (default: 20)")
    println("  --years=N                    Number of forcing years (default: 10)")
    println("  --backend=threads|cpu|gpu    Execution backend (default: threads)")
    println("  --checkpoint-path=PATH       Write restart checkpoints after completed years")
    println("  --checkpoint-year-stride=N   Checkpoint every N completed years (default: 1)")
    println("  --restart-from=PATH          Resume from a checkpoint and ignore setup options")
    println("  --restart-add-years=N        On restart, run N additional years beyond the checkpoint")
    println("  --restart-target-years=N     On restart, run until total completed years reaches N")
    println("  --restart-override-physics   Rebuild BESSI model physics from explicit CLI/env physics options")
    println()
    println("Physics options:")
    println("  --albedo=NAME                constant|dynamic (default: dynamic)")
    println("  --alpha-dry=VALUE            Dry-snow albedo constant")
    println("  --alpha-wet=VALUE            Wet-snow minimum albedo constant")
    println("  --alpha-ice=VALUE            Bare-ice albedo constant")
    println("  --max-lwc-albedo=VALUE       Liquid-water content scale for albedo darkening")
    println("  --densification=NAME         bessi|htessel (default: bessi)")
    println("  --fresh-snow-density=NAME    constant|parameterized (default: constant)")
    println("  --diurnal-shortwave          Enable adaptive diurnal shortwave substepping for BESSI")
    println("  --diurnal-shortwave-substeps Enable adaptive diurnal shortwave substepping for BESSI")
    println("  --diurnal-shortwave-threshold=VALUE")
    println("                                Peak-minus-mean shortwave threshold for substepping (default: 0)")
    println("  --diurnal-shortwave-max-substeps=N")
    println("                                Maximum diurnal substeps, 1 to 24 (default: 2)")
    println("  --diurnal-shortwave-min-air-temperature-c=VALUE")
    println("                                Minimum daily mean air temperature for substepping (default: -8)")
    println("  --diurnal-temperature-cycle  Enable daily sinusoidal air-temperature cycle inside substeps")
    println("  --diurnal-temperature-amplitude-c=VALUE")
    println("                                Half-amplitude of daily air-temperature cycle in degC (default: 5)")
    println("  --pdd-ddf-snow=VALUE         PDD snow degree-day factor in mmWE d-1 C-1 (default: 3)")
    println("  --pdd-ddf-ice=VALUE          PDD ice degree-day factor in mmWE d-1 C-1 (default: 8)")
    println("  --pdd-refreezing-fraction=X  PDD refreezing fraction (default: 0.6)")
    println()
    println("NetCDF variable groups:")
    println("  final, layers, history, monthly, step")
end

function bool_env_or_flag(args::Vector{String}, flag::String, env_name::String; default::Bool=false)
    has_flag(args, flag) && return true
    value = lowercase(strip(env_value(env_name, default ? "true" : "false")))
    return value in ("1", "true", "yes", "on")
end

function physical_constant_overrides(args::Vector{String})
    pairs = Pair{Symbol, Float64}[]
    for (arg_name, env_name, key) in (
        ("alpha-dry", "ALPHA_DRY", :alpha_dry),
        ("alpha-wet", "ALPHA_WET", :alpha_wet),
        ("alpha-ice", "ALPHA_ICE", :alpha_ice),
        ("max-lwc-albedo", "MAX_LWC_ALBEDO", :max_lwc_albedo),
    )
        value = arg_value(args, arg_name, env_value(env_name, ""))
        isempty(value) && continue
        push!(pairs, key => parse(Float64, value))
    end
    return pairs
end

has_arg_prefix(args::Vector{String}, prefix::AbstractString) = any(arg -> startswith(arg, prefix), args)

function _bessi_physical_constant_pairs(c)
    return Pair{Symbol, Float64}[
        :rho_s => c.rho_s,
        :rho_i => c.rho_i,
        :rho_w => c.rho_w,
        :rho_s_a => c.rho_s_a,
        :rho_s_b => c.rho_s_b,
        :rho_s_c => c.rho_s_c,
        :Ki => c.Ki,
        :ci => c.ci,
        :cw => c.cw,
        :Lm => c.Lm,
        :D_sh => c.D_sh,
        :alpha_dry => c.alpha_dry,
        :alpha_wet => c.alpha_wet,
        :alpha_ice => c.alpha_ice,
        :max_lwc_albedo => c.max_lwc_albedo,
        :ϵ_air => c.ϵ_air,
        :ϵ_snow => c.ϵ_snow,
        :σ => c.σ,
        :R => c.R,
        :T0 => c.T0,
        :seconds_per_day => c.seconds_per_day,
        :seconds_per_month => c.seconds_per_month,
        :seconds_per_year => c.seconds_per_year,
    ]
end

function _merge_pairs(base::Vector{Pair{Symbol, Float64}}, overrides::Vector{Pair{Symbol, Float64}})
    values = Dict{Symbol, Float64}(base)
    for (key, value) in overrides
        values[key] = value
    end
    return Pair{Symbol, Float64}[key => values[key] for (key, _) in base]
end

function _current_albedo_scheme(model::BESSIModel)
    return model.c.albedo_scheme == Chion.ALBEDO_CONSTANT ? ConstantAlbedo() : DynamicAlbedo()
end

function _current_densification_scheme(model::BESSIModel)
    return model.c.low_density_densification == Chion.LOW_DENSIFICATION_HTESSEL ? HTESSELDensification() : BESSIDensification()
end

function _current_fresh_snow_scheme(model::BESSIModel)
    return model.c.fresh_snow_density_scheme == Chion.FRESH_SNOW_DENSITY_PARAMETERIZED ? ParameterizedFreshSnowDensity() : ConstantFreshSnowDensity()
end

function _explicit_or_env(args::Vector{String}, name::AbstractString, env_name::AbstractString)
    has_arg_prefix(args, "--$(name)=") && return arg_value(args, name, "")
    value = env_value(env_name, "")
    return isempty(value) ? "" : value
end

function _update_domain_physics!(domain, model::BESSIModel)
    domain.c = model.c
    domain.mass_max = model.mass_max
    domain.mass_split = model.mass_split
    domain.mass_min = model.mass_min
    domain.rho_max = model.rho_max
    return domain
end

function _restart_physics_transform(args::Vector{String})
    return function (sim)
        model = sim.model
        model isa BESSIModel || error("--restart-override-physics currently supports BESSI checkpoints only.")

        ntot_value = _explicit_or_env(args, "ntot", "NTOT")
        if !isempty(ntot_value) && parse(Int, ntot_value) != model.Ntot
            error("Cannot override Ntot on restart: checkpoint has Ntot=$(model.Ntot), requested $(ntot_value). Start a fresh run to change layer count.")
        end

        albedo_value = _explicit_or_env(args, "albedo", "ALBEDO")
        densification_value = _explicit_or_env(args, "densification", "DENSIFICATION")
        fresh_snow_value = _explicit_or_env(args, "fresh-snow-density", "FRESH_SNOW_DENSITY")
        diurnal_shortwave_set = has_flag(args, "diurnal-shortwave") || has_flag(args, "diurnal-shortwave-substeps") ||
            lowercase(strip(env_value("DIURNAL_SHORTWAVE", ""))) in ("1", "true", "yes", "on") ||
            lowercase(strip(env_value("DIURNAL_SHORTWAVE_SUBSTEPS", ""))) in ("1", "true", "yes", "on")
        diurnal_temperature_set = has_flag(args, "diurnal-temperature-cycle") ||
            lowercase(strip(env_value("DIURNAL_TEMPERATURE_CYCLE", ""))) in ("1", "true", "yes", "on")

        diurnal_threshold_value = _explicit_or_env(args, "diurnal-shortwave-threshold", "DIURNAL_SHORTWAVE_THRESHOLD")
        diurnal_max_substeps_value = _explicit_or_env(args, "diurnal-shortwave-max-substeps", "DIURNAL_SHORTWAVE_MAX_SUBSTEPS")
        diurnal_min_temperature_value = _explicit_or_env(args, "diurnal-shortwave-min-air-temperature-c", "DIURNAL_SHORTWAVE_MIN_AIR_TEMPERATURE_C")
        diurnal_amplitude_value = _explicit_or_env(args, "diurnal-temperature-amplitude-c", "DIURNAL_TEMPERATURE_AMPLITUDE_C")

        physical_kwargs = _merge_pairs(_bessi_physical_constant_pairs(model.c), physical_constant_overrides(args))
        new_model = BESSIModel(
            model.grid;
            albedo=isempty(albedo_value) ? _current_albedo_scheme(model) : _albedo_scheme(lowercase(albedo_value)),
            densification=isempty(densification_value) ? _current_densification_scheme(model) : _densification_scheme(lowercase(densification_value)),
            fresh_snow_density=isempty(fresh_snow_value) ? _current_fresh_snow_scheme(model) : _fresh_snow_scheme(lowercase(fresh_snow_value)),
            Ntot=model.Ntot,
            mass_max=model.mass_max,
            mass_split=model.mass_split,
            mass_min=model.mass_min,
            rho_max=model.rho_max,
            density_init=model.density_init,
            temperature_init=model.temperature_init,
            diurnal_shortwave_substeps=diurnal_shortwave_set ? true : model.diurnal_shortwave_substeps,
            diurnal_shortwave_threshold=isempty(diurnal_threshold_value) ? model.diurnal_shortwave_threshold : parse(Float64, diurnal_threshold_value),
            diurnal_shortwave_max_substeps=isempty(diurnal_max_substeps_value) ? model.diurnal_shortwave_max_substeps : parse(Int, diurnal_max_substeps_value),
            diurnal_shortwave_min_air_temperature_c=isempty(diurnal_min_temperature_value) ? model.diurnal_shortwave_min_air_temperature - 273.15 : parse(Float64, diurnal_min_temperature_value),
            diurnal_temperature_cycle=diurnal_temperature_set ? true : model.diurnal_temperature_cycle,
            diurnal_temperature_amplitude_c=isempty(diurnal_amplitude_value) ? model.diurnal_temperature_amplitude : parse(Float64, diurnal_amplitude_value),
            physical_kwargs...,
        )

        sim.model = new_model
        _update_domain_physics!(sim.ref.domain, new_model)
        _update_domain_physics!(sim.now.domain, new_model)
        println("Restart physics override: rebuilt BESSI model from explicit CLI/env physics options.")
        println(
            "Restart model physics: albedo_scheme=$(new_model.c.albedo_scheme), " *
            "densification=$(new_model.c.low_density_densification), " *
            "fresh_snow=$(new_model.c.fresh_snow_density_scheme), " *
            "diurnal_shortwave=$(new_model.diurnal_shortwave_substeps), " *
            "diurnal_temperature_cycle=$(new_model.diurnal_temperature_cycle)",
        )
        return sim
    end
end

function _maybe_arg_or_env(args::Vector{String}, name::AbstractString, env_name::AbstractString)
    has_arg_prefix(args, "--$(name)=") && return arg_value(args, name, "")
    return env_value(env_name, "")
end

function _restart_options_transform(args::Vector{String})
    return function (options::Chion.RunOptions, checkpoint)
        completed_years = checkpoint.stepper.completed_years
        years_value = _maybe_arg_or_env(args, "years", "YEARS")
        add_years_value = arg_value(args, "restart-add-years", env_value("RESTART_ADD_YEARS", ""))
        target_years_value = arg_value(args, "restart-target-years", env_value("RESTART_TARGET_YEARS", ""))

        target_sources = count(!isempty, (years_value, target_years_value))
        if !isempty(add_years_value) && target_sources > 0
            error("Use either --restart-add-years=N or a target year option (--years=N / --restart-target-years=N), not both.")
        end

        years = options.years
        if !isempty(add_years_value)
            years = completed_years + parse(Int, add_years_value)
        elseif !isempty(target_years_value)
            years = parse(Int, target_years_value)
        elseif !isempty(years_value)
            years = parse(Int, years_value)
        end
        years > completed_years || error(
            "Restart target years=$(years) is not greater than checkpoint completed years=$(completed_years). " *
            "For another 200 years from this checkpoint, use --restart-add-years=200 or --years=$(completed_years + 200).",
        )

        output_dir = arg_value(args, "output-dir", options.output_dir)
        netcdf_path = arg_value(args, "netcdf-path", options.netcdf_path)
        write_outputs = has_flag(args, "no-output") ? false : options.write_outputs
        write_netcdf = has_flag(args, "no-nc") ? false : options.write_netcdf
        netcdf_variables = write_netcdf ?
            Chion.normalize_netcdf_variables(arg_value(args, "netcdf-vars", join(String.(options.netcdf_variables), ","))) :
            Symbol[]
        backend = arg_value(args, "backend", string(options.backend))
        history_year_stride = parse(Int, arg_value(args, "history-year-stride", string(options.history_year_stride)))

        if years != options.years
            println(
                "Restart year target override: checkpoint completed=$(completed_years), " *
                "checkpoint target=$(options.years), new target=$(years), additional=$(years - completed_years)",
            )
        end

        return Chion.RunOptions(
            name=options.name,
            input_label=options.input_label,
            output_dir=output_dir,
            netcdf_path=netcdf_path,
            write_outputs=write_outputs,
            write_netcdf=write_netcdf,
            netcdf_variables=netcdf_variables,
            years=years,
            backend=backend,
            history_year_stride=history_year_stride,
        )
    end
end

const RESTART_PHYSICS_ARG_NAMES = (
    "albedo",
    "alpha-dry",
    "alpha-wet",
    "alpha-ice",
    "max-lwc-albedo",
    "densification",
    "fresh-snow-density",
    "ntot",
    "diurnal-shortwave",
    "diurnal-shortwave-substeps",
    "diurnal-shortwave-threshold",
    "diurnal-shortwave-max-substeps",
    "diurnal-shortwave-min-air-temperature-c",
    "diurnal-temperature-cycle",
    "diurnal-temperature-amplitude-c",
)

function _restart_has_physics_args(args::Vector{String})
    return any(RESTART_PHYSICS_ARG_NAMES) do name
        has_flag(args, name) || has_arg_prefix(args, "--$(name)=")
    end
end

function main(args::Vector{String})
    if has_flag(args, "help")
        print_gris_api_help()
        return
    end

    restart_from = arg_value(args, "restart-from", "")
    if !isempty(restart_from)
        checkpoint_path = arg_value(args, "checkpoint-path", "")
        checkpoint_year_stride = parse(Int, arg_value(args, "checkpoint-year-stride", "1"))
        if !isempty(checkpoint_path) && abspath(checkpoint_path) == abspath(restart_from)
            println("Restart warning: --checkpoint-path equals --restart-from; this run will overwrite the input checkpoint as it progresses.")
        end
        override_physics = has_flag(args, "restart-override-physics")
        if !override_physics && _restart_has_physics_args(args)
            println(
                "Restart warning: physics options were provided, but checkpointed physics will be used. " *
                "Pass --restart-override-physics to rebuild BESSI physics from CLI/env options.",
            )
        end
        sim_transform = override_physics ? _restart_physics_transform(args) : identity
        integrator = restart_integrator(
            restart_from;
            io=stdout,
            sim_transform=sim_transform,
            options_transform=_restart_options_transform(args),
        )
        run!(
            integrator;
            checkpoint_path=checkpoint_path,
            checkpoint_year_stride=checkpoint_year_stride,
        )
        result = finalize!(integrator)
        println("Status: $(result.status)")
        return
    end

    default_forcing_file = default_gris_forcing_file_path()
    forcing_file = arg_value(args, "forcing-file", default_forcing_file)
    isempty(forcing_file) && error("Pass --forcing-file=PATH or set FORCING_PATH to the prepared forcing file.")

    loaded = load_forcing_file(forcing_file; air_temperature_name="TTZ")
    loaded = apply_gris_mask(
        loaded,
        forcing_file;
        threshold=parse(Float64, arg_value(args, "mask-threshold", "50.0")),
    )
    model_name = lowercase(arg_value(args, "model", env_value("MODEL", "bessi")))
    model = if model_name == "bessi"
        physical_kwargs = physical_constant_overrides(args)
        diurnal_shortwave = bool_env_or_flag(args, "diurnal-shortwave", "DIURNAL_SHORTWAVE")
        diurnal_shortwave_substeps = bool_env_or_flag(
            args,
            "diurnal-shortwave-substeps",
            "DIURNAL_SHORTWAVE_SUBSTEPS";
            default=diurnal_shortwave,
        )
        diurnal_temperature_cycle = bool_env_or_flag(
            args,
            "diurnal-temperature-cycle",
            "DIURNAL_TEMPERATURE_CYCLE",
        )
        build_model(
            :bessi,
            loaded.grid;
            Ntot=parse(Int, arg_value(args, "ntot", "20")),
            albedo=_albedo_scheme(lowercase(arg_value(args, "albedo", "dynamic"))),
            densification=_densification_scheme(lowercase(arg_value(args, "densification", "bessi"))),
            fresh_snow_density=_fresh_snow_scheme(lowercase(arg_value(args, "fresh-snow-density", "constant"))),
            diurnal_shortwave=diurnal_shortwave,
            diurnal_shortwave_substeps=diurnal_shortwave_substeps,
            diurnal_shortwave_threshold=parse(Float64, arg_value(args, "diurnal-shortwave-threshold", env_value("DIURNAL_SHORTWAVE_THRESHOLD", "0.0"))),
            diurnal_shortwave_max_substeps=parse(Int, arg_value(args, "diurnal-shortwave-max-substeps", env_value("DIURNAL_SHORTWAVE_MAX_SUBSTEPS", "2"))),
            diurnal_shortwave_min_air_temperature_c=parse(Float64, arg_value(args, "diurnal-shortwave-min-air-temperature-c", env_value("DIURNAL_SHORTWAVE_MIN_AIR_TEMPERATURE_C", "-8.0"))),
            diurnal_temperature_cycle=diurnal_temperature_cycle,
            diurnal_temperature_amplitude_c=parse(Float64, arg_value(args, "diurnal-temperature-amplitude-c", env_value("DIURNAL_TEMPERATURE_AMPLITUDE_C", "5.0"))),
            physical_kwargs...,
        )
    elseif model_name == "pdd"
        build_model(
            :pdd,
            loaded.grid;
            ddf_snow=parse(Float64, arg_value(args, "pdd-ddf-snow", env_value("PDD_DDF_SNOW", "3.0"))),
            ddf_ice=parse(Float64, arg_value(args, "pdd-ddf-ice", env_value("PDD_DDF_ICE", "8.0"))),
            refreezing_fraction=parse(Float64, arg_value(args, "pdd-refreezing-fraction", env_value("PDD_REFREEZING_FRACTION", "0.6"))),
        )
    else
        error("Unsupported --model=$(model_name). Use bessi or pdd.")
    end
    println("Model: $(model_name)")
    if model_name == "bessi"
        isempty(physical_kwargs) || println("Physical constant overrides: ", join(["$(k)=$(v)" for (k, v) in physical_kwargs], ", "))
        println(
            "Diurnal shortwave substeps: $(model.diurnal_shortwave_substeps) " *
            "(legacy flag requested: $(diurnal_shortwave), " *
            "threshold: $(model.diurnal_shortwave_threshold), " *
            "max substeps: $(model.diurnal_shortwave_max_substeps), " *
            "min air temp: $(model.diurnal_shortwave_min_air_temperature - 273.15) degC, " *
            "temperature cycle: $(model.diurnal_temperature_cycle), " *
            "temperature amplitude: $(model.diurnal_temperature_amplitude) degC)",
        )
    end
    write_netcdf = !has_flag(args, "no-nc")
    default_netcdf_vars = "all"
    netcdf_vars = write_netcdf ? arg_value(args, "netcdf-vars", default_netcdf_vars) : Symbol[]
    println("NetCDF variables: ", isempty(netcdf_vars) ? "(none)" : netcdf_vars)
    output_dir = arg_value(args, "output-dir", DEFAULT_GRIS_FORCING_FILE_OUTPUT_DIR)
    checkpoint_path = arg_value(args, "checkpoint-path", "")
    checkpoint_year_stride = parse(Int, arg_value(args, "checkpoint-year-stride", "1"))
    sim = Simulation(
        model;
        forcing=loaded.forcing,
        netcdf_variables=netcdf_vars,
        write_netcdf=write_netcdf,
        output_dir=output_dir,
        netcdf_path=arg_value(args, "netcdf-path", ""),
        write_outputs=!has_flag(args, "no-output"),
        years=parse(Int, arg_value(args, "years", "10")),
        backend=arg_value(args, "backend", "threads"),
    )
    result = run!(
        sim;
        io=stdout,
        checkpoint_path=checkpoint_path,
        checkpoint_year_stride=checkpoint_year_stride,
    )
    println("Status: $(result.status)")
    return
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
