#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Chion

env_bool(name, default) = lowercase(get(ENV, name, string(default))) in ("1", "true", "yes")

# Edit this block to configure the run.
const CONFIG = (
    forcing_file="/p/projects/ou/labs/ai/Nils/MAR3.14/MARv3.14.3-10km-daily-ERA5-1940-1980_daily_climatology.nc",
    model=:bessi,
    output_file=get(ENV, "CHION_OUTPUT_FILE", joinpath(@__DIR__, "..", "plots", "gris_itm_monthly.nc")),
    mask_threshold=50.0,
    years=parse(Int, get(ENV, "CHION_YEARS", "100")),
    backend=:gpu,
    write_netcdf=env_bool("CHION_WRITE_NETCDF", true),
    netcdf_variables=:monthly,
)

const MODEL_OPTIONS = (
    Ntot=20,
    albedo=:aging,
    densification=:bessi,
    fresh_snow_density=:constant,
    diurnal_shortwave_substeps=true,
    diurnal_shortwave_threshold=0.0,
    diurnal_shortwave_max_substeps=4,
)

function main()
    # 1. Load the valid Greenland columns from the prepared MAR forcing file.
    loaded = load_forcing_file(
        CONFIG.forcing_file;
        time_name="TIME",
        air_temperature_name="TTZ",
        mask_name="MSK",
        mask_threshold=CONFIG.mask_threshold,
    )
    grid, forcing = loaded.grid, loaded.forcing

    if CONFIG.model == :itm
        # The MAR file supplies surface elevation but not ice thickness or
        # annual PDD. The mask selects ice-sheet columns, so a positive ice
        # thickness explicitly selects ITM's ice branch; annual PDD is
        # diagnosed from this daily climatology and then held fixed per year.
        ncol, ntime = size(forcing.air_temperature)
        annual_pdd = vec(sum(max.(forcing.air_temperature .- 273.15, 0.0) .* reshape(forcing.dt_days, 1, ntime); dims=2))
        forcing = SnowpackForcing(
            time_values=forcing.time_values,
            dt_days=forcing.dt_days,
            air_temperature=forcing.air_temperature,
            snowfall_rate=forcing.snowfall_rate,
            rainfall_rate=forcing.rainfall_rate,
            shortwave_down=forcing.shortwave_down,
            latitude_deg=forcing.latitude_deg,
            surface_height=forcing.surface_height,
            ice_thickness=ones(ncol, ntime),
            annual_pdd=repeat(annual_pdd, 1, ntime),
        )
        model = ITMModel(grid)
    elseif CONFIG.model == :pdd
        model = PDDModel(grid; pdd_method=:simple, temperature_sigma=5.0)
    elseif CONFIG.model == :bessi
        model = BESSIModel(grid; MODEL_OPTIONS...)
    else
        error("Unsupported CONFIG.model=$(CONFIG.model). Use :bessi, :pdd, or :itm.")
    end
    # 3. Configure the simulation.
    simulation = Simulation(
        model;
        forcing=forcing,
        years=CONFIG.years,
        backend=CONFIG.backend,
        write_netcdf=CONFIG.write_netcdf,
        netcdf_variables=CONFIG.netcdf_variables,
        netcdf_path=CONFIG.output_file,
        # GPU year-summary reduction currently faults on this full spatial
        # domain. Monthly output is independent of these optional summaries.
        compute_year_metrics=false,
        name="gris_example",
    )

    # 4. Run and inspect simulation.now for the final model state.
    result = run!(simulation)
    println("Status: ", result.status)
    println("Output: ", result.netcdf_path)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
