#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Chion

# Edit this block to configure the run.
const CONFIG = (
    forcing_file="/p/projects/ou/labs/ai/Nils/MAR3.14/MARv3.14.3-10km-daily-ERA5-1940-1980_daily_climatology.nc",
    output_file=joinpath(@__DIR__, "..", "plots", "gris_example.nc"),
    mask_threshold=50.0,
    years=100,
    backend=:threads,
    write_netcdf=true,
    netcdf_variables=:monthly,
)

const MODEL_OPTIONS = (
    Ntot=20,
    albedo=:dynamic,
    densification=:bessi,
    fresh_snow_density=:constant,
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

    # 2. Configure the physical model.
    model = BESSIModel(grid; MODEL_OPTIONS...)

    # 3. Configure the simulation.
    simulation = Simulation(
        model;
        forcing=forcing,
        years=CONFIG.years,
        backend=CONFIG.backend,
        write_netcdf=CONFIG.write_netcdf,
        netcdf_variables=CONFIG.netcdf_variables,
        netcdf_path=CONFIG.output_file,
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
