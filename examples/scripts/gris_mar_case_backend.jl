#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Chion

const SM = Chion.SnowpackModel
const DEFAULT_GRIS_MAR_NC_PATH = let path = ""
    for candidate in (
        "/p/projects/ou/labs/ai/Nils/MARv3.14.3-10km-daily-ERA5-2025.nc",
        "/Users/niboch001/Downloads/MARv3.14.3-10km-daily-ERA5-2026.nc",
    )
        if isfile(candidate)
            path = candidate
            break
        end
    end
    path
end
const DEFAULT_GRIS_MAR_OUTPUT_DIR = joinpath(@__DIR__, "..", "plots", "gris_mar_case")

function load_gris_mar_problem(
    nc_path::AbstractString;
    mask_threshold::Float64=50.0,
    turbulent_flux_sign::Float64=1.0,
    ntot::Integer=20,
    physics::SM.SnowpackPhysicalConstants{Float64}=Chion.physics(),
)
    definition = Chion.mar_case(
        nc_path;
        mask_threshold=mask_threshold,
        turbulent_flux_sign=turbulent_flux_sign,
        physics=physics,
        ntot=ntot,
    )
    return (
        domain=definition.domain,
        forcing=definition.forcing,
        layout=definition.layout,
        wind_forcing_message=isempty(definition.notes) ? "" : definition.notes[1],
        nvalid=hasproperty(definition.metadata, :ncol) ? definition.metadata.ncol : SM.column_count(definition.domain),
        mask_threshold=mask_threshold,
        nc_path=String(nc_path),
    )
end

function run_gris_mar_case(
    nc_path::AbstractString;
    io::IO=stdout,
    mask_threshold::Float64=50.0,
    turbulent_flux_sign::Float64=1.0,
    ntot::Integer=20,
    run::Chion.RunConfig=Chion.RunConfig(),
    physics::SM.SnowpackPhysicalConstants{Float64}=Chion.physics(),
)
    definition = Chion.mar_case(
        nc_path;
        mask_threshold=mask_threshold,
        turbulent_flux_sign=turbulent_flux_sign,
        physics=physics,
        ntot=ntot,
    )
    case = Chion.build_case(
        definition;
        run=Chion.RunConfig(
            name=isempty(run.name) ? "GrIS MAR case" : run.name,
            input_label=isempty(run.input_label) ? abspath(nc_path) : run.input_label,
            output_dir=run.output_dir,
            netcdf_path=run.netcdf_path,
            write_outputs=run.write_outputs,
            write_netcdf=run.write_netcdf,
            netcdf_variables=run.netcdf_variables,
            cycles=run.cycles,
            history_stride=run.history_stride,
            backend=run.backend,
        ),
    )
    return Chion.run_case(case; io=io)
end
