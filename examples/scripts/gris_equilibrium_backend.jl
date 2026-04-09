#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Chion

const SM = Chion.SnowpackModel
const DEFAULT_GRIS_API_NC_PATH = let path = ""
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
const DEFAULT_OUT_DIR_EQUIL_API = joinpath(@__DIR__, "..", "plots", "gris_equilibrium_api")

function load_gris_equilibrium_problem(
    nc_path::AbstractString;
    mask_threshold::Float64=50.0,
    turbulent_flux_sign::Float64=1.0,
    ntot::Integer=20,
    physics::SM.SnowpackPhysicalConstants{Float64}=SM.SnowpackPhysicalConstants(),
)
    data = Chion.load_forcing(
        Chion.mar_forcing(
            nc_path;
            mask_threshold=mask_threshold,
            turbulent_flux_sign=turbulent_flux_sign,
        );
        physics=physics,
        ntot=ntot,
    )
    return (
        domain=data.domain,
        forcing=data.forcing,
        layout=data.layout,
        wind_forcing_message=isempty(data.notes) ? "" : data.notes[1],
        nvalid=hasproperty(data.metadata, :ncol) ? data.metadata.ncol : SM.column_count(data.domain),
        mask_threshold=mask_threshold,
        nc_path=String(nc_path),
    )
end

function run_gris_equilibrium_from_file(
    nc_path::AbstractString;
    io::IO=stdout,
    mask_threshold::Float64=50.0,
    turbulent_flux_sign::Float64=1.0,
    ntot::Integer=20,
    options::Chion.EquilibriumRunOptions=Chion.EquilibriumRunOptions(),
    physics::SM.SnowpackPhysicalConstants{Float64}=SM.SnowpackPhysicalConstants(),
)
    case = Chion.build_case(
        name=isempty(options.name) ? "GrIS equilibrium spin-up" : options.name,
        forcing=Chion.mar_forcing(
            nc_path;
            mask_threshold=mask_threshold,
            turbulent_flux_sign=turbulent_flux_sign,
        ),
        physics=physics,
        ntot=ntot,
        backend=options.backend,
        forcing_label=isempty(options.forcing_label) ? abspath(nc_path) : options.forcing_label,
        out_dir=options.out_dir,
        out_nc=options.out_nc,
        write_outputs=options.write_outputs,
        write_netcdf=options.write_netcdf,
        netcdf_variables=options.netcdf_variables,
        max_cycles=options.max_cycles,
        cycle_metrics_stride=options.cycle_metrics_stride,
    )
    return Chion.run_case(case; io=io)
end
