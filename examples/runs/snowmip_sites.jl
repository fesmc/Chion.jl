#!/usr/bin/env julia

# The ten ESM-SnowMIP reference sites (Ménard et al., 2019) with their in situ forcing, in the land
# configuration of Chion (snow on soil without a glacier-ice substrate, HTESSEL densification,
# parameterised fresh-snow density; the forcing is already hourly, so no diurnal substeps).
# Writes output/snowmip/runs/<site>.csv: SWE, depth and albedo after every forcing step.
# Figures and statistics: examples/figures/snowmip_figures.jl.
#
#   julia --project=. --threads 10 examples/runs/snowmip_sites.jl
using Chion

const DATA = joinpath(@__DIR__, "..", "..", "data", "ESM-SnowMIP_all")
const OUT = joinpath(@__DIR__, "..", "..", "output", "snowmip", "runs")

sites = [only(match(r"^met_insitu_([a-z]+)_", f).captures) => joinpath(DATA, f)
         for f in readdir(DATA) if startswith(f, "met_insitu_")]
Threads.@threads for (site, path) in sites
    forcing = Forcing(path; format=:alma)
    model = BESSIModel(forcing; Ntot=60, densification=:htessel, fresh_snow_density=:parameterized,
        diurnal_shortwave_substeps=false, diurnal_temperature_cycle=false, ice_substrate_layers=0)
    run!(Experiment(model, forcing; name=site, verbose=false, outputs=Output(joinpath(OUT, "$site.csv");
        variables=[:snow_mass, :thickness, :albedo], frequency=:step, kind=:snapshot)))
end
