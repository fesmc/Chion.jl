#!/usr/bin/env julia
"""Add a missing MAR SWU field as AL2 × SWD without loading a full year at once."""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using NCDatasets

function main(path)
    NCDataset(path, "a") do dataset
        haskey(dataset, "SWU") && return
        albedo = dataset["AL2"]
        downward = dataset["SWD"]
        swu = defVar(
            dataset,
            "SWU",
            Float32,
            dimnames(downward);
            attrib=Dict(
                "units" => "W/m2",
                "long_name" => "Derived upwelling shortwave (AL2 × SWD)",
            ),
        )
        ntime = size(downward, 3)
        for start_index in 1:31:ntime
            stop_index = min(start_index + 30, ntime)
            time_range = start_index:stop_index
            # MAR AL2 has two sectors; the forcing loader consistently uses
            # sector 1 as the surface albedo.
            swu[:, :, time_range] = albedo[:, :, 1, time_range] .* downward[:, :, time_range]
        end
    end
end

length(ARGS) == 1 || error("usage: add_mar_derived_swu.jl FILE")
main(only(ARGS))
