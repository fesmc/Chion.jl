#!/usr/bin/env julia

"""Render the MAR-forced Antarctic monthly BESSI SMB (aging albedo) as an MP4.

The input is the completed max-snow Antarctic experiment using ``albedo =
:aging``.  Every available monthly record is included, with timestamps read
directly from the NetCDF file.
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "visualization"))

using CairoMakie
using Dates
using GeoJSON
using NCDatasets
using Proj

const CONFIG = (
    input_file=joinpath(@__DIR__, "..", "plots", "antarctica_mar_max_snow_10cycles_monthly_gpu.nc"),
    output_file=joinpath(@__DIR__, "..", "plots", "antarctica_monthly_smb_aging_1979_2026.mp4"),
    coastline_file=joinpath(@__DIR__, "..", "data", "ne_50m_coastline.geojson"),
    framerate=12,
    smb_limit=400.0,
)

frame_title(time) = Dates.format(Date(time), "mm - yyyy")

"""Return Natural Earth coastlines in the Antarctic MAR grid coordinates.

MAR Antarctic x/y coordinates are south-polar stereographic EPSG:3031 in km.
NaN points delimit individual coastline segments for Makie.
"""
function projected_coastline(path)
    isfile(path) || error("Coastline file does not exist: $path")
    projection = Proj.Transformation("EPSG:4326", "EPSG:3031"; always_xy=true)
    points = Point2f[]

    function append_segment!(coordinates)
        for (lon, lat) in coordinates
            easting_m, northing_m = projection((Float64(lon), Float64(lat)))
            if isfinite(easting_m) && isfinite(northing_m)
                push!(points, Point2f(easting_m / 1000, northing_m / 1000))
            end
        end
        push!(points, Point2f(NaN, NaN))
    end

    for feature in GeoJSON.read(path).features
        coordinates = feature.geometry.coordinates
        if first(coordinates) isa Tuple
            append_segment!(coordinates)
        else
            foreach(append_segment!, coordinates)
        end
    end
    return points
end

function monthly_field(smb, record)
    # NetCDF storage is (x, y, t), matching Makie's x-horizontal convention.
    return Float32.(coalesce.(smb[:, :, record], NaN))
end

function main()
    isfile(CONFIG.input_file) || error("Input file does not exist: $(CONFIG.input_file)")
    mkpath(dirname(CONFIG.output_file))

    NCDataset(CONFIG.input_file) do dataset
        haskey(dataset, "smb_ice") || error("Input file has no `smb_ice` variable.")
        smb = dataset["smb_ice"]
        x = Float64.(dataset["x"][:])
        y = Float64.(dataset["y"][:])
        times = DateTime.(dataset["t"][:])
        nframes = length(times)
        size(smb, 3) == nframes || error("SMB and time dimensions disagree.")
        coastline = projected_coastline(CONFIG.coastline_file)
        vmax = CONFIG.smb_limit

        CairoMakie.activate!()
        first_frame = Observable(monthly_field(smb, 1))
        title = Observable(frame_title(times[1]))
        figure = Figure(
            size=(900, 900),
            fontsize=22,
            figure_padding=(8, 8, 10, 8),
            backgroundcolor=:white,
        )
        colgap!(figure.layout, 6)
        axis = Axis(
            figure[1, 1];
            title,
            titlegap=5,
            aspect=DataAspect(),
            backgroundcolor=:white,
            limits=(minimum(x), maximum(x), minimum(y), maximum(y)),
        )
        hidedecorations!(axis)
        hidespines!(axis)
        levels = range(-vmax, vmax; length=21)
        filled_contours = contourf!(
            axis, x, y, first_frame;
            levels,
            colormap=cgrad([:firebrick3, :white, :dodgerblue3]),
            nan_color=:transparent,
        )
        lines!(axis, coastline; color=:black, linewidth=2.0)
        Colorbar(
            figure[1, 2], filled_contours;
            label="SMB (mm w.e. month⁻¹)",
            ticks=-400:200:400,
            width=30,
            height=Relative(0.70),
            valign=:center,
        )

        record(figure, CONFIG.output_file, 1:nframes; framerate=CONFIG.framerate) do frame
            first_frame[] = monthly_field(smb, frame)
            title[] = frame_title(times[frame])
        end
        println("Animation: $(CONFIG.output_file)")
        println("Frames: $nframes ($(frame_title(first(times))) to $(frame_title(last(times)))); fixed colour limit: ±$(vmax) mm w.e. month⁻¹")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
