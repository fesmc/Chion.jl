#!/usr/bin/env julia

"""Render the transient MAR-forced Greenland monthly BESSI SMB as an MP4.

The daily GrIS run includes a repeated 500-year spin-up before the actual
1940--1980 forcing sequence.  This animation uses only that final, transient
forcing segment (492 monthly records).
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "visualization"))

using CairoMakie
using Dates
using GeoJSON
using NCDatasets
using Proj
using Statistics

const CONFIG = (
    input_file=get(ENV, "CHION_ANIMATION_INPUT_FILE", joinpath(@__DIR__, "..", "plots", "gris_mar_daily_monthly_cpu.nc")),
    output_file=get(ENV, "CHION_ANIMATION_OUTPUT_FILE", joinpath(@__DIR__, "..", "plots", "gris_monthly_smb_1940_1980.mp4")),
    coastline_file=joinpath(@__DIR__, "..", "data", "ne_50m_coastline.geojson"),
    first_year=parse(Int, get(ENV, "CHION_ANIMATION_FIRST_YEAR", "1940")),
    last_year=parse(Int, get(ENV, "CHION_ANIMATION_LAST_YEAR", "1980")),
    all_records=lowercase(get(ENV, "CHION_ANIMATION_ALL_RECORDS", "false")) in ("1", "true", "yes"),
    variable=Symbol(get(ENV, "CHION_ANIMATION_VARIABLE", "smb_ice")),
    framerate=parse(Int, get(ENV, "CHION_ANIMATION_FRAMERATE", "12")),
    smb_limit=parse(Float64, get(ENV, "CHION_ANIMATION_SMB_LIMIT", "400.0")),
)

frame_title(index) = begin
    year = CONFIG.first_year + (index - 1) ÷ 12
    month = mod(index - 1, 12) + 1
    date = Dates.format(Date(year, month, 1), "mm - yyyy")
    return "$date"
end

"""Return Natural Earth coastline points in the plotted (x, y) grid coordinates.

The MAR grid is EPSG:3413 in kilometres, with x (easting) horizontal and y
(northing) vertical.  `Point2f(NaN, NaN)` separates individual coastline
segments.
"""
function projected_coastline(path)
    isfile(path) || error("Coastline file does not exist: $path")
    projection = Proj.Transformation("EPSG:4326", "EPSG:3413"; always_xy=true)
    points = Point2f[]

    function append_segment!(coordinates)
        for (lon, lat) in coordinates
            easting_m, northing_m = projection((Float64(lon), Float64(lat)))
            if isfinite(easting_m) && isfinite(northing_m)
                # Natural Earth / Proj coordinates are metres; MAR uses km.
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
            # The Natural Earth file is mostly LineStrings, with one
            # MultiLineString that needs a separator between each component.
            foreach(append_segment!, coordinates)
        end
    end
    return points
end

function finite_values(field)
    values = Float64.(coalesce.(field, NaN))
    return values[isfinite.(values)]
end

function colour_limit(smb, first_record, nframes)
    sample = Float64[]
    for record in first_record:(first_record + nframes - 1)
        append!(sample, finite_values(smb[:, :, record]))
    end
    isempty(sample) && error("No finite SMB values found in the selected forcing period.")
    return max(1.0, quantile(abs.(sample), 0.99))
end

function monthly_field(smb, record)
    # NetCDF storage is (x, y, t), matching Makie's x-horizontal, y-vertical
    # grid convention for contourf.
    return Float32.(coalesce.(smb[:, :, record], NaN))
end

function main()
    isfile(CONFIG.input_file) || error("Input file does not exist: $(CONFIG.input_file)")
    mkpath(dirname(CONFIG.output_file))

    NCDataset(CONFIG.input_file) do dataset
        variable_name = String(CONFIG.variable)
        haskey(dataset, variable_name) || error("Input file has no `$variable_name` variable.")
        smb = dataset[variable_name]
        x = Float64.(dataset["x"][:])
        y = Float64.(dataset["y"][:])
        nrecords = size(smb, 3)
        expected_frames = (CONFIG.last_year - CONFIG.first_year + 1) * 12
        nframes = CONFIG.all_records ? nrecords : expected_frames
        nrecords >= nframes || error("Input has $nrecords monthly records, need $nframes.")
        first_record = CONFIG.all_records ? 1 : nrecords - nframes + 1
        vmax = CONFIG.smb_limit
        coastline = projected_coastline(CONFIG.coastline_file)

        CairoMakie.activate!()
        first_frame = Observable(monthly_field(smb, first_record))
        title = Observable(frame_title(1))
        figure = Figure(
            size=(780, 900),
            fontsize=22,
            figure_padding=(8, 8, 10, 8),
            backgroundcolor=:white,
        )
        colgap!(figure.layout, 6)
        axis = Axis(
            figure[1, 1];
            title=title,
            titlegap=5,
            aspect=DataAspect(),
            backgroundcolor=:white,
            limits=(minimum(x), maximum(x), minimum(y), maximum(y)),
        )
        hidedecorations!(axis)
        hidespines!(axis)
        levels = range(-vmax, vmax; length=21)
        filled_contours = contourf!(
            axis,
            x,
            y,
            first_frame;
            levels,
            colormap=cgrad([:firebrick3, :white, :dodgerblue3]),
            nan_color=:transparent,
        )
        lines!(axis, coastline; color=:black, linewidth=2.5)
        Colorbar(
            figure[1, 2],
            filled_contours;
            label=CONFIG.variable == :surface_smb ? "Surface Mass Balance" : "Ice SMB (mm w.e. month⁻¹)",
            ticks=-400:200:400,
            width=30,
            height=Relative(0.70),
            valign=:center,
        )

        record(figure, CONFIG.output_file, 1:nframes; framerate=CONFIG.framerate) do frame
            first_frame[] = monthly_field(smb, first_record + frame - 1)
            title[] = frame_title(frame)
        end
        println("Animation: $(CONFIG.output_file)")
        println("Frames: $nframes; fixed colour limit: ±$(round(vmax; digits=1)) mm w.e. month⁻¹")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
