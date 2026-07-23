#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Printf
using Chion
using Plots

const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "..", "plots", "diurnal_shortwave")

# Present-day (0 ka) row of ZB18a(1,1), from Kocken and Zeebe's
# paleoinsolation data set: https://github.com/japhir/paleoinsolation.
# The source longitude of perihelion is shifted by pi before use, matching
# the conversion in src/paleoinsolation.f90 in that repository.
const ZB18A_1_1_PRESENT_DAY = (
    eccentricity=1.670545044954422e-2,
    obliquity_rad=4.090928042223287e-1,
    lpx_rad=mod(1.796246057579526 - pi, 2pi),
)

arg_value(args, name, default) = begin
    prefix = "--$(name)="
    for arg in args
        startswith(arg, prefix) && return arg[length(prefix)+1:end]
    end
    return default
end

has_flag(args, name) = any(==("--$(name)"), args)

parse_list(txt, ::Type{T}) where {T} =
    T[parse(T, strip(part)) for part in split(txt, ",") if !isempty(strip(part))]

function print_help()
    println("Usage:")
    println("  julia --project=. examples/scripts/plot_diurnal_shortwave_cycle.jl [options]")
    println()
    println("Options:")
    println("  --latitudes=LIST             Comma-separated degrees north, one PDF panel each (default: 40,50,60,70)")
    println("  --solar-longitudes=LIST      Comma-separated degrees, e.g. 0,90,180,270 (default: 90)")
    println("  --shortwave-mean=VALUE       Target daily mean W m^-2 for normalized curves (default: 200)")
    println("  --solar-constant=VALUE       Paleo solar constant W m^-2 (default: 1360.7)")
    println("  --samples=N                  Samples across one day (default: 289)")
    println("  --substeps=LIST              Substep-average curves to plot (default: 4,8,12,24)")
    println("  --output-dir=PATH            Output directory (default: examples/plots/diurnal_shortwave)")
end

# Daily TOA insolation matches Kocken and Zeebe (2026),
# doi:10.1029/2025PA005287, src/insolation.f90.  The instantaneous form
# below is the corresponding solar-zenith-angle expression.
function paleo_insolation(
    eccentricity,
    obliquity_rad,
    lpx_rad;
    longitude_rad=pi / 2,
    latitude_rad=65 * pi / 180,
    solar_constant=1360.7,
    hour_angle=nothing,
)
    true_anomaly = longitude_rad - lpx_rad
    rho = (1 - eccentricity^2) / (1 + eccentricity * cos(true_anomaly))
    sin_delta = sin(obliquity_rad) * sin(longitude_rad)
    cos_delta = sqrt(max(0.0, 1 - sin_delta^2))
    sin_lat_sin_delta = sin(latitude_rad) * sin_delta
    cos_lat_cos_delta = cos(latitude_rad) * cos_delta

    if isnothing(hour_angle)
        cos_h0 = clamp(-sin_lat_sin_delta / cos_lat_cos_delta, -1.0, 1.0)
        sin_h0 = sqrt(max(0.0, 1 - cos_h0^2))
        h0 = acos(cos_h0)
        return solar_constant / (pi * rho^2) *
               (h0 * sin_lat_sin_delta + cos_lat_cos_delta * sin_h0)
    end

    return max(
        0.0,
        solar_constant / rho^2 *
        (sin_lat_sin_delta + cos_lat_cos_delta * cos(hour_angle)),
    )
end

function chion_instantaneous(shortwave_mean, latitude_deg, solar_longitude_deg, hour_angle)
    terms = Chion._diurnal_shortwave_integral_terms(latitude_deg, solar_longitude_deg)
    if terms.daylight_integral <= eps(Float64) ||
       terms.sunset_hour_angle <= 0 ||
       abs(hour_angle) > terms.sunset_hour_angle
        return 0.0
    end
    scale = shortwave_mean * 2pi / terms.daylight_integral
    return max(0.0, scale * (terms.sin_lat_sin_dec + terms.cos_lat_cos_dec * cos(hour_angle)))
end

function paleo_normalized(shortwave_mean, eccentricity, obliquity_rad, lpx_rad, latitude_deg, solar_longitude_deg, hour_angle; solar_constant=1360.7)
    daily_mean = paleo_insolation(
        eccentricity,
        obliquity_rad,
        lpx_rad;
        longitude_rad=deg2rad(solar_longitude_deg),
        latitude_rad=deg2rad(latitude_deg),
        solar_constant=solar_constant,
    )
    daily_mean <= 0 && return 0.0
    instantaneous = paleo_insolation(
        eccentricity,
        obliquity_rad,
        lpx_rad;
        longitude_rad=deg2rad(solar_longitude_deg),
        latitude_rad=deg2rad(latitude_deg),
        solar_constant=solar_constant,
        hour_angle=hour_angle,
    )
    return shortwave_mean * instantaneous / daily_mean
end

function substep_intervals(count)
    count >= 1 || error("Each substep count must be at least 1.")
    width = 2pi / count
    return [
        (-pi + (i - 1) * width, i == count ? pi : -pi + i * width)
        for i in 1:count
    ]
end

function substep_average_curve(shortwave_mean, latitude_deg, solar_longitude_deg, intervals, hour_angle)
    for (lo, hi) in intervals
        if (hour_angle >= lo && hour_angle <= hi) || (isapprox(hour_angle, pi) && isapprox(hi, pi))
            return Chion._diurnal_shortwave_interval_average(
                shortwave_mean,
                latitude_deg,
                solar_longitude_deg,
                lo,
                hi,
            )
        end
    end
    return 0.0
end

function write_csv(path, rows, substep_counts)
    open(path, "w") do io
        headers = ["substep_average_$(count)_w_m2" for count in substep_counts]
        println(io, join(["hour_angle_rad", "hour", "chion_w_m2", "paleo_normalized_w_m2", headers...], ','))
        for row in rows
            values = [row.hour_angle, row.hour, row.chion, row.paleo, row.substeps...]
            println(io, join((@sprintf("%.10f", value) for value in values), ','))
        end
    end
end

panel_label(index) = index <= 26 ? "($(Char('a' + index - 1)))" : "($(index))"

function write_pdf(path, cases, substep_counts, shortwave_mean; solar_longitude_deg)
    panels = Any[]
    substep_colors = [:purple, :green, :magenta, :orange]
    substep_linestyles = [:dash, :dot, :dashdot, :dashdotdot]
    for (index, case) in enumerate(cases)
        rows = case.rows
        hours = [row.hour for row in rows]
        ymax = 1.32 * maximum(max(row.chion, row.paleo, shortwave_mean, maximum(row.substeps)) for row in rows)
        panel = plot(
            title="$(case.latitude_deg)°N",
            xlabel="Hour from solar noon",
            ylabel="Shortwave Radiation (W m⁻²)",
            xlims=(-12, 12),
            xticks=-12:4:12,
            ylims=(0, ymax),
            grid=true,
            legend=index == 1 ? :topright : false,
            legendfontsize=7,
            left_margin=15Plots.mm,
            bottom_margin=6Plots.mm,
        )
        annotate!(panel, -14.4, 1.04 * ymax, text(panel_label(index), 12, :black, :left))
        hline!(panel, [shortwave_mean]; label=index == 1 ? "target daily mean" : false, color=:black, linestyle=:dash, linewidth=2)
        plot!(panel, hours, [row.paleo for row in rows]; label=index == 1 ? "paleo normalized" : false, color=:blue, linewidth=2.5)
        for substep_index in eachindex(substep_counts)
            plot!(panel, hours, [row.substeps[substep_index] for row in rows]; label=index == 1 ? "$(substep_counts[substep_index])-substep average" : false, color=substep_colors[mod1(substep_index, length(substep_colors))], linestyle=substep_linestyles[mod1(substep_index, length(substep_linestyles))], linewidth=1.8)
        end
        push!(panels, panel)
    end
    ncols = min(2, length(panels))
    nrows = cld(length(panels), ncols)
    figure = plot(panels...; layout=(nrows, ncols), size=(1150, 390 * nrows))
    savefig(figure, path)
end

function run_case(; latitude_deg, solar_longitude_deg, shortwave_mean, solar_constant, samples, substep_counts, output_dir)
    hour_angles = collect(range(-pi, pi; length=samples))
    intervals = [substep_intervals(count) for count in substep_counts]
    eccentricity = ZB18A_1_1_PRESENT_DAY.eccentricity
    obliquity_rad = ZB18A_1_1_PRESENT_DAY.obliquity_rad
    lpx_rad = ZB18A_1_1_PRESENT_DAY.lpx_rad

    rows = map(hour_angles) do h
        hour = 12 * h / pi
        (
            hour_angle=h,
            hour=hour,
            chion=chion_instantaneous(shortwave_mean, latitude_deg, solar_longitude_deg, h),
            paleo=paleo_normalized(shortwave_mean, eccentricity, obliquity_rad, lpx_rad, latitude_deg, solar_longitude_deg, h; solar_constant),
            substeps=[substep_average_curve(shortwave_mean, latitude_deg, solar_longitude_deg, interval, h) for interval in intervals],
        )
    end

    mean_chion = sum(row.chion for row in rows) / length(rows)
    mean_paleo = sum(row.paleo for row in rows) / length(rows)
    max_abs_diff = maximum(abs(row.chion - row.paleo) for row in rows)
    paleo_daily = paleo_insolation(
        eccentricity,
        obliquity_rad,
        lpx_rad;
        longitude_rad=deg2rad(solar_longitude_deg),
        latitude_rad=deg2rad(latitude_deg),
        solar_constant=solar_constant,
    )

    slug = @sprintf("lat_%+.2f_lon_%06.2f", latitude_deg, solar_longitude_deg)
    slug = replace(slug, "+" => "p", "-" => "m", "." => "p")
    csv_path = joinpath(output_dir, "$(slug).csv")
    write_csv(csv_path, rows, substep_counts)

    return (
        latitude_deg=latitude_deg,
        solar_longitude_deg=solar_longitude_deg,
        rows=rows,
        csv_path=csv_path,
        substeps=join(substep_counts, ";"),
        mean_chion=mean_chion,
        mean_paleo=mean_paleo,
        max_abs_diff=max_abs_diff,
        paleo_daily=paleo_daily,
    )
end

function main(args)
    if has_flag(args, "help")
        print_help()
        return
    end

    output_dir = arg_value(args, "output-dir", DEFAULT_OUTPUT_DIR)
    mkpath(output_dir)

    latitudes = parse_list(arg_value(args, "latitudes", "40,50,60,70"), Float64)
    solar_longitudes = parse_list(arg_value(args, "solar-longitudes", "90"), Float64)
    shortwave_mean = parse(Float64, arg_value(args, "shortwave-mean", "200"))
    solar_constant = parse(Float64, arg_value(args, "solar-constant", "1360.7"))
    samples = parse(Int, arg_value(args, "samples", "289"))
    substep_counts = parse_list(arg_value(args, "substeps", "4,8,12,24"), Int)
    isempty(substep_counts) && error("`--substeps` must contain at least one count.")

    summaries = NamedTuple[]
    for latitude in latitudes, solar_longitude in solar_longitudes
        push!(
            summaries,
            run_case(
                latitude_deg=latitude,
                solar_longitude_deg=solar_longitude,
                shortwave_mean=shortwave_mean,
                solar_constant=solar_constant,
                samples=samples,
                substep_counts=substep_counts,
                output_dir=output_dir,
            ),
        )
    end

    pdf_paths = Dict{Float64,String}()
    for solar_longitude in solar_longitudes
        solar_cases = filter(case -> case.solar_longitude_deg == solar_longitude, summaries)
        longitude_slug = replace(@sprintf("lon_%06.2f", solar_longitude), "+" => "p", "-" => "m", "." => "p")
        pdf_path = joinpath(output_dir, "diurnal_shortwave_$(longitude_slug).pdf")
        write_pdf(pdf_path, solar_cases, substep_counts, shortwave_mean; solar_longitude_deg=solar_longitude)
        pdf_paths[solar_longitude] = pdf_path
        println("Wrote plot: ", pdf_path)
    end

    summary_path = joinpath(output_dir, "summary.csv")
    open(summary_path, "w") do io
        println(io, "latitude_deg,solar_longitude_deg,substeps,mean_chion,mean_paleo,max_abs_diff,paleo_daily_w_m2,csv_path,pdf_path")
        for row in summaries
            @printf(
                io,
                "%.6f,%.6f,%s,%.10f,%.10f,%.10f,%.10f,%s,%s\n",
                row.latitude_deg,
                row.solar_longitude_deg,
                row.substeps,
                row.mean_chion,
                row.mean_paleo,
                row.max_abs_diff,
                row.paleo_daily,
                row.csv_path,
                pdf_paths[row.solar_longitude_deg],
            )
        end
    end

    println("Wrote summary: ", summary_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
