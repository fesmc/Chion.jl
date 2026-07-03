#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Printf
using Chion

const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "..", "plots", "diurnal_shortwave")

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
    println("  --latitudes=LIST             Comma-separated degrees north (default: 65)")
    println("  --solar-longitudes=LIST      Comma-separated degrees, e.g. 0,90,180,270 (default: 90)")
    println("  --shortwave-mean=VALUE       Target daily mean W m^-2 for normalized curves (default: 200)")
    println("  --eccentricity=VALUE         Paleo orbital eccentricity (default: 0.0)")
    println("  --obliquity-deg=VALUE        Paleo obliquity degrees (default: 23.439291)")
    println("  --lpx-deg=VALUE              Paleo longitude of perihelion degrees (default: 0)")
    println("  --solar-constant=VALUE       Paleo solar constant W m^-2 (default: 1360.7)")
    println("  --samples=N                  Samples across one day (default: 289)")
    println("  --threshold=VALUE            Chion substep threshold W m^-2 (default: 0)")
    println("  --max-substeps=N             Chion max adaptive substeps, 1 to 24 (default: 3)")
    println("  --air-temperature-c=VALUE    Daily mean air temperature for substep criterion (default: 0)")
    println("  --output-dir=PATH            Output directory (default: examples/plots/diurnal_shortwave)")
end

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

function paleo_normalized(shortwave_mean, eccentricity, obliquity_rad, lpx_rad, latitude_deg, solar_longitude_deg, hour_angle)
    daily_mean = paleo_insolation(
        eccentricity,
        obliquity_rad,
        lpx_rad;
        longitude_rad=deg2rad(solar_longitude_deg),
        latitude_rad=deg2rad(latitude_deg),
    )
    daily_mean <= 0 && return 0.0
    instantaneous = paleo_insolation(
        eccentricity,
        obliquity_rad,
        lpx_rad;
        longitude_rad=deg2rad(solar_longitude_deg),
        latitude_rad=deg2rad(latitude_deg),
        hour_angle=hour_angle,
    )
    return shortwave_mean * instantaneous / daily_mean
end

function substep_intervals(latitude_deg, solar_longitude_deg, shortwave_mean, threshold, max_substeps, air_temperature_c)
    count = Chion._diurnal_shortwave_substep_count(
        1.0,
        shortwave_mean,
        air_temperature_c + 273.15,
        -8.0 + 273.15,
        latitude_deg,
        solar_longitude_deg,
        threshold,
        max_substeps,
    )
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

function write_csv(path, rows)
    open(path, "w") do io
        println(io, "hour_angle_rad,hour,chion_w_m2,paleo_normalized_w_m2,substep_average_w_m2")
        for row in rows
            @printf(
                io,
                "%.10f,%.10f,%.10f,%.10f,%.10f\n",
                row.hour_angle,
                row.hour,
                row.chion,
                row.paleo,
                row.substep,
            )
        end
    end
end

function svg_polyline(points, width, height, xmin, xmax, ymin, ymax)
    coords = String[]
    for (x, y) in points
        px = 70 + (x - xmin) / (xmax - xmin) * (width - 100)
        py = height - 55 - (y - ymin) / (ymax - ymin) * (height - 95)
        push!(coords, @sprintf("%.2f,%.2f", px, py))
    end
    return join(coords, " ")
end

function write_svg(path, rows; title, subtitle)
    width = 980
    height = 540
    xmin, xmax = -12.0, 12.0
    ymax = maximum(max(row.chion, row.paleo, row.substep) for row in rows)
    ymax = max(1.0, ceil(ymax / 50) * 50)
    ymin = 0.0

    chion_points = [(row.hour, row.chion) for row in rows]
    paleo_points = [(row.hour, row.paleo) for row in rows]
    substep_points = [(row.hour, row.substep) for row in rows]

    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""")
        println(io, """<rect width="100%" height="100%" fill="white"/>""")
        println(io, """<text x="70" y="36" font-family="Arial" font-size="20" font-weight="700">$title</text>""")
        println(io, """<text x="70" y="60" font-family="Arial" font-size="13" fill="#555">$subtitle</text>""")
        println(io, """<line x1="70" y1="$(height - 55)" x2="$(width - 30)" y2="$(height - 55)" stroke="#333"/>""")
        println(io, """<line x1="70" y1="80" x2="70" y2="$(height - 55)" stroke="#333"/>""")

        for hour in -12:4:12
            x = 70 + (hour - xmin) / (xmax - xmin) * (width - 100)
            println(io, """<line x1="$x" y1="80" x2="$x" y2="$(height - 55)" stroke="#eee"/>""")
            println(io, """<text x="$x" y="$(height - 32)" text-anchor="middle" font-family="Arial" font-size="12">$hour</text>""")
        end
        for y in range(ymin, ymax; length=6)
            py = height - 55 - (y - ymin) / (ymax - ymin) * (height - 95)
            println(io, """<line x1="70" y1="$py" x2="$(width - 30)" y2="$py" stroke="#eee"/>""")
            println(io, """<text x="58" y="$(py + 4)" text-anchor="end" font-family="Arial" font-size="12">$(round(Int, y))</text>""")
        end

        println(io, """<polyline points="$(svg_polyline(substep_points, width, height, xmin, xmax, ymin, ymax))" fill="none" stroke="#777" stroke-width="2.5" stroke-dasharray="6 5"/>""")
        println(io, """<polyline points="$(svg_polyline(paleo_points, width, height, xmin, xmax, ymin, ymax))" fill="none" stroke="#2b6cb0" stroke-width="3"/>""")
        println(io, """<polyline points="$(svg_polyline(chion_points, width, height, xmin, xmax, ymin, ymax))" fill="none" stroke="#c2410c" stroke-width="2" stroke-dasharray="3 4"/>""")
        println(io, """<text x="$(width - 250)" y="95" font-family="Arial" font-size="13" fill="#2b6cb0">paleo normalized</text>""")
        println(io, """<text x="$(width - 250)" y="117" font-family="Arial" font-size="13" fill="#c2410c">Chion instantaneous</text>""")
        println(io, """<text x="$(width - 250)" y="139" font-family="Arial" font-size="13" fill="#777">Chion substep averages</text>""")
        println(io, """<text x="$(width / 2)" y="$(height - 8)" text-anchor="middle" font-family="Arial" font-size="13">Hour from solar noon</text>""")
        println(io, """<text x="18" y="$(height / 2)" transform="rotate(-90 18 $(height / 2))" text-anchor="middle" font-family="Arial" font-size="13">Shortwave (W m^-2)</text>""")
        println(io, "</svg>")
    end
end

function run_case(; latitude_deg, solar_longitude_deg, shortwave_mean, eccentricity, obliquity_deg, lpx_deg, solar_constant, samples, threshold, max_substeps, air_temperature_c, output_dir)
    hour_angles = collect(range(-pi, pi; length=samples))
    intervals = substep_intervals(latitude_deg, solar_longitude_deg, shortwave_mean, threshold, max_substeps, air_temperature_c)
    obliquity_rad = deg2rad(obliquity_deg)
    lpx_rad = deg2rad(lpx_deg)

    rows = map(hour_angles) do h
        hour = 12 * h / pi
        (
            hour_angle=h,
            hour=hour,
            chion=chion_instantaneous(shortwave_mean, latitude_deg, solar_longitude_deg, h),
            paleo=paleo_normalized(shortwave_mean, eccentricity, obliquity_rad, lpx_rad, latitude_deg, solar_longitude_deg, h),
            substep=substep_average_curve(shortwave_mean, latitude_deg, solar_longitude_deg, intervals, h),
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
    svg_path = joinpath(output_dir, "$(slug).svg")
    write_csv(csv_path, rows)
    write_svg(
        svg_path,
        rows;
        title=@sprintf("Diurnal shortwave: %.2f deg lat, %.2f deg solar longitude", latitude_deg, solar_longitude_deg),
        subtitle=@sprintf(
            "target daily mean %.1f W m^-2; paleo e=%.5f, obliquity=%.4f deg, lpx=%.2f deg; paleo daily TOA %.1f W m^-2; max diff %.3g W m^-2",
            shortwave_mean,
            eccentricity,
            obliquity_deg,
            lpx_deg,
            paleo_daily,
            max_abs_diff,
        ),
    )

    return (
        latitude_deg=latitude_deg,
        solar_longitude_deg=solar_longitude_deg,
        csv_path=csv_path,
        svg_path=svg_path,
        substeps=length(intervals),
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

    latitudes = parse_list(arg_value(args, "latitudes", "65"), Float64)
    solar_longitudes = parse_list(arg_value(args, "solar-longitudes", "90"), Float64)
    shortwave_mean = parse(Float64, arg_value(args, "shortwave-mean", "200"))
    eccentricity = parse(Float64, arg_value(args, "eccentricity", "0.0"))
    obliquity_deg = parse(Float64, arg_value(args, "obliquity-deg", "23.439291"))
    lpx_deg = parse(Float64, arg_value(args, "lpx-deg", "0.0"))
    solar_constant = parse(Float64, arg_value(args, "solar-constant", "1360.7"))
    samples = parse(Int, arg_value(args, "samples", "289"))
    threshold = parse(Float64, arg_value(args, "threshold", "0.0"))
    max_substeps = parse(Int, arg_value(args, "max-substeps", "3"))
    air_temperature_c = parse(Float64, arg_value(args, "air-temperature-c", "0.0"))

    summaries = NamedTuple[]
    for latitude in latitudes, solar_longitude in solar_longitudes
        push!(
            summaries,
            run_case(
                latitude_deg=latitude,
                solar_longitude_deg=solar_longitude,
                shortwave_mean=shortwave_mean,
                eccentricity=eccentricity,
                obliquity_deg=obliquity_deg,
                lpx_deg=lpx_deg,
                solar_constant=solar_constant,
                samples=samples,
                threshold=threshold,
                max_substeps=max_substeps,
                air_temperature_c=air_temperature_c,
                output_dir=output_dir,
            ),
        )
    end

    summary_path = joinpath(output_dir, "summary.csv")
    open(summary_path, "w") do io
        println(io, "latitude_deg,solar_longitude_deg,substeps,mean_chion,mean_paleo,max_abs_diff,paleo_daily_w_m2,csv_path,svg_path")
        for row in summaries
            @printf(
                io,
                "%.6f,%.6f,%d,%.10f,%.10f,%.10f,%.10f,%s,%s\n",
                row.latitude_deg,
                row.solar_longitude_deg,
                row.substeps,
                row.mean_chion,
                row.mean_paleo,
                row.max_abs_diff,
                row.paleo_daily,
                row.csv_path,
                row.svg_path,
            )
        end
    end

    println("Wrote summary: ", summary_path)
    for row in summaries
        println("Wrote plot: ", row.svg_path)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
