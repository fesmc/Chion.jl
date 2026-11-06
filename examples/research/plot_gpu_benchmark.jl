#!/usr/bin/env julia

using Plots
using Statistics

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const BENCHMARK_DIR = joinpath(ROOT, "data", "benchmark")
const OUTPUT_DIR = joinpath(ROOT, "examples", "plots")

const DATASETS = [
    (label="GPU A100", backend="gpu", file="snowpack-columns-gpu_a100.csv"),
    (label="GPU H100", backend="gpu", file="snowpack-columns_gpu_h100.csv"),
    (label="CPU (8 threads)", backend="threads", file="snowpack-columns-cpu8.csv"),
    (label="CPU (16 threads)", backend="threads", file="snowpack-columns-cpu16.csv"),
    (label="CPU (32 threads)", backend="threads", file="snowpack-columns-cpu32.csv"),
    (label="CPU (64 threads)", backend="threads", file="snowpack-columns-cpu64.csv"),
    (label="CPU (128 threads)", backend="threads", file="snowpack-columns-cpu128.csv"),
]

"""Read the first complete set of three repetitions for each configuration."""
function read_timings(path, expected_backend)
    lines = readlines(path)
    header = split(first(lines), ',')
    @assert header == ["backend", "threads", "columns", "ntot", "repetition", "seconds", "status"]

    # The A100 file includes appended duplicate batches.  Retain the first value
    # for repetitions 1:3 so every configuration contributes exactly three runs.
    timings = Dict{Tuple{Int, Int}, Dict{Int, Float64}}()
    for line in Iterators.drop(lines, 1)
        backend, _, columns, ntot, repetition, seconds, status = split(line, ',')
        backend == expected_backend && status == "complete" || continue
        key = (parse(Int, columns), parse(Int, ntot))
        runs = get!(timings, key, Dict{Int, Float64}())
        run = parse(Int, repetition)
        run in 1:3 && !haskey(runs, run) && (runs[run] = parse(Float64, seconds))
    end

    incomplete = [key for (key, runs) in timings if length(runs) != 3]
    isempty(incomplete) || error("Expected three complete runs in $path; incomplete: $incomplete")
    return Dict(key => mean(values(runs)) for (key, runs) in timings)
end

function shared_configurations(timings)
    columns = reduce(intersect, [Set(first(key) for key in keys(data)) for data in timings])
    ntots = reduce(intersect, [Set(last(key) for key in keys(data)) for data in timings])
    return sort!(collect(columns)), sort!(collect(ntots))
end

function ylimits(values)
    return (minimum(values) / 1.2, maximum(values) * 1.25)
end

function panel_title(index, label)
    return "($(Char('a' + index - 1)))  $label"
end

seconds_label(seconds) = "$(round(Int, seconds)) s"

"""Place labels above nearby values with enough vertical separation to remain legible."""
function label_positions(values)
    positions = similar(values, Float64)
    previous = 0.0
    for index in sortperm(values)
        positions[index] = max(values[index] * 1.12, previous * 1.18)
        previous = positions[index]
    end
    return positions
end

function common_panel(; title, columns, ntots, data, limits)
    return plot(
        xlabel="Number of columns",
        ylabel="Time (s)",
        title=title,
        titlelocation=(-0.06, 1.0),
        titlefonthalign=:left,
        titlefontvalign=:bottom,
        xscale=:log2,
        xlims=(columns[1] / 1.2, columns[end] * 1.3),
        xticks=(columns, string.(columns)),
        xrotation=25,
        yscale=:log10,
        ylims=limits,
        legend=:topleft,
        grid=true,
        guidefont=font(14, "helvetica bold"),
        tickfont=font(11, "helvetica bold"),
        titlefont=font(16, "helvetica bold"),
        legendfont=font(11, "helvetica bold"),
    )
end

function save_figure(figure, basename)
    mkpath(OUTPUT_DIR)
    for extension in ("png", "pdf")
        path = joinpath(OUTPUT_DIR, "$basename.$extension")
        savefig(figure, path)
        println("Wrote $path")
    end
end

function plot_benchmarks()
    timings = [read_timings(joinpath(BENCHMARK_DIR, dataset.file), dataset.backend)
               for dataset in DATASETS]
    columns, ntots = shared_configurations(timings)
    all_values = [data[(ncol, ntot)] for data in timings for ntot in ntots for ncol in columns]
    limits = ylimits(all_values)

    panels = Any[]
    for (index, (dataset, data)) in enumerate(zip(DATASETS, timings))
        panel = common_panel(
            title=panel_title(index, dataset.label), columns=columns, ntots=ntots,
            data=data, limits=limits,
        )
        for ntot in ntots
            values = [data[(ncol, ntot)] for ncol in columns]
            plot!(panel, columns, values;
                marker=:circle, markerstrokecolor=:auto, linewidth=3, label="Ntot = $ntot")
            annotate!(panel, columns[end] * 1.065, values[end],
                text(seconds_label(values[end]), 9, "helvetica bold", :left))
        end
        push!(panels, panel)
    end
    figure = plot(panels...; layout=(4, 2), size=(1400, 1900),
        left_margin=10Plots.mm, bottom_margin=12Plots.mm)
    save_figure(figure, "gpu_cpu_benchmark_a100_h100")

    comparison_ntot = 20
    comparison_series = [[data[(ncol, comparison_ntot)] for ncol in columns] for data in timings]
    comparison_values = reduce(vcat, comparison_series)
    annotation_values = reduce(vcat, [label_positions([series[index] for series in comparison_series])
                                      for index in eachindex(columns)])
    comparison = plot(
        xlabel="Number of columns",
        ylabel="Time (s)",
        title="(a)",
        titlelocation=(-0.09, 1.0),
        titlefonthalign=:left,
        titlefontvalign=:bottom,
        xscale=:log2,
        yscale=:log10,
        xlims=(columns[1] / 1.2, columns[end] * 1.7),
        xticks=(columns, string.(columns)),
        xrotation=25,
        ylims=(minimum(comparison_values) / 1.2,
            max(maximum(comparison_values) * 1.25, maximum(annotation_values) * 1.08)),
        legend=:topleft,
        grid=true,
        size=(850, 600),
        guidefont=font(14, "helvetica bold"),
        tickfont=font(11, "helvetica bold"),
        titlefont=font(16, "helvetica bold"),
        legendfont=font(9, "helvetica bold"),
        left_margin=10Plots.mm,
        top_margin=10Plots.mm,
        bottom_margin=12Plots.mm,
    )
    markers = (:circle, :diamond, :utriangle, :rtriangle, :square, :star5, :hexagon)
    for (dataset, values, marker) in zip(DATASETS, comparison_series, markers)
        plot!(comparison, columns, values;
            marker, markerstrokecolor=:auto, linewidth=3, label=dataset.label)
    end
    for (index, ncol) in enumerate(columns)
        positions = label_positions([series[index] for series in comparison_series])
        for (seconds, label_y) in zip((series[index] for series in comparison_series), positions)
            xoffset = index == length(columns) ? 1.07 : 1.04
            annotate!(comparison, ncol * xoffset, label_y,
                text(seconds_label(seconds), 8, "helvetica bold", :left))
        end
    end
    max_columns = columns[end]
    ntot_series = [[data[(max_columns, ntot)] for ntot in ntots] for data in timings]
    ntot_values = reduce(vcat, ntot_series)
    ntot_annotation_values = reduce(vcat, [label_positions([series[index] for series in ntot_series])
                                           for index in eachindex(ntots)])
    ntot_comparison = plot(
        xlabel="Number of maximum layers",
        ylabel="Time (s)",
        title="(b)",
        titlelocation=(-0.09, 1.0),
        titlefonthalign=:left,
        titlefontvalign=:bottom,
        xticks=ntots,
        xlims=(first(ntots) - 0.5, last(ntots) + 1.5),
        yscale=:log10,
        ylims=(minimum(ntot_values) / 1.2,
            max(maximum(ntot_values) * 1.25, maximum(ntot_annotation_values) * 1.08)),
        legend=false,
        grid=true,
        size=(850, 600),
        guidefont=font(14, "helvetica bold"),
        tickfont=font(11, "helvetica bold"),
        titlefont=font(16, "helvetica bold"),
        left_margin=10Plots.mm,
        top_margin=10Plots.mm,
        bottom_margin=12Plots.mm,
    )
    for (dataset, values, marker) in zip(DATASETS, ntot_series, markers)
        plot!(ntot_comparison, ntots, values;
            marker, markerstrokecolor=:auto, linewidth=3, label=dataset.label)
    end
    for (index, ntot) in enumerate(ntots)
        positions = label_positions([series[index] for series in ntot_series])
        for (seconds, label_y) in zip((series[index] for series in ntot_series), positions)
            annotate!(ntot_comparison, ntot + 0.35, label_y,
                text(seconds_label(seconds), 8, "helvetica bold", :left))
        end
    end
    comparison_figure = plot(comparison, ntot_comparison; layout=(1, 2), size=(1500, 650),
        left_margin=10Plots.mm, bottom_margin=12Plots.mm)
    save_figure(comparison_figure, "benchmark_gpu_cpu_comparisons")
end

plot_benchmarks()
