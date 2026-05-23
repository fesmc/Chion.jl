#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using NCDatasets
using Plots
using Printf
using Statistics

const DEFAULT_MAR_FILE =
    "/p/projects/ou/labs/ai/Nils/MAR3.14/MARv3.14.3-10km-daily-ERA5-1940-1980_daily_climatology.nc"
const DEFAULT_CHION_FILE =
    joinpath(@__DIR__, "..", "plots", "gris_forcing_file_simulation", "chion_run_final_state.nc")
const DEFAULT_OUTPUT_DIR =
    joinpath(@__DIR__, "..", "plots", "gris_mar_monthly_comparison")

struct ComparisonSpec
    key::Symbol
    title::String
    chion_var::String
    mar_var::String
    mar_dim::Symbol
    temporal::Symbol
    units::String
end

arg_value(args, name, default="") = begin
    prefix = "--$(name)="
    for arg in args
        startswith(arg, prefix) && return arg[length(prefix)+1:end]
    end
    return default
end

has_flag(args, name) = any(==("--$(name)"), args)

function print_help()
    println("Usage:")
    println("  julia --project=. examples/scripts/compare_gris_mar_monthly.jl [options]")
    println()
    println("Options:")
    println("  --mar-file=PATH          MAR daily climatology NetCDF")
    println("  --chion-file=PATH        Chion monthly NetCDF")
    println("  --output-dir=PATH        Directory for PNG/CSV outputs")
    println("  --chion-period=last      last|mean|all|N, where N is a 1-based repeated year")
    println("  --mar-sector=N           1-based MAR SECTOR index for SMB/RU/albedo (default: 1)")
    println("  --mar-runoff-var=NAME    RU or RU2 (default: RU)")
    println("  --mar-sublimation-var=NAME  MAR sublimation variable (default: SU)")
    println("  --mar-albedo-var=NAME    AL1 or AL2 (default: AL2)")
    println("  --mask-threshold=VALUE   MAR/Chion domain mask threshold (default: 50)")
    println("  --no-maps                Skip spatial map and scatter plots")
end

function clean_array(A)
    out = Array{Float64}(undef, size(A))
    @inbounds for I in eachindex(A)
        value = A[I]
        if ismissing(value)
            out[I] = NaN
        else
            x = Float64(value)
            out[I] = isfinite(x) && abs(x) < 1.0e18 ? x : NaN
        end
    end
    return out
end

read_clean(ds, name, inds...) = clean_array(Array(ds[name][inds...]))

function nc_attr(ds, name::AbstractString, default)
    return haskey(ds.attrib, name) ? ds.attrib[name] : default
end

function records_written(ds, varname::AbstractString)
    ntime = size(ds[varname], 1)
    raw = nc_attr(ds, "records_written", "")
    isempty(String(raw)) && return ntime
    parsed = tryparse(Int, String(raw))
    return isnothing(parsed) ? ntime : min(parsed, ntime)
end

function grid_shape(ds)
    haskey(ds, "domain_mask") && return size(ds["domain_mask"])
    haskey(ds, "x") && haskey(ds, "y") && return (length(ds["x"][:]), length(ds["y"][:]))
    for name in keys(ds)
        ndims(ds[name]) >= 3 && return (size(ds[name], 2), size(ds[name], 3))
    end
    error("Could not infer Chion grid shape from NetCDF variables.")
end

month_sequence(nrecord::Int) = [mod(i - 1, 12) + 1 for i in 1:nrecord]

function finite_mask(A)
    mask = falses(size(A))
    @inbounds for I in eachindex(A)
        mask[I] = isfinite(A[I])
    end
    return mask
end

function read_mar_daily(ds, spec::ComparisonSpec, sector::Int)
    if spec.mar_dim === :sector
        return read_clean(ds, spec.mar_var, :, :, sector, :)
    elseif spec.mar_dim === :single_sector
        return read_clean(ds, spec.mar_var, :, :, 1, :)
    elseif spec.mar_dim === :none
        return read_clean(ds, spec.mar_var, :, :, :)
    else
        error("Unsupported MAR dimension selector $(spec.mar_dim)")
    end
end

function aggregate_mar_monthly(ds, spec::ComparisonSpec, sector::Int)
    months = Int.(vec(clean_array(Array(ds["MM"][:]))))
    daily = read_mar_daily(ds, spec, sector)
    nx, ny, ntime = size(daily)
    ntime == length(months) || error("MAR $(spec.mar_var) time length $(ntime) does not match MM length $(length(months)).")

    monthly = fill(NaN, 12, nx, ny)
    counts = zeros(Int, nx, ny)
    accum = zeros(Float64, nx, ny)
    for month in 1:12
        fill!(accum, 0.0)
        fill!(counts, 0)
        for t in 1:ntime
            months[t] == month || continue
            @inbounds for j in 1:ny, i in 1:nx
                value = daily[i, j, t]
                if isfinite(value)
                    accum[i, j] += value
                    counts[i, j] += 1
                end
            end
        end
        @inbounds for j in 1:ny, i in 1:nx
            if counts[i, j] > 0
                monthly[month, i, j] = spec.temporal === :mean ? accum[i, j] / counts[i, j] : accum[i, j]
            end
        end
    end
    return monthly
end

function selected_chion_indices(months::Vector{Int}, period::AbstractString)
    nmonth = length(months)
    if period == "all"
        return collect(1:nmonth)
    elseif period == "last"
        nmonth >= 12 || error("Chion file has fewer than 12 monthly records.")
        return collect((nmonth - 11):nmonth)
    elseif all(isdigit, period)
        year = parse(Int, period)
        first_idx = (year - 1) * 12 + 1
        last_idx = year * 12
        1 <= first_idx <= last_idx <= nmonth || error("--chion-period=$(year) is outside available records 1:$(div(nmonth, 12)).")
        return collect(first_idx:last_idx)
    else
        error("--chion-period must be last, mean, all, or a 1-based year index.")
    end
end

function read_chion_monthly(ds, name::AbstractString, period::AbstractString)
    nrecord = records_written(ds, name)
    months = haskey(ds, "month_of_year") ?
        Int.(vec(clean_array(Array(ds["month_of_year"][1:nrecord])))) :
        month_sequence(nrecord)
    nmonth = length(months)
    nx, ny = grid_shape(ds)

    if period == "mean"
        accum = zeros(Float64, 12, nx, ny)
        counts = zeros(Int, 12, nx, ny)
        for t in 1:nmonth
            m = months[t]
            field = read_clean(ds, name, t, :, :)
            @inbounds for j in 1:ny, i in 1:nx
                value = field[i, j]
                if isfinite(value)
                    accum[m, i, j] += value
                    counts[m, i, j] += 1
                end
            end
        end
        out = fill(NaN, 12, nx, ny)
        @inbounds for I in eachindex(out)
            counts[I] > 0 && (out[I] = accum[I] / counts[I])
        end
        return out, collect(1:12)
    end

    idxs = selected_chion_indices(months, period)
    return read_clean(ds, name, idxs, :, :), months[idxs]
end

function max_abs_finite_difference(a, b)
    maxdiff = 0.0
    found = false
    @inbounds for I in eachindex(a, b)
        av = a[I]
        bv = b[I]
        if isfinite(av) && isfinite(bv)
            maxdiff = max(maxdiff, abs(av - bv))
            found = true
        end
    end
    return found ? maxdiff : NaN
end

function validate_chion_smb_variables(chion_ds, period::AbstractString)
    haskey(chion_ds, "monthly_smb") || return nothing
    haskey(chion_ds, "monthly_net_ice_sheet_forcing") || return nothing

    climatic_smb, _ = read_chion_monthly(chion_ds, "monthly_smb", period)
    net_ice_forcing, _ = read_chion_monthly(chion_ds, "monthly_net_ice_sheet_forcing", period)
    maxdiff = max_abs_finite_difference(climatic_smb, net_ice_forcing)
    if isfinite(maxdiff) && maxdiff <= 1.0e-6
        error(
            "The selected Chion file has identical `monthly_smb` and " *
            "`monthly_net_ice_sheet_forcing` over --chion-period=$(period). " *
            "This usually means the NetCDF was produced before `monthly_smb` " *
            "was changed to climatic SMB. Regenerate the Chion monthly output " *
            "with the current code before comparing both diagnostics.",
        )
    end
    return nothing
end

function repeat_mar_to_months(mar_monthly, months::Vector{Int})
    out = Array{Float64}(undef, length(months), size(mar_monthly, 2), size(mar_monthly, 3))
    for (k, month) in enumerate(months)
        out[k, :, :] .= mar_monthly[month, :, :]
    end
    return out
end

function weight_matrix(mar_ds, chion_ds, threshold::Float64)
    nx, ny = grid_shape(chion_ds)
    area = haskey(mar_ds, "AREA") ? read_clean(mar_ds, "AREA", :, :) : ones(nx, ny)
    mar_mask = haskey(mar_ds, "MSK") ? read_clean(mar_ds, "MSK", :, :) : fill(threshold, size(area))
    chion_mask = haskey(chion_ds, "domain_mask") ? read_clean(chion_ds, "domain_mask", :, :) : fill(threshold, size(area))
    weights = zeros(Float64, size(area))
    @inbounds for I in eachindex(area)
        if isfinite(area[I]) && area[I] > 0 &&
           isfinite(mar_mask[I]) && mar_mask[I] >= threshold &&
           isfinite(chion_mask[I]) && chion_mask[I] >= threshold
            weights[I] = area[I]
        end
    end
    return weights
end

function apply_mask!(A, weights)
    if ndims(A) == 3
        @inbounds for t in axes(A, 1), j in axes(A, 3), i in axes(A, 2)
            weights[i, j] > 0 || (A[t, i, j] = NaN)
        end
    elseif ndims(A) == 2
        @inbounds for j in axes(A, 2), i in axes(A, 1)
            weights[i, j] > 0 || (A[i, j] = NaN)
        end
    else
        error("Unsupported array rank $(ndims(A)) for spatial mask.")
    end
    return A
end

function weighted_mean(field, weights)
    total = 0.0
    wsum = 0.0
    @inbounds for I in eachindex(field, weights)
        value = field[I]
        w = weights[I]
        if isfinite(value) && w > 0
            total += value * w
            wsum += w
        end
    end
    return wsum > 0 ? total / wsum : NaN
end

function weighted_series(A, weights)
    out = Vector{Float64}(undef, size(A, 1))
    for t in axes(A, 1)
        out[t] = weighted_mean(view(A, t, :, :), weights)
    end
    return out
end

function reduce_period(A, temporal::Symbol)
    nx, ny = size(A, 2), size(A, 3)
    out = fill(NaN, nx, ny)
    counts = zeros(Int, nx, ny)
    accum = zeros(Float64, nx, ny)
    for t in axes(A, 1)
        @inbounds for j in 1:ny, i in 1:nx
            value = A[t, i, j]
            if isfinite(value)
                accum[i, j] += value
                counts[i, j] += 1
            end
        end
    end
    @inbounds for j in 1:ny, i in 1:nx
        if counts[i, j] > 0
            out[i, j] = temporal === :mean ? accum[i, j] / counts[i, j] : accum[i, j]
        end
    end
    return out
end

function paired_stats(chion, mar, weights)
    pairs_chion = Float64[]
    pairs_mar = Float64[]
    sq = 0.0
    abs_sum = 0.0
    bias_num = 0.0
    wsum = 0.0
    @inbounds for I in eachindex(chion, mar, weights)
        c = chion[I]
        m = mar[I]
        w = weights[I]
        if isfinite(c) && isfinite(m) && w > 0
            d = c - m
            bias_num += d * w
            sq += d * d * w
            abs_sum += abs(d) * w
            wsum += w
            push!(pairs_chion, c)
            push!(pairs_mar, m)
        end
    end
    corr = length(pairs_chion) > 1 && std(pairs_chion) > 0 && std(pairs_mar) > 0 ? cor(pairs_chion, pairs_mar) : NaN
    return (
        n=length(pairs_chion),
        chion_mean=weighted_mean(chion, weights),
        mar_mean=weighted_mean(mar, weights),
        bias=wsum > 0 ? bias_num / wsum : NaN,
        mae=wsum > 0 ? abs_sum / wsum : NaN,
        rmse=wsum > 0 ? sqrt(sq / wsum) : NaN,
        corr=corr,
    )
end

function write_monthly_csv(path, rows)
    open(path, "w") do io
        println(io, "variable,month_index,month_of_year,chion_mean,mar_mean,bias,units")
        for row in rows
            @printf(io, "%s,%d,%d,%.10g,%.10g,%.10g,%s\n",
                row.variable, row.month_index, row.month_of_year,
                row.chion_mean, row.mar_mean, row.bias, row.units)
        end
    end
end

function write_stats_csv(path, rows)
    open(path, "w") do io
        println(io, "variable,n,chion_mean,mar_mean,bias,mae,rmse,corr,units")
        for row in rows
            s = row.stats
            @printf(io, "%s,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%s\n",
                row.variable, s.n, s.chion_mean, s.mar_mean, s.bias, s.mae, s.rmse, s.corr, row.units)
        end
    end
end

function save_timeseries_plot(path, spec, months, chion_series, mar_series)
    p = plot(1:length(months), chion_series;
        label="Chion", lw=2.5, marker=:circle, xlabel="Selected Chion month",
        ylabel=spec.units, title="$(spec.title) monthly domain mean", framestyle=:box)
    plot!(p, 1:length(months), mar_series; label="MAR", lw=2.5, marker=:diamond)
    plot!(p, 1:length(months), chion_series .- mar_series; label="Chion - MAR", lw=2, ls=:dash)
    xticks!(p, 1:length(months), string.(months))
    savefig(p, path)
end

function finite_extrema(arrays...)
    lo = Inf
    hi = -Inf
    for A in arrays
        @inbounds for value in A
            if isfinite(value)
                lo = min(lo, value)
                hi = max(hi, value)
            end
        end
    end
    if !isfinite(lo) || !isfinite(hi)
        return (NaN, NaN)
    elseif lo == hi
        pad = max(abs(lo), 1.0) * 0.05
        return (lo - pad, hi + pad)
    end
    return (lo, hi)
end

function save_map_plot(path, spec, x, y, chion_map, mar_map)
    diff = chion_map .- mar_map
    common_clim = finite_extrema(chion_map, mar_map)
    diff_abs = maximum(abs, filter(isfinite, vec(diff)); init=0.0)
    diff_clim = diff_abs > 0 ? (-diff_abs, diff_abs) : (-1.0, 1.0)
    p1 = heatmap(x, y, chion_map'; title="Chion", xlabel="x", ylabel="y", colorbar_title=spec.units, aspect_ratio=:equal, clim=common_clim)
    p2 = heatmap(x, y, mar_map'; title="MAR", xlabel="x", ylabel="y", colorbar_title=spec.units, aspect_ratio=:equal, clim=common_clim)
    p3 = heatmap(x, y, diff'; title="Chion - MAR", xlabel="x", ylabel="y", colorbar_title=spec.units, aspect_ratio=:equal, clim=diff_clim, c=:balance)
    savefig(plot(p1, p2, p3; layout=(1, 3), size=(1500, 430), plot_title="$(spec.title) period map"), path)
end

function save_monthly_map_plot(path, spec, x, y, months, chion_monthly, mar_monthly)
    common_clim = finite_extrema(chion_monthly, mar_monthly)
    diff = chion_monthly .- mar_monthly
    diff_abs = maximum(abs, filter(isfinite, vec(diff)); init=0.0)
    diff_clim = diff_abs > 0 ? (-diff_abs, diff_abs) : (-1.0, 1.0)
    panels = Any[]
    for k in eachindex(months)
        month = months[k]
        push!(panels, heatmap(
            x,
            y,
            chion_monthly[k, :, :]';
            title="M$(month) Chion",
            aspect_ratio=:equal,
            ticks=false,
            colorbar=false,
            clim=common_clim,
        ))
        push!(panels, heatmap(
            x,
            y,
            mar_monthly[k, :, :]';
            title="M$(month) MAR",
            aspect_ratio=:equal,
            ticks=false,
            colorbar=false,
            clim=common_clim,
        ))
        push!(panels, heatmap(
            x,
            y,
            diff[k, :, :]';
            title="M$(month) Chion - MAR",
            aspect_ratio=:equal,
            ticks=false,
            colorbar=false,
            clim=diff_clim,
            c=:balance,
        ))
    end
    subtitle = @sprintf(
        "%s monthly maps; Chion/MAR clim %.4g..%.4g %s, diff clim %.4g..%.4g %s",
        spec.title,
        common_clim[1],
        common_clim[2],
        spec.units,
        diff_clim[1],
        diff_clim[2],
        spec.units,
    )
    savefig(plot(
        panels...;
        layout=(length(months), 3),
        size=(1500, max(260, 230 * length(months))),
        plot_title=subtitle,
        margin=2 * Plots.mm,
    ), path)
end

function save_scatter_plot(path, spec, chion_map, mar_map, weights, stats)
    chion_values = Float64[]
    mar_values = Float64[]
    @inbounds for I in eachindex(chion_map, mar_map, weights)
        if weights[I] > 0 && isfinite(chion_map[I]) && isfinite(mar_map[I])
            push!(chion_values, chion_map[I])
            push!(mar_values, mar_map[I])
        end
    end
    isempty(chion_values) && return
    lo = min(minimum(chion_values), minimum(mar_values))
    hi = max(maximum(chion_values), maximum(mar_values))
    p = scatter(mar_values, chion_values;
        label=false, markersize=1.5, markerstrokewidth=0, alpha=0.35,
        xlabel="MAR ($(spec.units))", ylabel="Chion ($(spec.units))",
        title=@sprintf("%s spatial comparison: bias %.3g, RMSE %.3g, r %.3g", spec.title, stats.bias, stats.rmse, stats.corr),
        framestyle=:box)
    plot!(p, [lo, hi], [lo, hi]; label="1:1", lw=2, color=:black)
    savefig(p, path)
end

function paired_spatial_values(chion_map, mar_map, weights)
    chion_values = Float64[]
    mar_values = Float64[]
    @inbounds for I in eachindex(chion_map, mar_map, weights)
        if weights[I] > 0 && isfinite(chion_map[I]) && isfinite(mar_map[I])
            push!(chion_values, chion_map[I])
            push!(mar_values, mar_map[I])
        end
    end
    return chion_values, mar_values
end

function save_monthly_scatter_plot(path, spec, months, chion_monthly, mar_monthly, weights)
    all_chion = Float64[]
    all_mar = Float64[]
    monthly_pairs = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, length(months))
    monthly_stats = Vector{NamedTuple}(undef, length(months))
    for k in eachindex(months)
        chion_values, mar_values = paired_spatial_values(
            view(chion_monthly, k, :, :),
            view(mar_monthly, k, :, :),
            weights,
        )
        monthly_pairs[k] = (chion_values, mar_values)
        monthly_stats[k] = paired_stats(view(chion_monthly, k, :, :), view(mar_monthly, k, :, :), weights)
        append!(all_chion, chion_values)
        append!(all_mar, mar_values)
    end
    isempty(all_chion) && return

    lo = min(minimum(all_chion), minimum(all_mar))
    hi = max(maximum(all_chion), maximum(all_mar))
    if lo == hi
        pad = max(abs(lo), 1.0) * 0.05
        lo -= pad
        hi += pad
    end

    panels = Any[]
    for k in eachindex(months)
        chion_values, mar_values = monthly_pairs[k]
        stats = monthly_stats[k]
        p = scatter(
            mar_values,
            chion_values;
            label=false,
            markersize=1.0,
            markerstrokewidth=0,
            alpha=0.25,
            xlabel="MAR",
            ylabel="Chion",
            title=@sprintf("M%d bias %.3g RMSE %.3g r %.3g", months[k], stats.bias, stats.rmse, stats.corr),
            framestyle=:box,
            xlim=(lo, hi),
            ylim=(lo, hi),
            aspect_ratio=:equal,
        )
        plot!(p, [lo, hi], [lo, hi]; label=false, lw=1.5, color=:black)
        push!(panels, p)
    end

    n = length(panels)
    ncols = min(4, n)
    nrows = cld(n, ncols)
    savefig(plot(
        panels...;
        layout=(nrows, ncols),
        size=(320 * ncols, 300 * nrows),
        plot_title="$(spec.title) monthly spatial scatter ($(spec.units))",
        margin=2 * Plots.mm,
    ), path)
end

function main(args)
    if has_flag(args, "help")
        print_help()
        return
    end

    mar_file = arg_value(args, "mar-file", DEFAULT_MAR_FILE)
    chion_file = arg_value(args, "chion-file", DEFAULT_CHION_FILE)
    output_dir = arg_value(args, "output-dir", DEFAULT_OUTPUT_DIR)
    chion_period = lowercase(arg_value(args, "chion-period", "last"))
    mar_sector = parse(Int, arg_value(args, "mar-sector", "1"))
    mar_runoff_var = uppercase(arg_value(args, "mar-runoff-var", "RU"))
    mar_sublimation_var = uppercase(arg_value(args, "mar-sublimation-var", "SU"))
    mar_albedo_var = uppercase(arg_value(args, "mar-albedo-var", "AL2"))
    mask_threshold = parse(Float64, arg_value(args, "mask-threshold", "50.0"))
    make_maps = !has_flag(args, "no-maps")

    isfile(mar_file) || error("MAR file not found: $(mar_file)")
    isfile(chion_file) || error("Chion file not found: $(chion_file)")
    mkpath(output_dir)
    default(fmt=:png)

    specs = ComparisonSpec[
        ComparisonSpec(:melt, "Melt", "melt", "ME", :single_sector, :sum, "mmWE"),
        ComparisonSpec(:runoff, "Runoff", "runoff", mar_runoff_var, :sector, :sum, "mmWE"),
        ComparisonSpec(:refreezing, "Refreezing", "refreezing", "RZ", :single_sector, :sum, "mmWE"),
        ComparisonSpec(:sublimation, "Sublimation mass loss", "sublimation", mar_sublimation_var, :sector, :sum, "mmWE"),
        ComparisonSpec(:net_ice_sheet_forcing, "Net ice sheet forcing", "smb_ice", "SMB", :sector, :sum, "mmWE"),
        ComparisonSpec(:albedo, "Albedo", "albedo", mar_albedo_var, :sector, :mean, "1"),
    ]

    monthly_rows = NamedTuple[]
    stats_rows = NamedTuple[]
    NCDataset(mar_file) do mar_ds
        NCDataset(chion_file) do chion_ds
            validate_chion_smb_variables(chion_ds, chion_period)
            weights = weight_matrix(mar_ds, chion_ds, mask_threshold)
            nx, ny = grid_shape(chion_ds)
            x = haskey(chion_ds, "x") ? clean_array(Array(chion_ds["x"][:])) : collect(1:nx)
            y = haskey(chion_ds, "y") ? clean_array(Array(chion_ds["y"][:])) : collect(1:ny)

            for spec in specs
                if !haskey(chion_ds, spec.chion_var)
                    @warn "Skipping $(spec.key): Chion variable $(spec.chion_var) is missing."
                    continue
                end
                if !haskey(mar_ds, spec.mar_var)
                    @warn "Skipping $(spec.key): MAR variable $(spec.mar_var) is missing."
                    continue
                end

                println("Comparing $(spec.title): Chion $(spec.chion_var) vs MAR $(spec.mar_var)")
                chion_monthly, months = read_chion_monthly(chion_ds, spec.chion_var, chion_period)
                mar_monthly = aggregate_mar_monthly(mar_ds, spec, mar_sector)
                mar_selected = repeat_mar_to_months(mar_monthly, months)
                apply_mask!(chion_monthly, weights)
                apply_mask!(mar_selected, weights)

                chion_series = weighted_series(chion_monthly, weights)
                mar_series = weighted_series(mar_selected, weights)
                for k in eachindex(months)
                    push!(monthly_rows, (
                        variable=String(spec.key),
                        month_index=k,
                        month_of_year=months[k],
                        chion_mean=chion_series[k],
                        mar_mean=mar_series[k],
                        bias=chion_series[k] - mar_series[k],
                        units=spec.units,
                    ))
                end

                save_timeseries_plot(
                    joinpath(output_dir, "timeseries_$(spec.key).png"),
                    spec,
                    months,
                    chion_series,
                    mar_series,
                )

                chion_map = reduce_period(chion_monthly, spec.temporal)
                mar_map = reduce_period(mar_selected, spec.temporal)
                stats = paired_stats(chion_map, mar_map, weights)
                push!(stats_rows, (variable=String(spec.key), stats=stats, units=spec.units))

                if make_maps
                    save_map_plot(joinpath(output_dir, "maps_$(spec.key).png"), spec, x, y, chion_map, mar_map)
                    save_monthly_map_plot(
                        joinpath(output_dir, "monthly_maps_$(spec.key).png"),
                        spec,
                        x,
                        y,
                        months,
                        chion_monthly,
                        mar_selected,
                    )
                    save_scatter_plot(joinpath(output_dir, "scatter_$(spec.key).png"), spec, chion_map, mar_map, weights, stats)
                    save_monthly_scatter_plot(
                        joinpath(output_dir, "monthly_scatter_$(spec.key).png"),
                        spec,
                        months,
                        chion_monthly,
                        mar_selected,
                        weights,
                    )
                end
            end
        end
    end

    isempty(monthly_rows) && error("No comparison variables were found in the selected Chion/MAR files; not writing empty output tables.")
    write_monthly_csv(joinpath(output_dir, "monthly_domain_means.csv"), monthly_rows)
    write_stats_csv(joinpath(output_dir, "period_spatial_stats.csv"), stats_rows)

    println("Wrote comparison outputs to $(abspath(output_dir))")
    println("Main tables:")
    println("  ", joinpath(output_dir, "monthly_domain_means.csv"))
    println("  ", joinpath(output_dir, "period_spatial_stats.csv"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
