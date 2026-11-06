#!/usr/bin/env julia

# Ice-sheet mass budget of a Chion run against MAR, year by year: reads budget_yearly.csv and the
# configuration of a run made with `julia --project -m Chion <config>.toml`, sums the MAR fields
# (SMB, ME, RU, RZ, −SU) over the same columns and writes mar_comparison.csv and .pdf next to it.
#
#   julia --project=. examples/figures/compare_with_mar.jl output/greenland_parameterized
using Chion, NCDatasets, Plots, Statistics, TOML, Dates

const MAR_TERMS = (surface_smb=("SMB", 1), melt=("ME", 1), runoff=("RU", 1), refreezing=("RZ", 1), vapor_mass=("SU", -1))

"""
    mar_comparison(run_dir) -> Dict

Chion and MAR totals (Gt) per record of budget_yearly.csv: keys `:year`, `:period_end`,
`Symbol("chion_\$v")` and `Symbol("mar_\$v")`. Cached in run_dir/mar_comparison.csv.
"""
function mar_comparison(run_dir)
    budget = joinpath(run_dir, "budget_yearly.csv")
    cache = joinpath(run_dir, "mar_comparison.csv")
    isfile(cache) && mtime(cache) >= mtime(budget) && return _read_comparison(cache)
    cfg = TOML.parsefile(only(filter(endswith(".toml"), readdir(run_dir; join=true))))
    f = cfg["forcing"]
    files = get(f, "files", Dict())
    forcing = Forcing(f["path"]; files=Dict(parse(Int, k) => v for (k, v) in files), area_file=get(f, "area_file", nothing))
    grid, area = forcing_grid(forcing), forcing_area(forcing)
    lines = readlines(budget)
    header = split(lines[1], ',')
    rows = [split(l, ',') for l in lines[2:end]]
    years = [year(DateTime(r[1])) for r in rows]
    cached = Dict{Tuple{Int, String}, Float64}()               # equilibrium runs repeat years
    function mar_total(y, name, sign)
        get!(cached, (y, name)) do
            path = get(files, string(y), replace(f["path"], "{year}" => y))
            NCDataset(path) do ds
                haskey(ds, name) || return NaN                 # e.g. no RZ in MAR Antarctica
                times = DateTime.(ds["TIME"][:])
                cols = grid.is .+ (grid.js .- 1) .* length(grid.x)
                x_name, y_name = Chion._mar_xy_names(ds)
                x = read_forcing_columns(ds, name, cols, length(times), length(grid.y), length(grid.x); x_name, y_name)
                sign * sum(vec(sum(x .* reshape(Chion.infer_dt_days(times), 1, :); dims=2)) .* area) * 1e-6
            end
        end
    end
    open(cache, "w") do io
        println(io, "year,period_end,", join(("chion_$(k)_gt,mar_$(k)_gt" for k in keys(MAR_TERMS)), ','))
        for (i, y) in enumerate(years)
            chion(k) = parse(Float64, rows[i][findfirst(==("$k (Gt)"), header)])
            println(io, y, ',', rows[i][2], ',', join(("$(chion(k)),$(mar_total(y, MAR_TERMS[k]...))" for k in keys(MAR_TERMS)), ','))
        end
    end
    return _read_comparison(cache)
end

function _read_comparison(path)
    lines = readlines(path)
    header = split(lines[1], ',')
    rows = [split(l, ',') for l in lines[2:end]]
    t = Dict{Symbol, Any}(:year => [parse(Int, r[1]) for r in rows], :period_end => [DateTime(r[2]) for r in rows])
    for (j, h) in enumerate(header)
        j > 2 && (t[Symbol(replace(h, "_gt" => ""))] = [parse(Float64, r[j]) for r in rows])
    end
    return t
end

"""Plot Chion against MAR for every budget term; incomplete years are left out."""
function plot_mar_comparison(run_dir, t=mar_comparison(run_dir))
    full = [Date(e) >= Date(y, 12, 30) for (y, e) in zip(t[:year], t[:period_end])]
    panels = map(collect(keys(MAR_TERMS))) do k
        c, m = t[Symbol("chion_$k")], t[Symbol("mar_$k")]
        ok = isfinite.(m) .& full
        p = plot(t[:year][full], [m[full] c[full]]; label=["MAR" "Chion"], c=[:black "#2a78d6"], lw=1.5, framestyle=:box,
            ylabel="$(replace(String(k), '_' => ' ')) (Gt yr⁻¹)", legend=k == :surface_smb ? :bottomleft : false)
        any(ok) && title!(p, "bias $(round(mean(c[ok] .- m[ok]); digits=1)), " *
            "RMSE $(round(sqrt(mean((c[ok] .- m[ok]) .^ 2)); digits=1)) Gt yr⁻¹"; titlefontsize=9)
        p
    end
    savefig(plot(panels...; layout=(3, 2), size=(1100, 900), left_margin=5Plots.mm), joinpath(run_dir, "mar_comparison.pdf"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_dir = abspath(get(ARGS, 1, "output/greenland_parameterized"))
    plot_mar_comparison(run_dir)
    println("Wrote ", joinpath(run_dir, "mar_comparison.csv"))
end
