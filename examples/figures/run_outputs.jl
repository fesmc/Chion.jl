# Readers for the outputs of the runs in run_paper_experiments.sh (configs and examples/runs):
# budget_yearly.csv, yearly.nc and warming.csv of examples/runs/cmip6_delta.jl, and firn_profiles.nc.
using Chion, NCDatasets, TOML, Dates

"""Grid and cell areas (km²) of a run, from the configuration file kept in its directory."""
const _RUN_GRIDS = Dict{String, Any}()
function run_grid(dir)
    cfg_path = only(filter(endswith(".toml"), readdir(dir; join=true)))
    f = TOML.parsefile(cfg_path)["forcing"]
    get!(_RUN_GRIDS, f["path"]) do
        forcing = Forcing(f["path"]; area_file=get(f, "area_file", nothing))
        (grid=forcing_grid(forcing), area=forcing_area(forcing))
    end
end

"""A budget CSV (`period_start, period_end, var (Gt), …`) as Dict(:year, :<var>_gt, …)."""
function read_budget(path)
    lines = readlines(path)
    header = split(lines[1], ',')
    rows = [split(l, ',') for l in lines[2:end]]
    t = Dict{Symbol, Any}(:year => [year(DateTime(r[1])) for r in rows])
    for (j, h) in enumerate(header)
        m = match(r"^(\w+) \(Gt\)$", h)
        m === nothing || (t[Symbol(m.captures[1], "_gt")] = [parse(Float64, r[j]) for r in rows])
    end
    return t
end

"""
Yearly totals of a cmip6_delta.jl run (`nothing` if incomplete), with the keys of the earlier
projection_budgets.csv: :year, :surface_smb_gt, :runoff_gt, …, :mean_delta_t_jja/_annual/_djf and,
for uniform-warming runs, the phase in `t[:_text][:phase]` ("control" or "warm").
"""
function read_projection(dir; min_years=120)
    if isfile(joinpath(dir, "warm_budget_yearly.csv"))
        c, w = read_budget(joinpath(dir, "control_budget_yearly.csv")), read_budget(joinpath(dir, "warm_budget_yearly.csv"))
        t = Dict{Symbol, Any}(k => vcat(c[k], w[k]) for k in keys(c))
        t[:_text] = Dict(:phase => vcat(fill("control", length(c[:year])), fill("warm", length(w[:year]))))
        return t
    end
    isfile(joinpath(dir, "warming.csv")) || (isfile(joinpath(dir, "budget_yearly.csv")) && min_years == 0) ||
        return nothing
    t = read_budget(joinpath(dir, "budget_yearly.csv"))
    length(t[:year]) >= min_years || return nothing
    if isfile(joinpath(dir, "warming.csv"))
        lines = readlines(joinpath(dir, "warming.csv"))
        w = [parse.(Float64, split(l, ',')) for l in lines[2:end]]
        t[:mean_delta_t_jja] = getindex.(w, 3); t[:mean_delta_t_annual] = getindex.(w, 4); t[:mean_delta_t_djf] = getindex.(w, 5)
    end
    if isfile(joinpath(dir, "yearly.nc"))
        cols = read_projection_columns(dir; min_years=0)
        haskey(cols, :fac) && (t[:fac_km3] = vec(sum(cols.fac .* cols.area; dims=1)) .* 1e-3)   # m × km² → km³
    end
    return t
end

"""Per-column yearly fields of a cmip6_delta.jl run: (years, i, j, area, smb, runoff, melt, refreezing, fac)."""
function read_projection_columns(dir; min_years=120)
    path = joinpath(dir, "yearly.nc")
    isfile(path) || return nothing
    g = run_grid(dir)
    NCDataset(path) do nc
        years = year.(nc["time"][:])
        length(years) >= min_years || return nothing
        field(name) = (x = nc[name].var[:, :, :]; Float64[x[g.grid.is[c], g.grid.js[c], k] for c in eachindex(g.grid.is), k in eachindex(years)])
        names = (smb="surface_smb", runoff="runoff", melt="melt", refreezing="refreezing", fac="firn_air_content")
        present = filter(k -> haskey(nc, names[k]), keys(names))
        return (; years, i=g.grid.is, j=g.grid.js, area=g.area, (k => field(names[k]) for k in present)...)
    end
end

"""Static MAR field `name` (e.g. LAT, LON, SH) on the columns of a run."""
function run_static(dir, name)
    g = run_grid(dir).grid
    f = TOML.parsefile(only(filter(endswith(".toml"), readdir(dir; join=true))))["forcing"]
    forcing = Forcing(f["path"])
    path = forcing.source.files[first(forcing_years(forcing))]
    field = NCDataset(path) do ds
        raw, dims = Chion._read_variable_data(ds, name)
        Chion._as_y_x(raw, dims, length(g.y), length(g.x), name)
    end
    return [Float64(field[g.js[c], g.is[c]]) for c in eachindex(g.js)]
end

"""Firn profiles of a configuration run, in the layout of the earlier firn_profiles.nc:
(lat, lon, elev, time, zr, zt, rho, temp) with rho (column, depth bin, snapshot) on the bin centres
zr and temp (column, depth, snapshot) at zt = 0.5, 1, …, 20 m; `time` is the snapshot time."""
function read_firn_profiles(dir)
    g = run_grid(dir).grid
    NCDataset(joinpath(dir, "firn_profiles.nc")) do ds
        nt = length(ds["time"])
        column(x) = Float32[coalesce(x[g.is[c], g.js[c], k, t], NaN32) for c in eachindex(g.is), k in axes(x, 3), t in 1:nt]
        zt = Float64.(ds["depth"][:]); keep = findall(>(0), zt)                    # temperature below the surface
        (; lat=run_static(dir, "LAT"), lon=run_static(dir, "LON"), elev=run_static(dir, "SH"),
           time=DateTime.(ds["time_bnds"][2, :]), zr=Float64.(ds["depth_bin"][:]), zt=zt[keep],
           rho=column(ds["density_profile"].var[:, :, :, :]), temp=column(ds["temperature_profile"].var[:, :, keep, :]))
    end
end
