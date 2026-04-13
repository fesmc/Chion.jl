#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Dates
using HDF5
import Plots
using Printf
using Statistics
using Chion

const DEFAULT_DATA_DIR = joinpath(@__DIR__, "..", "..", "data", "ESM-SnowMIP_all")
const DEFAULT_OUT_DIR = joinpath(@__DIR__, "..", "plots")

arg_value(args, name, default="") = something(findfirst(a -> startswith(a, "--$name="), args), 0) > 0 ?
    split(args[findfirst(a -> startswith(a, "--$name="), args)], "=", limit=2)[2] : default
has_flag(args, name) = any(==("--$name"), args)
fmtf(x) = isfinite(x) ? @sprintf("%.6f", x) : "NaN"

function print_help()
    println("Usage:")
    println("  julia --project=. examples/scripts/validate_esm_snowmip_site.jl [options]")
    println()
    println("Inputs:  --met-nc=PATH --obs-nc=PATH")
    println("         or --site=snb --forcing=insitu|gswp3c [--data-dir=PATH]")
    println("Vars:    --met-time-var=time --obs-time-var=time --tair-var=Tair --rain-var=Rainf")
    println("         --snow-var=Snowf --sw-var=SWdown --lw-var=LWdown --wind-var=Wind")
    println("         --obs-depth-var=auto --obs-swe-var=auto")
    println("Model:   --densification=bessi|htessel --albedo=dynamic|constant|legacy|bessi")
    println("         --fresh-snow-density=constant|parameterized|bessi|htessel")
    println("         --initial-density=300 --ntot=60 --no-init-from-obs --no-radiation-forcing")
    println("Output:  --out-dir=PATH --slug=NAME --no-plot --help")
end

function site_files(data_dir, site, forcing)
    site = lowercase(strip(String(site)))
    forcing = lowercase(strip(String(forcing)))
    met = obs = nothing
    for path in readdir(data_dir; join=true)
        name = lowercase(basename(path))
        occursin(Regex("^met_$(forcing)_$(site)_\\d{4}_\\d{4}\\.nc\$"), name) && (met = path)
        occursin(Regex("^obs_insitu_$(site)_\\d{4}_\\d{4}\\.nc\$"), name) && (obs = path)
    end
    isnothing(met) && error("No met file found for site=$(site), forcing=$(forcing) in $(abspath(data_dir)).")
    isnothing(obs) && error("No obs file found for site=$(site) in $(abspath(data_dir)).")
    return (met=met, obs=obs)
end

read_names(path) = h5open(path, "r") do f
    Set(String(k) for k in keys(f))
end

function read_vector(path, var)
    h5open(path, "r") do f
        haskey(f, var) || error("Variable '$var' not found in $(abspath(path)).")
        dset = f[var]
        vals = vec(Float64.(read(dset)))
        attrs = attributes(dset)
        fill = "_FillValue" in keys(attrs) ? read(attrs["_FillValue"]) : nothing
        fill = fill isa AbstractArray ? Float64(fill[1]) : (isnothing(fill) ? nothing : Float64(fill))
        @inbounds for i in eachindex(vals)
            if !isfinite(vals[i]) || vals[i] <= -9e18 || (!isnothing(fill) && vals[i] == fill)
                vals[i] = NaN
            end
        end
        vals
    end
end

function read_time(path, var)
    vals = read_vector(path, var)
    units = h5open(path, "r") do f
        String(read(attributes(f[var])["units"]))
    end
    m = match(r"(?i)^(seconds|minutes|hours|days)\s+since\s+(.+)$", units)
    isnothing(m) && error("Unsupported time units '$units'.")
    ref_txt = strip(replace(replace(m.captures[2], "T" => " "), "UTC" => ""))
    ref = nothing
    for fmt in (dateformat"y-m-d H:M:S.s", dateformat"y-m-d H:M:S", dateformat"y-m-d H:M", dateformat"y-m-d")
        ref = tryparse(DateTime, ref_txt, fmt)
        !isnothing(ref) && break
    end
    isnothing(ref) && error("Could not parse time reference '$ref_txt'.")
    scale = Dict("days" => 86_400_000.0, "hours" => 3_600_000.0, "minutes" => 60_000.0, "seconds" => 1_000.0)[lowercase(m.captures[1])]
    out = Vector{DateTime}(undef, length(vals))
    @inbounds for i in eachindex(vals)
        out[i] = isfinite(vals[i]) ? ref + Millisecond(round(Int, vals[i] * scale)) : (i == 1 ? ref : out[i - 1])
    end
    out
end

function align_obs(target_dates, obs_dates, obs_depth, obs_swe)
    lookup = Dict(d => i for (i, d) in pairs(obs_dates))
    depth = fill(NaN, length(target_dates))
    swe = fill(NaN, length(target_dates))
    for i in eachindex(target_dates)
        j = get(lookup, target_dates[i], 0)
        j == 0 && continue
        depth[i] = obs_depth[j]
        swe[i] = obs_swe[j]
    end
    return depth, swe
end

function carry_forward(v, fallback; nonnegative=false)
    out = copy(v)
    last = fallback
    @inbounds for i in eachindex(out)
        ok = isfinite(out[i]) && (!nonnegative || out[i] >= 0)
        last = ok ? out[i] : last
        out[i] = ok ? out[i] : last
    end
    out
end

function pick_name(names, candidates)
    for name in candidates
        name in names && return name
    end
    return nothing
end

function load_inputs(cfg)
    paths = isempty(cfg.met_nc) ? site_files(cfg.data_dir, cfg.site, cfg.forcing) : (met=abspath(cfg.met_nc), obs=abspath(cfg.obs_nc))
    obs_names = read_names(paths.obs)
    obs_depth_var = cfg.obs_depth_var == "auto" ? pick_name(obs_names, ["snd_auto", "snd_can_auto", "snd_gap_auto", "snd_gap1_auto", "snd_gap2_auto", "snd_man"]) : cfg.obs_depth_var
    obs_swe_var = cfg.obs_swe_var == "auto" ? pick_name(obs_names, ["snw_auto", "snw_man"]) : cfg.obs_swe_var
    isnothing(obs_depth_var) && error("Could not resolve observation depth variable.")
    met_dates = read_time(paths.met, cfg.met_time_var)
    obs_dates = read_time(paths.obs, cfg.obs_time_var)
    obs_depth, obs_swe = align_obs(
        met_dates,
        obs_dates,
        read_vector(paths.obs, obs_depth_var),
        isnothing(obs_swe_var) ? fill(NaN, length(obs_dates)) : read_vector(paths.obs, obs_swe_var),
    )
    return (
        met_path = paths.met,
        obs_path = paths.obs,
        dates = met_dates,
        tair = read_vector(paths.met, cfg.tair_var),
        rain = read_vector(paths.met, cfg.rain_var),
        snow = read_vector(paths.met, cfg.snow_var),
        sw = read_vector(paths.met, cfg.sw_var),
        lw = read_vector(paths.met, cfg.lw_var),
        wind = read_vector(paths.met, cfg.wind_var),
        obs_depth = obs_depth,
        obs_swe = obs_swe,
        obs_depth_var = obs_depth_var,
        obs_swe_var = obs_swe_var,
    )
end

parse_densification(text) = begin
    scheme = Symbol(lowercase(strip(String(text))))
    scheme in (:bessi, :htessel) || error("Unsupported densification scheme '$text'.")
    scheme
end

function build_case(inp, cfg, physics)
    tair = carry_forward(inp.tair, 268.0)
    wind = carry_forward(inp.wind, 5.0; nonnegative=true)
    snow = max.(replace(copy(inp.snow), NaN => 0.0), 0.0)
    rain = max.(replace(copy(inp.rain), NaN => 0.0), 0.0)
    sw = max.(replace(copy(inp.sw), NaN => 0.0), 0.0)
    lw = replace(copy(inp.lw), NaN => 0.0)
    dt_days = [i < length(inp.dates) ? Dates.value(inp.dates[i + 1] - inp.dates[i]) / 86_400_000.0 : Dates.value(inp.dates[i] - inp.dates[i - 1]) / 86_400_000.0 for i in eachindex(inp.dates)]
    init_mass = (!cfg.no_init_from_obs && isfinite(inp.obs_depth[1]) && inp.obs_depth[1] > 0) ? cfg.initial_density * inp.obs_depth[1] : 0.0
    definition = Chion.prescribed_case(
        physics = physics,
        ntot = cfg.ntot,
        nx = 1,
        ny = 1,
        dt_days = dt_days,
        air_temperature_c = tair .- 273.15,
        snowfall_mm_day = snow .* 86_400.0,
        rainfall_mm_day = rain .* 86_400.0,
        shortwave_down = sw,
        wind_speed = wind,
        q_lw_down = lw,
        has_q_lw_down = isfinite.(inp.lw),
        time_values = inp.dates,
        initial_surface_mass = init_mass,
        initial_density = cfg.initial_density,
        initial_temperature_c = tair[1] - 273.15,
        input_label = inp.met_path,
    )
    Chion.build_case(
        definition;
        name = "ESM-SnowMIP validation",
        backend = :cpu,
        write_outputs = false,
        write_netcdf = false,
        cycles = 1,
    )
end

function current_state(domain)
    state = Chion.get_state(domain, 1)
    depth = max(state["total_thickness"], 0.0)
    return (
        depth = depth,
        swe = state["total_mass"],
        rho = depth > 0 ? state["total_mass"] / depth : NaN,
        cover = state["snow_cover"],
    )
end

function run_timeseries(case; use_radiation=true)
    domain = deepcopy(case.definition.domain)
    forcing = case.definition.forcing
    n = length(forcing.time_values)
    sim = (
        depth = fill(NaN, n),
        swe = fill(NaN, n),
        rho = fill(NaN, n),
        cover = fill(NaN, n),
        tair = fill(NaN, n),
        wind = fill(NaN, n),
        sw = fill(NaN, n),
        lw = fill(NaN, n),
        rain = fill(NaN, n),
        snow = fill(NaN, n),
        cum_snow = zeros(Float64, n),
    )
    ws = Chion.StepWorkspace(domain)
    s0 = current_state(domain)
    sim.depth[1] = s0.depth; sim.swe[1] = s0.swe; sim.rho[1] = s0.rho; sim.cover[1] = s0.cover
    last_tair = isfinite(forcing.air_temperature[1, 1]) ? forcing.air_temperature[1, 1] : 268.0
    last_wind = isfinite(forcing.wind_speed[1, 1]) ? max(forcing.wind_speed[1, 1], 0.0) : 5.0
    for i in 1:n
        T = forcing.air_temperature[1, i]; T = isfinite(T) ? (last_tair = T) : last_tair
        W = forcing.wind_speed[1, i]; W = (isfinite(W) && W >= 0) ? (last_wind = W) : last_wind
        S = isfinite(forcing.snowfall_rate[1, i]) ? max(forcing.snowfall_rate[1, i], 0.0) : 0.0
        R = isfinite(forcing.rainfall_rate[1, i]) ? max(forcing.rainfall_rate[1, i], 0.0) : 0.0
        sw_i = use_radiation && isfinite(forcing.shortwave_down[1, i]) ? max(forcing.shortwave_down[1, i], 0.0) : NaN
        lw_i = use_radiation && forcing.has_q_lw_down[1, i] && isfinite(forcing.q_lw_down[1, i]) ? forcing.q_lw_down[1, i] : NaN
        sim.tair[i] = T - 273.15; sim.wind[i] = W; sim.snow[i] = S; sim.rain[i] = R; sim.sw[i] = sw_i; sim.lw[i] = lw_i
        i == 1 && continue
        dt = forcing.dt_days[i - 1]
        sim.cum_snow[i] = sim.cum_snow[i - 1] + S * dt * 86_400.0
        Chion.step!(domain, 1, T, S + R, dt; workspace=ws, snowfall_rate=S, rainfall_rate=R, wind_speed=W, shortwave_down=isfinite(sw_i) ? sw_i : nothing, q_lw_down=isfinite(lw_i) ? lw_i : nothing)
        st = current_state(domain)
        sim.depth[i] = st.depth; sim.swe[i] = st.swe; sim.rho[i] = st.rho; sim.cover[i] = st.cover
    end
    sim
end

function metric(obs, sim)
    valid = isfinite.(obs) .& isfinite.(sim)
    n = count(valid)
    n == 0 && return (n=0, bias=NaN, mae=NaN, rmse=NaN, corr=NaN)
    delta = sim[valid] .- obs[valid]
    corr = n > 1 && std(obs[valid]) > 0 && std(sim[valid]) > 0 ? cor(obs[valid], sim[valid]) : NaN
    return (n=n, bias=mean(delta), mae=mean(abs.(delta)), rmse=sqrt(mean(delta .^ 2)), corr=corr)
end

bulk_density(swe, depth) = [isfinite(m) && isfinite(h) && h > 0 ? m / h : NaN for (m, h) in zip(swe, depth)]

function write_plot(path, slug, dates, obs_depth, sim_depth, obs_swe, sim_swe, sim, metrics)
    P = Plots
    try P.default(fmt=:svg) catch end
    obs_rho = bulk_density(obs_swe, obs_depth)
    p1 = P.scatter(dates, obs_depth; label="obs", color=:steelblue, ms=2.5, markerstrokewidth=0, ylabel="m", title=@sprintf("Snow Depth RMSE=%.3f", metrics.depth.rmse), framestyle=:box)
    P.plot!(p1, dates, sim_depth; label="sim", color=:firebrick, lw=2)
    p2 = P.scatter(dates, obs_swe; label="obs", color=:steelblue, ms=2.5, markerstrokewidth=0, ylabel="kg m^-2", title=@sprintf("Snow SWE RMSE=%.3f", metrics.swe.rmse), framestyle=:box)
    P.plot!(p2, dates, sim_swe; label="sim", color=:firebrick, lw=2)
    p3 = P.scatter(dates, obs_rho; label="obs", color=:steelblue, ms=2.5, markerstrokewidth=0, ylabel="kg m^-3", title="Bulk Density", framestyle=:box)
    P.plot!(p3, dates, sim.rho; label="sim", color=:firebrick, lw=2)
    p4 = P.scatter(obs_depth, sim_depth; label="samples", color=:firebrick, alpha=0.35, xlabel="obs depth [m]", ylabel="sim depth [m]", title="Depth Scatter", framestyle=:box)
    vals = [filter(isfinite, obs_depth); filter(isfinite, sim_depth)]
    lo, hi = isempty(vals) ? (0.0, 1.0) : (minimum(vals), maximum(vals))
    P.plot!(p4, [lo, hi], [lo, hi]; label="1:1", color=:black, ls=:dash)
    p5 = P.plot(dates, sim.tair; label="Tair", color=:teal, lw=2, ylabel="C", title="Air Temperature", framestyle=:box)
    p6 = P.plot(dates, sim.sw; label="SWdown", color=:purple, lw=2, ylabel="W m^-2", title="Radiation", framestyle=:box)
    P.plot!(p6, dates, sim.lw; label="LWdown", color=:dodgerblue, lw=2)
    P.savefig(P.plot(p1, p2, p3, p4, p5, p6; layout=(3, 2), size=(1500, 1100), plot_title=slug), path)
    path
end

function write_outputs(out_dir, slug, dates, obs_depth, sim_depth, obs_swe, sim_swe, sim, metrics; plot_on=true)
    mkpath(out_dir)
    csv = joinpath(out_dir, slug * "_timeseries.csv")
    txt = joinpath(out_dir, slug * "_metrics.txt")
    svg = joinpath(out_dir, slug * "_evolution.svg")
    open(csv, "w") do io
        println(io, "date,obs_depth_m,sim_depth_m,obs_swe_kgm2,sim_swe_kgm2,sim_bulk_density_kgm3,sim_snow_cover,tair_C,swdown_Wm2,lwdown_Wm2,rainf_kgm2s,snowf_kgm2s,cum_snow_mmwe")
        for i in eachindex(dates)
            println(io, join([
                Dates.format(dates[i], dateformat"yyyy-mm-ddTHH:MM:SS"),
                fmtf(obs_depth[i]), fmtf(sim_depth[i]), fmtf(obs_swe[i]), fmtf(sim_swe[i]),
                fmtf(sim.rho[i]), fmtf(sim.cover[i]), fmtf(sim.tair[i]), fmtf(sim.sw[i]),
                fmtf(sim.lw[i]), fmtf(sim.rain[i]), fmtf(sim.snow[i]), fmtf(sim.cum_snow[i]),
            ], ","))
        end
    end
    open(txt, "w") do io
        println(io, "General API  : prescribed_case + build_case + step!")
        println(io, @sprintf("Snow depth : n=%d bias=%.5f mae=%.5f rmse=%.5f corr=%.5f", metrics.depth.n, metrics.depth.bias, metrics.depth.mae, metrics.depth.rmse, metrics.depth.corr))
        println(io, @sprintf("Snow SWE   : n=%d bias=%.5f mae=%.5f rmse=%.5f corr=%.5f", metrics.swe.n, metrics.swe.bias, metrics.swe.mae, metrics.swe.rmse, metrics.swe.corr))
    end
    plot_on && write_plot(svg, slug, dates, obs_depth, sim_depth, obs_swe, sim_swe, sim, metrics)
    return (csv=csv, txt=txt, svg=svg)
end

function main(args)
    has_flag(args, "help") && return print_help()
    cfg = (
        met_nc = arg_value(args, "met-nc"),
        obs_nc = arg_value(args, "obs-nc"),
        site = arg_value(args, "site", "snb"),
        forcing = arg_value(args, "forcing", "insitu"),
        data_dir = arg_value(args, "data-dir", DEFAULT_DATA_DIR),
        met_time_var = arg_value(args, "met-time-var", "time"),
        obs_time_var = arg_value(args, "obs-time-var", "time"),
        tair_var = arg_value(args, "tair-var", "Tair"),
        rain_var = arg_value(args, "rain-var", "Rainf"),
        snow_var = arg_value(args, "snow-var", "Snowf"),
        sw_var = arg_value(args, "sw-var", "SWdown"),
        lw_var = arg_value(args, "lw-var", "LWdown"),
        wind_var = arg_value(args, "wind-var", "Wind"),
        obs_depth_var = arg_value(args, "obs-depth-var", "auto"),
        obs_swe_var = arg_value(args, "obs-swe-var", "auto"),
        initial_density = parse(Float64, arg_value(args, "initial-density", "300.0")),
        ntot = parse(Int, arg_value(args, "ntot", "60")),
        no_init_from_obs = has_flag(args, "no-init-from-obs"),
        no_radiation = has_flag(args, "no-radiation-forcing"),
    )
    xor(isempty(cfg.met_nc), isempty(cfg.obs_nc)) && error("Pass both --met-nc and --obs-nc, or neither.")

    physics = Chion.physics(
        albedo = Symbol(lowercase(arg_value(args, "albedo", "dynamic"))),
        densification = parse_densification(arg_value(args, "densification", arg_value(args, "low-density-densification", "bessi"))),
        fresh_snow_density = Symbol(lowercase(arg_value(args, "fresh-snow-density", "constant"))),
    )
    inp = load_inputs(cfg)
    case = build_case(inp, cfg, physics)
    sim = run_timeseries(case; use_radiation=!cfg.no_radiation)
    metrics = (depth=metric(inp.obs_depth, sim.depth), swe=metric(inp.obs_swe, sim.swe))
    slug = arg_value(args, "slug", isempty(cfg.met_nc) ? "$(cfg.site)_$(cfg.forcing)" : splitext(basename(cfg.met_nc))[1])
    out = write_outputs(arg_value(args, "out-dir", DEFAULT_OUT_DIR), slug, inp.dates, inp.obs_depth, sim.depth, inp.obs_swe, sim.swe, sim, metrics; plot_on=!has_flag(args, "no-plot"))

    println("Validation complete.")
    println("General API : prescribed_case + build_case + step!")
    println("Met file    : $(abspath(inp.met_path))")
    println("Obs file    : $(abspath(inp.obs_path))")
    println("CSV output  : $(abspath(out.csv))")
    println("Metrics file: $(abspath(out.txt))")
    has_flag(args, "no-plot") || println("Plot output : $(abspath(out.svg))")
    println(@sprintf("Snow-depth RMSE = %.5f m (n=%d)", metrics.depth.rmse, metrics.depth.n))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
