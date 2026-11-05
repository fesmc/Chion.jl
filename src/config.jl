"""
Run an `Experiment` from a TOML file:

    julia --project -m Chion greenland.toml          # Julia ≥ 1.11
    julia --project -e 'using Chion; run_config("greenland.toml")'

See `examples/configs/` for complete files. Relative output paths are taken relative to
`output_dir` (default: the directory of the TOML file).
"""

using TOML

const _SYMBOL_KWARGS = (:albedo, :seb_scheme, :turbulent_flux_scheme, :densification, :fresh_snow_density,
    :longwave_scheme, :low_density_densification)

function _config_years(v)
    v isa AbstractString && (p = parse.(Int, split(v, ':')); return p[1]:p[end])
    v isa Integer && return v:v
    v isa AbstractVector && length(v) == 2 && return Int(v[1]):Int(v[2])
    return collect(Int, v)
end

function _config_frequency(v)
    v isa AbstractString || error("`frequency` must be a string, e.g. \"monthly\" or \"5 days\".")
    s = lowercase(strip(v))
    s in ("step", "daily", "monthly", "yearly") && return Symbol(s)
    s == "annual" && return :yearly
    m = match(r"^(\d+)\s*(hour|day|week|month|year)s?$", s)
    m === nothing && error("Cannot read the output frequency \"$v\".")
    n = parse(Int, m.captures[1])
    unit = m.captures[2]
    return unit == "hour" ? Hour(n) : unit == "day" ? Day(n) : unit == "week" ? Week(n) : unit == "month" ? Month(n) : Year(n)
end

function _config_model_kwargs(d)
    kw = Pair{Symbol, Any}[]
    for (k, v) in d
        k == "type" && continue
        key = Symbol(k)
        val = key == :number_type ? (v == "Float32" ? Float32 : Float64) :
            v isa AbstractString && key in _SYMBOL_KWARGS ? Symbol(v) :
            v isa AbstractVector ? Tuple(v) : v
        push!(kw, key => val)
    end
    return kw
end

function _config_forcing(d, base)
    path = d["path"]
    isabspath(path) || (path = joinpath(base, path))
    kw = Pair{Symbol, Any}[]
    haskey(d, "format") && push!(kw, :format => Symbol(d["format"]))
    haskey(d, "surface_fluxes") && push!(kw, :surface_fluxes => Symbol(d["surface_fluxes"]))
    haskey(d, "files") && push!(kw, :files => Dict(parse(Int, k) => v for (k, v) in d["files"]))
    haskey(d, "mask") && push!(kw, :mask => (d["mask"][1], Float64(d["mask"][2])))
    haskey(d, "wind") && push!(kw, :wind => d["wind"])
    haskey(d, "area") && push!(kw, :area => d["area"] isa AbstractString ? d["area"] : nothing)
    haskey(d, "area_file") && push!(kw, :area_file => d["area_file"])
    haskey(d, "variables") && push!(kw, :variables => NamedTuple(Symbol(k) => v for (k, v) in d["variables"]))
    for (k, v) in get(d, "options", Dict())
        push!(kw, Symbol(k) => v)
    end
    return Forcing(path; kw...)
end

function _config_spinup(d, main, base)
    d === nothing && return NoSpinup()
    t = lowercase(get(d, "type", "cycled"))
    t == "none" && return NoSpinup()
    t == "cycled" && return CycledSpinup(_config_years(d["years"]); cycles=get(d, "cycles", 5))
    if t == "climatology"
        # The settings of [forcing] (format, surface fluxes, …) with the path of the climatology file.
        sub = merge(filter(kv -> first(kv) != "files", main), get(d, "forcing", Dict{String, Any}()))
        haskey(d, "path") && (sub["path"] = d["path"])
        return ClimatologySpinup(_config_forcing(sub, base); years=get(d, "years", 200))
    end
    error("Unknown spin-up type \"$t\"; use none, cycled or climatology.")
end

function _config_output(d, outdir)
    file = d["file"]
    isabspath(file) || (file = joinpath(outdir, file))
    vars = get(d, "variables", "default")
    kw = Pair{Symbol, Any}[:variables => vars isa AbstractString ? Symbol(vars) : Symbol.(vars),
        :frequency => _config_frequency(get(d, "frequency", "monthly")),
        :aggregate => Symbol(get(d, "aggregate", "columns")), :kind => Symbol(get(d, "kind", "mean"))]
    haskey(d, "period") && push!(kw, :period => Tuple(Int.(d["period"])))
    # depths = {start = 0, stop = 20, step = 0.5} or a list of values
    haskey(d, "depths") && push!(kw, :depths => (dd = d["depths"]; dd isa AbstractDict ?
        range(dd["start"], dd["stop"]; step=dd["step"]) : Float64.(dd)))
    haskey(d, "snow_cover_threshold") && push!(kw, :snow_cover_threshold => d["snow_cover_threshold"])
    return Output(file; kw...)
end

"""Value of a `--set` assignment: TOML syntax (`0.8`, `true`, `[1940, 1969]`, `"x"`), else a string."""
function _set_value(v::AbstractString)
    try
        return TOML.parse("v = " * v)["v"]
    catch
        return String(v)
    end
end

"""Apply `"model.alpha_dry" => 0.84`-style overrides to a parsed configuration."""
function _apply_overrides!(cfg, set)
    for (key, value) in set
        parts = split(String(key), '.')
        d = cfg
        for p in parts[1:end-1]
            d = get!(d, p, Dict{String, Any}())
        end
        d[parts[end]] = value isa AbstractString ? _set_value(value) : value
    end
    return cfg
end

"""
    run_config(path; backend=nothing, output_dir=nothing, set=(), dry_run=false) -> ExperimentResult

Build the `Experiment` described by the TOML file `path` and run it. `backend` and `output_dir`
override the values in the file; `set` overrides any entry, e.g.
`set=["model.diurnal_shortwave_max_substeps" => 4, "forcing.surface_fluxes" => "prescribed"]`.
`dry_run=true` returns the `Experiment` without running it.
"""
function run_config(path::AbstractString; backend=nothing, output_dir=nothing, set=(), dry_run::Bool=false)
    cfg = _apply_overrides!(TOML.parsefile(path), set)
    base = dirname(abspath(path))
    outdir = output_dir === nothing ? get(cfg, "output_dir", base) : output_dir
    isabspath(outdir) || (outdir = joinpath(output_dir === nothing ? base : pwd(), outdir))
    forcing = _config_forcing(cfg["forcing"], base)
    mcfg = get(cfg, "model", Dict{String, Any}())
    mtype = lowercase(get(mcfg, "type", "bessi"))
    builder = mtype == "bessi" ? BESSIModel : mtype == "pdd" ? PDDModel : mtype == "itm" ? ITMModel :
        error("Unknown model type \"$mtype\".")
    model = builder(forcing; _config_model_kwargs(mcfg)...)
    outputs = [_config_output(o, outdir) for o in get(cfg, "output", Any[])]
    exp = Experiment(model, forcing;
        # `repeat = n` runs the years n times (e.g. to record every cycle of an equilibrium run)
        years=haskey(cfg, "years") ? repeat(collect(_config_years(cfg["years"])), get(cfg, "repeat", 1)) : nothing,
        spinup=_config_spinup(get(cfg, "spinup", nothing), cfg["forcing"], base),
        outputs, backend=backend === nothing ? Symbol(get(cfg, "backend", "threads")) : backend,
        name=get(cfg, "name", splitext(basename(path))[1]))
    dry_run && return exp
    exp.backend == :gpu && !cuda_available() && return Base.invokelatest(_run_with_gpu, path, cfg, set, exp, outdir)
    _keep_config(path, cfg, set, outdir)
    return run!(exp)
end

"""Keep the configuration with the results (with the overrides applied, if any)."""
function _keep_config(path, cfg, set, outdir)
    mkpath(outdir)
    copy_to = joinpath(outdir, basename(path))
    if isempty(set)
        abspath(path) == abspath(copy_to) || cp(path, copy_to; force=true)
    else
        open(io -> (println(io, "# ", basename(path), " with --set ", join(("$k=$v" for (k, v) in set), " "));
            TOML.print(io, cfg)), copy_to, "w")
    end
end

"""Load CUDA (a weak dependency that activates the GPU backend) and run `exp`."""
function _run_with_gpu(path, cfg, set, exp, outdir)
    Core.eval(Main, :(using CUDA))
    Base.invokelatest(cuda_available) || error("The GPU backend is not available (CUDA failed to load or no GPU).")
    _keep_config(path, cfg, set, outdir)
    return Base.invokelatest(run!, exp)
end

"""
Command line: `julia --project -m Chion config.toml [--backend=gpu] [--output-dir=DIR] [--set key=value …]`,
e.g. `--set model.alpha_dry=0.84 --set forcing.surface_fluxes=prescribed`.
"""
function main(args::Vector{String})
    isempty(args) && (println(stderr, "usage: julia --project -m Chion config.toml [--backend=threads|gpu] " *
        "[--output-dir=DIR] [--set section.key=value ...]"); return 1)
    backend = nothing; output_dir = nothing
    files = String[]
    set = Pair{String, Any}[]
    i = 0
    while (i += 1) <= length(args)
        a = args[i]
        if a == "--set" || startswith(a, "--set=")
            kv = a == "--set" ? args[i += 1] : split(a, '='; limit=2)[2]
            k, v = split(kv, '='; limit=2)
            push!(set, String(k) => String(v))
        elseif startswith(a, "--backend=")
            backend = Symbol(split(a, '='; limit=2)[2])
        elseif startswith(a, "--output-dir=")
            output_dir = split(a, '='; limit=2)[2]
        else
            push!(files, a)
        end
    end
    for f in files
        result = run_config(f; backend, output_dir, set)
        println(result)
    end
    return 0
end

@static if isdefined(Base, Symbol("@main"))
    (@main)(args) = main(args)
end
