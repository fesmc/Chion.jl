"""
Run a `Simulation` from a TOML file:

    julia --project -m Chion greenland.toml          # Julia ≥ 1.11
    julia --project -e 'using Chion; run_config("greenland.toml")'

See `examples/configs/` for complete files. Relative output paths are taken relative to
`output_dir` (default: the directory of the TOML file).
"""

using TOML

# ---------------------------------------------------------------------------
# Scheme names: TOML strings are converted once, by dispatch, into scheme instances
# ---------------------------------------------------------------------------

# The abstract type of every scheme field of the parameter structs and schemes.
const _SCHEME_SLOTS = (
    fresh_snow_density_scheme = AbstractFreshSnowDensityScheme,
    albedo_scheme = AbstractAlbedoScheme,
    densification_scheme = AbstractDensificationScheme,
    surface_energy_balance_scheme = AbstractSurfaceEnergyBalanceScheme,
    turbulent_flux_scheme = AbstractTurbulentFluxScheme,
    longwave_scheme = AbstractLongwaveScheme,
    diurnal_cycle = AbstractDiurnalCycle,
    snow_bands = AbstractSemixSnowAlbedo,
    temperature_cycle = AbstractDiurnalTemperatureCycle,
    method = AbstractPDDMethod,
)

# One table of the scheme names: it defines the `_scheme_type(kind, ::Val{name})` methods and the
# names listed in errors.
const _SCHEME_TABLE = (
    (AbstractFreshSnowDensityScheme, "fresh-snow density scheme",
        (constant = ConstantFreshSnowDensity, parameterized = ParameterizedFreshSnowDensity)),
    (AbstractAlbedoScheme, "albedo scheme",
        (constant = ConstantAlbedo, dynamic = DynamicAlbedo, prescribed = PrescribedAlbedo,
         aging = AgingAlbedo, semix = SemixAlbedo)),
    (AbstractDensificationScheme, "densification scheme", (bessi = BESSIDensification, htessel = HTESSELDensification)),
    (AbstractSurfaceEnergyBalanceScheme, "surface-energy-balance scheme",
        (semix = SemixSurfaceEnergyBalance, bessi = BESSISurfaceEnergyBalance)),
    (AbstractTurbulentFluxScheme, "turbulent flux scheme", (semix = SemixTurbulentFluxes, bessi = BESSITurbulentFluxes)),
    (AbstractLongwaveScheme, "longwave scheme", (cloud_proxy = CloudProxyLongwave, graybody = GraybodyLongwave)),
    (AbstractSemixSnowAlbedo, "SEMIX snow albedo", (dang = DangSnowAlbedo, warren_wiscombe = WarrenWiscombeSnowAlbedo)),
    (AbstractDiurnalCycle, "diurnal cycle", (substeps = DiurnalSubsteps, none = NoDiurnalCycle)),
    (AbstractDiurnalTemperatureCycle, "diurnal temperature cycle",
        (sinusoidal = DiurnalTemperatureCycle, none = NoDiurnalTemperatureCycle)),
    (AbstractPDDMethod, "PDD method", (simple = SimplePDD, pism = PISMPDD)),
)

for (kind, description, schemes) in _SCHEME_TABLE
    @eval _scheme_names(::Type{$kind}) = $(keys(schemes))
    @eval _scheme_description(::Type{$kind}) = $description
    for (name, scheme) in pairs(schemes)
        @eval _scheme_type(::Type{$kind}, ::Val{$(QuoteNode(name))}) = $scheme
    end
end

function _scheme_type(::Type{K}, ::Val{name}) where {K, name}
    names = join(("\"" * String(n) * "\"" for n in _scheme_names(K)), ", ")
    throw(ArgumentError("Unknown $(_scheme_description(K)) \"$name\". Use $names."))
end

"""
    _scheme_from_config(kind, value)

The scheme of abstract type `kind` that a TOML value describes: a name (`"semix"`) for the scheme
with its default coefficients, or a table with the name under `type` and coefficients
(`{type = "aging", cold_timescale_days = 15}`), nested schemes included.
"""
_scheme_from_config(::Type{K}, name::AbstractString) where {K} = _scheme_type(K, Val(Symbol(name)))()
function _scheme_from_config(::Type{K}, table::AbstractDict) where {K}
    haskey(table, "type") || throw(ArgumentError("A $(_scheme_description(K)) table needs a `type`."))
    T = _scheme_type(K, Val(Symbol(table["type"])))
    kwargs = [Symbol(k) => _config_value(Symbol(k), v) for (k, v) in table if k != "type"]
    isempty(kwargs) && return T()
    fieldcount(T) == 0 && throw(ArgumentError("$(nameof(T)) has no coefficients."))
    return T(; kwargs...)
end

# A TOML value of a parameter or scheme field: schemes by name or table, arrays as tuples.
_config_value(name::Symbol, value) =
    haskey(_SCHEME_SLOTS, name) ? _scheme_from_config(_SCHEME_SLOTS[name], value) :
    value isa AbstractVector ? Tuple(value) : value

"""
    _scheme_name(scheme)

The TOML name of `scheme` (e.g. `:semix` for `SemixSurfaceEnergyBalance()`), or its type name.
"""
function _scheme_name(scheme)
    for (kind, _, schemes) in _SCHEME_TABLE
        scheme isa kind || continue
        for (name, T) in pairs(schemes)
            scheme isa T && return name
        end
    end
    return nameof(typeof(scheme))
end

# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------

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

const _NAMED_PARAMS = (grisv1 = ParamsGrISv1, aisv1 = ParamsAISv1)
const _MODEL_TYPES = (bessi = (BESSIModel, Params), pdd = (PDDModel, PDDParams), itm = (ITMModel, ITMParams))

_config_number_type(name) = name == "Float32" ? Float32 : name == "Float64" ? Float64 :
    throw(ArgumentError("`number_type` must be \"Float64\" or \"Float32\"."))

"""The model type and its parameters from the `[model]` table."""
function _config_params(d)
    mtype = Symbol(lowercase(get(d, "type", "bessi")))
    haskey(_MODEL_TYPES, mtype) || throw(ArgumentError("Unknown model type \"$mtype\"; use bessi, pdd or itm."))
    M, P = _MODEL_TYPES[mtype]
    NF = _config_number_type(get(d, "number_type", "Float64"))
    kwargs = [Symbol(k) => _config_value(Symbol(k), v) for (k, v) in d if !(k in ("type", "number_type", "setup"))]
    setup = Symbol(lowercase(get(d, "setup", "grisv1")))
    if P === Params
        haskey(_NAMED_PARAMS, setup) || throw(ArgumentError("Unknown setup \"$setup\"; use GrISv1 or AISv1."))
        return M, _NAMED_PARAMS[setup](; NF, kwargs...)
    end
    haskey(d, "setup") && throw(ArgumentError("`setup` applies to BESSI models only."))
    return M, P(; NF, kwargs...)
end

"""The layering from the `[grid]` table."""
_config_layering(d) = Layering(; (Symbol(k) => _config_value(Symbol(k), v) for (k, v) in d)...)

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
    run_config(path; backend = nothing, output_dir = nothing, set = (), dry_run = false) -> Simulation

Build the [`Simulation`](@ref) described by the TOML file `path` and run it. `backend` (`:threads`
or `:gpu`) and `output_dir` override the values in the file; `set` overrides any entry, e.g.
`set = ["model.alpha_dry" => 0.84, "model.diurnal_cycle.max_substeps" => 4]`. `dry_run = true`
returns the simulation without running it.

The file has the run settings (`name`, `backend`, `years`, `repeat`, `output_dir`) and the tables
`[forcing]`, `[grid]` (the [`Layering`](@ref)), `[model]` (`type`, `setup`, `number_type` and any
field of the parameters; a scheme as a name or a table with `type`), `[spinup]` and `[[output]]`;
see `examples/configs/`.
"""
function run_config(path::AbstractString; backend=nothing, output_dir=nothing, set=(), dry_run::Bool=false)
    cfg = _apply_overrides!(TOML.parsefile(path), set)
    base = dirname(abspath(path))
    outdir = output_dir === nothing ? get(cfg, "output_dir", base) : output_dir
    isabspath(outdir) || (outdir = joinpath(output_dir === nothing ? base : pwd(), outdir))
    backend_name = Symbol(backend === nothing ? get(cfg, "backend", "threads") : backend)
    backend_name === :gpu && !cuda_available() && Core.eval(Main, :(using CUDA))
    return Base.invokelatest(_run_config, path, cfg, set, base, outdir, backend_name, dry_run)
end

function _run_config(path, cfg, set, base, outdir, backend_name, dry_run)
    forcing = _config_forcing(cfg["forcing"], base)
    M, params = _config_params(get(cfg, "model", Dict{String, Any}()))
    grid = Grid(forcing; NF = number_type(params), layering = _config_layering(get(cfg, "grid", Dict{String, Any}())))
    model = M(grid; forcing, params, backend = _ka_backend_from(backend_name))
    outputs = [_config_output(o, outdir) for o in get(cfg, "output", Any[])]
    sim = Simulation(model;
        # `repeat = n` runs the years n times (e.g. to record every cycle of an equilibrium run)
        years = haskey(cfg, "years") ? repeat(collect(_config_years(cfg["years"])), get(cfg, "repeat", 1)) : nothing,
        spinup = _config_spinup(get(cfg, "spinup", nothing), cfg["forcing"], base),
        outputs, name = get(cfg, "name", splitext(basename(path))[1]))
    dry_run && return sim
    _keep_config(path, cfg, set, outdir)
    return run!(sim)
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

"""
Command line: `julia --project -m Chion config.toml [--backend=gpu] [--output-dir=DIR] [--set key=value …]`,
e.g. `--set model.alpha_dry=0.84 --set model.albedo_scheme=semix --set forcing.surface_fluxes=prescribed`.
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
        sim = run_config(f; backend, output_dir, set)
        println(sim.name, ": ", length(sim.years), " years in ", round(sim.io.wall_seconds; digits=1), " s")
    end
    return 0
end

@static if isdefined(Base, Symbol("@main"))
    (@main)(args) = main(args)
end
