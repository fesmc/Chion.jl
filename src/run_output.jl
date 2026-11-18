#=
Run-level output of a Simulation, as Laddie's `OutputConfig`: one run directory with the field
outputs, a log, the run metadata (a TOML file that `run_config` reruns) and restart files.
=#

using Serialization: serialize, deserialize

"""
$(TYPEDEF)

Where and what a [`Simulation`](@ref) writes.

    OutputConfig(; dir = nothing, outputs = Output[], log = true, metadata = true, restart_every = 0)

With a run directory `dir`, the relative paths of the [`Output`](@ref)s are taken relative to it,
and the simulation writes there:
- `log.txt`: one line per forcing year (column-mean SMB, melt, runoff; wall time),
- `chion_run.toml`: the grid, the parameters and the run settings in the format of
  [`run_config`](@ref) (rerunning it with the forcing section added reproduces the run),
- `restart_<year>.jls` every `restart_every` years and at the end of the run (`0`: only at the end),
  read by `Simulation(model; restart = path)`.

Without `dir` only the outputs are written, at their paths.

# Fields
$(TYPEDFIELDS)
"""
struct OutputConfig
    "run directory, or `nothing`"
    dir::Union{Nothing, String}
    "field outputs, one file each, with paths resolved against `dir`"
    outputs::Vector{Output}
    "whether to write `log.txt`"
    log::Bool
    "whether to write `chion_run.toml`"
    metadata::Bool
    "years between restart files (`0`: only at the end); restarts need `dir`"
    restart_every::Int
end

function OutputConfig(; dir::Union{Nothing, AbstractString} = nothing, outputs = Output[], log::Bool = true,
        metadata::Bool = true, restart_every::Integer = 0)
    _require(restart_every >= 0, "`restart_every` must be non-negative.")
    run_dir = dir === nothing ? nothing : abspath(String(dir))
    specs = [_with_path(o, run_dir === nothing || isabspath(o.path) ? o.path : joinpath(run_dir, o.path))
             for o in (outputs isa Output ? [outputs] : collect(Output, outputs))]
    paths = [o.path for o in specs]
    allunique(paths) || throw(ArgumentError("Two outputs write to the same file."))
    return OutputConfig(run_dir, specs, log, metadata, Int(restart_every))
end

"""An `OutputConfig` from an `Output` or a list of them (no run directory)."""
_output_config(output::OutputConfig) = output
_output_config(outputs) = OutputConfig(; outputs)

_with_path(o::Output, path) = Output(path, (getfield(o, name) for name in fieldnames(Output)[2:end])...)

_has_dir(output::OutputConfig) = output.dir !== nothing

# ---------------------------------------------------------------------------
# Log
# ---------------------------------------------------------------------------

function _log_line(sim, line)
    output = sim.output
    (_has_dir(output) && output.log) || return nothing
    open(joinpath(output.dir, "log.txt"), "a") do io
        println(io, Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "  ", line)
    end
    return nothing
end

function _log_year(sim, year, wall_seconds)
    (_has_dir(sim.output) && sim.output.log) || return nothing
    stats = smbstats(sim.model).mean_kg_m2
    values = join(("$(name) $(round(value; sigdigits = 6))" for (name, value) in pairs(stats)), ", ")
    _log_line(sim, "year $year done ($(sim.clock.year_index - 1) of $(length(sim.years))), " *
        "$(round(wall_seconds; digits = 2)) s; cumulative column means (kg m⁻²): $values")
    return nothing
end

# ---------------------------------------------------------------------------
# Run metadata (in the format of run_config)
# ---------------------------------------------------------------------------

# A parameter value as TOML: numbers as they are, schemes as tables with their name under `type`.
_toml_value(value::Integer) = value
_toml_value(value::AbstractFloat) = Float64(value)
_toml_value(value::Tuple) = collect(value)
function _toml_value(scheme)
    table = Dict{String, Any}("type" => String(_scheme_name(scheme)))
    for name in fieldnames(typeof(scheme))
        table[String(name)] = _toml_value(getfield(scheme, name))
    end
    return table
end

_model_type_name(::BESSIModel) = "bessi"
_model_type_name(::PDDModel) = "pdd"
_model_type_name(::ITMModel) = "itm"

function _git_revision()
    try
        return readchomp(pipeline(`git -C $(pkgdir(@__MODULE__)) describe --always --dirty`; stderr = devnull))
    catch
        return "unknown"
    end
end

"""The TOML dictionary of a simulation: run settings, `[grid]` and `[model]` as `run_config` reads them."""
function _run_metadata(sim)
    model = sim.model
    params = model.params
    model_table = Dict{String, Any}("type" => _model_type_name(model), "number_type" => string(number_type(params)))
    for name in fieldnames(typeof(params))
        model_table[String(name)] = _toml_value(getfield(params, name))
    end
    layering = model.grid.layering
    grid_table = Dict{String, Any}(String(name) => _toml_value(getfield(layering, name)) for name in fieldnames(Layering))
    return Dict{String, Any}(
        "name" => sim.name,
        "years" => sim.years,
        "chion_version" => string(pkgversion(@__MODULE__)),
        "chion_revision" => _git_revision(),
        "created" => string(now()),
        "forcing_description" => model.forcing.source.description,
        "spinup_description" => describe(sim.spinup),
        "columns" => ncols(model),
        "grid" => grid_table,
        "boundary" => _boundary_metadata(model),
        "model" => model_table,
    )
end

_boundary_metadata(model::BESSIModel) =
    Dict{String, Any}("surface" => _toml_value(model.boundary.surface), "bottom" => _toml_value(model.boundary.bottom))
_boundary_metadata(model) = Dict{String, Any}()

function _write_metadata(sim)
    (_has_dir(sim.output) && sim.output.metadata) || return nothing
    open(joinpath(sim.output.dir, "chion_run.toml"), "w") do io
        println(io, "# Chion run \"", sim.name, "\": grid and parameters in the format of run_config.")
        println(io, "# Add the [forcing] section (and [spinup], [[output]]) to rerun it.")
        TOML.print(io, _run_metadata(sim); sorted = true)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Restarts
# ---------------------------------------------------------------------------

"""
    write_restart(sim, path)

Write the host state of `sim`'s model and the next forcing year to `path` (Julia `Serialization`;
read it with `Simulation(model; restart = path)` in the same Chion version).
"""
function write_restart(sim, path::AbstractString)
    sim.clock.step == 1 || throw(ArgumentError("Restarts are written at the start of a forcing year."))
    next_year = finished(sim) ? nothing : sim.years[sim.clock.year_index]
    serialize(path, (; version = 1, state = cpu_state(sim.model), next_year, last_year = _last_run_year(sim),
        model = string(nameof(typeof(sim.model))), number_type = string(number_type(sim.model))))
    return path
end

_last_run_year(sim) = sim.clock.year_index > 1 ? sim.years[sim.clock.year_index - 1] : nothing

function _write_scheduled_restart(sim)
    output = sim.output
    _has_dir(output) || return nothing
    year = _last_run_year(sim)
    year === nothing && return nothing
    due = finished(sim) || (output.restart_every > 0 && (sim.clock.year_index - 1) % output.restart_every == 0)
    due || return nothing
    path = write_restart(sim, joinpath(output.dir, "restart_$(year).jls"))
    _log_line(sim, "restart written: $(basename(path))")
    return nothing
end

"""Load a restart file into the model of `sim`; the spin-up is then skipped."""
function _read_restart!(model, path::AbstractString, years)
    data = deserialize(path)
    data.model == string(nameof(typeof(model))) ||
        throw(ArgumentError("The restart is of a $(data.model), the model a $(nameof(typeof(model)))."))
    set_state!(model, data.state)
    data.next_year === nothing || data.next_year == first(years) ||
        @warn "The restart continues at $(data.next_year), the simulation starts at $(first(years))."
    return nothing
end
