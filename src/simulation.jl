#=
Simulation: a model with its time integration — the clock, the forcing years, the spin-up and the
outputs. `run!` continues from where the previous call stopped, `time_step!` takes one forcing
step; the time step is set by the forcing calendar.
=#

# ---------------------------------------------------------------------------
# Spin-up
# ---------------------------------------------------------------------------

"""
Abstract supertype of the spin-up that [`run!`](@ref) performs before the first forcing year.

Available subtypes: [`NoSpinup`](@ref), [`CycledSpinup`](@ref), [`ClimatologySpinup`](@ref).
"""
abstract type AbstractSpinup end

"""$(TYPEDEF)

No spin-up: start from the model's state."""
struct NoSpinup <: AbstractSpinup end

"""
$(TYPEDEF)

Run the forcing years `years` in order, `cycles` times, without output.

    CycledSpinup(years; cycles = 5)

# Fields
$(TYPEDFIELDS)
"""
struct CycledSpinup <: AbstractSpinup
    "forcing years of one cycle"
    years::Vector{Int}
    "number of cycles"
    cycles::Int
end
function CycledSpinup(years; cycles::Integer = 5)
    _require(cycles >= 1, "`cycles` must be at least 1.")
    return CycledSpinup(collect(Int, years), Int(cycles))
end

"""
$(TYPEDEF)

Run the years of another forcing (e.g. a day-by-day climatology) `years` times, without output.

    ClimatologySpinup(forcing::Forcing; years = 200)

# Fields
$(TYPEDFIELDS)
"""
struct ClimatologySpinup{F <: Forcing} <: AbstractSpinup
    "the forcing of the spin-up"
    forcing::F
    "number of repetitions of its years"
    years::Int
end
ClimatologySpinup(forcing::Forcing; years::Integer = 200) = ClimatologySpinup(forcing, Int(years))

describe(::NoSpinup) = "none"
describe(s::CycledSpinup) = "$(s.cycles) cycles of the forcing years $(first(s.years))–$(last(s.years))"
describe(s::ClimatologySpinup) = "$(s.years) repetitions of $(s.forcing.description)"

_validate_spinup(::AbstractSpinup, source) = nothing
function _validate_spinup(spinup::CycledSpinup, source)
    missing_years = setdiff(spinup.years, forcing_years(source))
    isempty(missing_years) ||
        throw(ArgumentError("The forcing has no data for the spin-up years $(join(sort(missing_years), ", "))."))
    return nothing
end

# ---------------------------------------------------------------------------
# Debug options
# ---------------------------------------------------------------------------

"""
$(TYPEDEF)

Debug options of a [`Simulation`](@ref).

Set `check_nans = true` to check the prognostic state of the active columns for NaN and Inf after
every launch of the column step (one step when outputs accumulate every step, else up to the next
output period). The run then errors at once, naming the first variable to blow up, the column and
the date. Costs a host copy of the state per launch.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct DebugConfig
    "check the state for NaN and Inf after every launch of the column step (default `false`)"
    check_nans::Bool = false
end

_checked_fields(::BESSIModel) = (:mass, :mass_w, :density, :temperature, :ice_temperature, :Tsrf, :albedo,
    :smb_ice, :runoff, :melt, :refreezing, :vapor_mass)
_checked_fields(model::AbstractSnowModel) = fieldnames(typeof(model.state))

function _check_finite(debug::DebugConfig, model, time)
    debug.check_nans || return nothing
    active = model.cache.active
    for name in _checked_fields(model)
        values = _primal.(Array(getproperty(model.state, name)))
        bad = values isa AbstractMatrix ? findfirst(c -> active[c] && !all(isfinite, view(values, :, c)), axes(values, 2)) :
            findfirst(c -> active[c] && !isfinite(values[c]), eachindex(values))
        bad === nothing ||
            error("Non-finite `$name` in column $bad at $(Dates.format(time, "yyyy-mm-dd HH:MM")) (DebugConfig(check_nans = true)).")
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Clock and simulation
# ---------------------------------------------------------------------------

"""
$(TYPEDEF)

Position of a [`Simulation`](@ref) in its forcing years. Never reset by `run!`, so successive
calls continue where the previous one stopped.

# Fields
$(TYPEDFIELDS)
"""
mutable struct Clock
    "index in `sim.years` of the current year (`length(sim.years) + 1` when finished)"
    year_index::Int
    "next forcing step of the current year"
    step::Int
    "forcing steps taken (spin-up excluded)"
    iteration::Int
    "start of the next forcing step (the end of the last one when finished)"
    time::DateTime
end

"""
$(TYPEDEF)

Run-time state of the outputs: the open recorders and whether the spin-up is done.

# Fields
$(TYPEDFIELDS)
"""
mutable struct IOState
    "one recorder per `Output`, opened at the first forcing step"
    recorders::Vector{OutputRecorder}
    "whether the spin-up has run"
    spun_up::Bool
    "whether the outputs are closed (the simulation finished)"
    closed::Bool
    "wall-clock seconds spent in `run!` and `time_step!`"
    wall_seconds::Float64
end

"""
$(TYPEDEF)

A model with its forcing years, spin-up and outputs. Build with
[`Simulation(model; …)`](@ref Simulation(::AbstractSnowModel)) and advance with [`run!`](@ref) or
[`time_step!`](@ref).

# Fields
$(TYPEDFIELDS)
"""
struct Simulation{M <: AbstractSnowModel, SP <: AbstractSpinup}
    "the model being integrated"
    model::M
    "position in the forcing years"
    clock::Clock
    "forcing years, in the order they are run (they may repeat)"
    years::Vector{Int}
    "spin-up before the first year"
    spinup::SP
    "run directory, field outputs, log, metadata and restarts"
    output::OutputConfig
    "run-time state of the outputs"
    io::IOState
    "name of the run (in logs and output metadata)"
    name::String
    "whether `run!` prints a line per year"
    verbose::Bool
    "debug options"
    debug::DebugConfig
end

"""
    Simulation(model; years = all forcing years, spinup = NoSpinup(), output = OutputConfig(),
               restart = nothing, name = "chion", verbose = true, debug = DebugConfig())

Wrap `model` for a run over the forcing `years` (in order; they may repeat, e.g.
`repeat(1940:1969, 5)`), after `spinup`, writing `output` (an [`OutputConfig`](@ref), or one
[`Output`](@ref) or a list of them without a run directory). The model's state is the initial
state; `restart` (a file from [`write_restart`](@ref)) replaces it and skips the spin-up.

```julia
forcing = Forcing("/data/MAR/MARv3.14.3-10km-daily-ERA5-{year}.nc"; surface_fluxes = :parameterized)
model = BESSIModel(Grid(forcing); forcing, params = Params())
sim = Simulation(model; years = 1940:2026, spinup = CycledSpinup(1940:1969; cycles = 5),
    output = OutputConfig(dir = "output/gris", restart_every = 10,
        outputs = [Output("monthly.nc"; variables = :monthly, frequency = :monthly)]))
run!(sim)
```
"""
function Simulation(model::AbstractSnowModel; years = nothing, spinup::AbstractSpinup = NoSpinup(),
        output = OutputConfig(), restart::Union{Nothing, AbstractString} = nothing,
        name::AbstractString = "chion", verbose::Bool = true, debug::DebugConfig = DebugConfig())
    source = model.forcing.source
    run_years = years === nothing ? forcing_years(source) : collect(Int, years)
    isempty(run_years) && throw(ArgumentError("`years` is empty."))
    missing_years = setdiff(run_years, forcing_years(source))
    isempty(missing_years) || throw(ArgumentError("The forcing has no data for $(join(sort(missing_years), ", "))."))
    _validate_spinup(spinup, source)
    config = _output_config(output)
    restart === nothing || _read_restart!(model, restart, run_years)
    update_forcing!(model, first(run_years))
    clock = Clock(1, 1, 0, first(model.forcing.data.time_values))
    io = IOState(OutputRecorder[], restart !== nothing, false, 0.0)
    return Simulation(model, clock, run_years, spinup, config, io, String(name), verbose, debug)
end

"""Whether all forcing years of `sim` have been run."""
finished(sim::Simulation) = sim.clock.year_index > length(sim.years)

"""
    to_backend(sim, backend)

A simulation of `to_backend(sim.model, backend)` with the same years, spin-up and outputs, before
it has started.
"""
function to_backend(sim::Simulation, backend)
    sim.clock.iteration == 0 && !sim.io.spun_up ||
        throw(ArgumentError("Move a simulation to another backend before it starts (or move its model)."))
    return Simulation(to_backend(sim.model, backend); years = sim.years, spinup = sim.spinup,
        output = sim.output, name = sim.name, verbose = sim.verbose, debug = sim.debug)
end

_log(sim::Simulation, message) = sim.verbose && (println(stdout, "[", sim.name, "] ", message); flush(stdout))

# ---------------------------------------------------------------------------
# Spin-up
# ---------------------------------------------------------------------------

"""
    spin_up!(sim)

Run the spin-up of `sim` (once; `run!` and `time_step!` call it before the first step).
"""
function spin_up!(sim::Simulation)
    sim.io.spun_up && return sim
    _log(sim, "spin-up: $(describe(sim.spinup))")
    _log_line(sim, "spin-up: $(describe(sim.spinup))")
    _spin_up!(sim.model, sim.spinup, sim)
    sim.io.spun_up = true
    update_forcing!(sim.model, sim.years[min(sim.clock.year_index, length(sim.years))])
    return sim
end

_spin_up!(model, ::NoSpinup, sim) = nothing

function _spin_up!(model, spinup::CycledSpinup, sim)
    for cycle in 1:spinup.cycles
        for year in spinup.years
            update_forcing!(model, year)
            _step!(model, 1, forcing_steps(model))
        end
        _log(sim, "spin-up cycle $cycle of $(spinup.cycles) done")
    end
    return nothing
end

function _spin_up!(model::M, spinup::ClimatologySpinup, sim) where {M}
    years = forcing_years(spinup.forcing)
    forcing = AtmosphericForcing(spinup.forcing, first(years), model.params, _backend(model))
    _validate_forcing(model.params, forcing.data)
    # The same state and cache, driven by the spin-up forcing.
    spinup_model = _replace(model; forcing)
    for _ in 1:spinup.years, year in years
        update_forcing!(spinup_model, year)
        _step!(spinup_model, 1, forcing_steps(spinup_model))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Time stepping
# ---------------------------------------------------------------------------

function _open_outputs!(sim::Simulation)
    outputs = sim.output.outputs
    (isempty(outputs) || !isempty(sim.io.recorders)) && return nothing
    model = sim.model
    area = forcing_area(model.forcing.source)
    metadata = _output_metadata(sim)
    try
        for spec in outputs
            mkpath(dirname(spec.path))
            push!(sim.io.recorders, OutputRecorder(spec, model, model.state, model.grid, area, metadata))
        end
    catch
        _close_outputs!(sim, :failed)
        rethrow()
    end
    return nothing
end

function _output_metadata(sim::Simulation)
    return Dict{String, Any}(
        "title" => sim.name,
        "source_model" => "Chion",
        "model" => string(nameof(typeof(sim.model))),
        "forcing" => sim.model.forcing.source.description,
        "spinup" => describe(sim.spinup),
        "years" => "$(first(sim.years))–$(last(sim.years))",
        "created" => string(now()),
    )
end

function _close_outputs!(sim::Simulation, status::Symbol)
    for recorder in sim.io.recorders
        _close_writer!(recorder.writer, status)
    end
    empty!(sim.io.recorders)
    sim.io.closed = true
    return nothing
end

"""
    _advance!(sim, last_step)

Step the current year of `sim` from `sim.clock.step` to `last_step` and record the outputs.
Without outputs that need every step, the steps up to the next end of an output period run in one
launch sequence.
"""
function _advance!(sim::Simulation, last_step::Int)
    model = sim.model
    clock = sim.clock
    year = sim.years[clock.year_index]
    update_forcing!(model, year)
    data = model.forcing.data
    times = data.time_values
    nsteps = length(times)
    recorders = sim.io.recorders
    # Output periods continue into the next run year only when it is the next calendar year.
    next_first = clock.year_index < length(sim.years) && sim.years[clock.year_index + 1] == year + 1 ?
        DateTime(year + 1) : nothing
    closes_at(rec, k) = k < nsteps ?
        _bucket(rec.spec.frequency, times[k], k) != _bucket(rec.spec.frequency, times[k + 1], k + 1) :
        next_first === nothing || _bucket(rec.spec.frequency, times[k], k) != _bucket(rec.spec.frequency, next_first, k + 1)
    every_step = any(rec -> rec.needs_every_step, recorders)
    seconds_per_day = Float64(_primal(model.params.seconds_per_day))
    k = clock.step
    while k <= last_step
        for rec in recorders
            rec.period_start === nothing && _start_period!(rec, model.state, times[k])
        end
        stop = _block_end(recorders, every_step, closes_at, k, last_step)
        _step!(model, k, stop)
        _check_finite(sim.debug, model, times[stop])
        for rec in recorders
            _accumulate_forcing!(rec, data, k:stop, seconds_per_day)
            if every_step
                _accumulate_step!(rec, model.state, Float64(data.dt_days[k]))
            else
                rec.days += sum(Float64(data.dt_days[j]) for j in k:stop)
            end
        end
        closing = [rec for rec in recorders if closes_at(rec, stop)]
        if !isempty(closing)
            host_state = any(rec -> rec.needs_host_state, closing) ? cpu_state(model) : nothing
            for rec in closing
                _close_period!(rec, model.state, host_state, times[stop])
                rec.period_start = nothing
            end
        end
        clock.iteration += stop - k + 1
        k = stop + 1
    end
    clock.step = k
    if k > nsteps
        clock.year_index += 1
        clock.step = 1
        clock.time = times[end] + Millisecond(round(Int, Float64(data.dt_days[end]) * 86_400_000))
    else
        clock.time = times[k]
    end
    return nothing
end

# Last step of the next launch sequence: the whole range without outputs, one step when an output
# accumulates every step, else the next end of an output period.
function _block_end(recorders, every_step::Bool, closes_at, k::Int, last_step::Int)
    isempty(recorders) && return last_step
    every_step && return k
    for j in k:last_step
        any(rec -> closes_at(rec, j), recorders) && return j
    end
    return last_step
end

function _start!(sim::Simulation)
    if _has_dir(sim.output) && sim.clock.iteration == 0 && !sim.io.spun_up
        mkpath(sim.output.dir)
        _write_metadata(sim)
        _log_line(sim, "run \"$(sim.name)\" of a $(nameof(typeof(sim.model))), years $(first(sim.years))–$(last(sim.years))")
    end
    spin_up!(sim)
    sim.io.closed && throw(ArgumentError("The simulation has finished; build a new one to run on."))
    _open_outputs!(sim)
    return nothing
end

function _finish_if_done!(sim::Simulation)
    finished(sim) && !sim.io.closed && _close_outputs!(sim, :complete)
    return nothing
end

"""
    time_step!(sim)

Take one forcing step of `sim` (after the spin-up, on the first call), recording the outputs.
"""
function time_step!(sim::Simulation)
    finished(sim) && throw(ArgumentError("The simulation has finished."))
    t0 = time()
    _start!(sim)
    _advance!(sim, sim.clock.step)
    _finish_if_done!(sim)
    sim.io.wall_seconds += time() - t0
    return sim
end

"""
    run!(sim; years = nothing)

Spin up (on the first call), then run the remaining forcing years of `sim`, or the next `years`
of them, writing the outputs. Successive calls continue where the previous one stopped; the
outputs are closed when the last year is done. Returns `sim`.
"""
function run!(sim::Simulation; years::Union{Nothing, Integer} = nothing)
    t0 = time()
    _start!(sim)
    remaining = length(sim.years) - sim.clock.year_index + 1
    nyears = years === nothing ? remaining : min(Int(years), remaining)
    status = :incomplete
    try
        for _ in 1:nyears
            year = sim.years[sim.clock.year_index]
            t_year = time()
            update_forcing!(sim.model, year)
            _advance!(sim, forcing_steps(sim.model))
            _log(sim, "year $year done ($(sim.clock.year_index - 1) of $(length(sim.years)))")
            _log_year(sim, year, time() - t_year)
            _write_scheduled_restart(sim)
            GC.gc(false)
        end
        status = :complete
    finally
        status === :complete || _close_outputs!(sim, :failed)
    end
    _finish_if_done!(sim)
    update_diagnostics!(sim.model)
    sim.io.wall_seconds += time() - t0
    return sim
end

function Base.show(io::IO, ::MIME"text/plain", sim::Simulation)
    println(io, "Simulation \"", sim.name, "\" of a ", nameof(typeof(sim.model)), "{", number_type(sim.model), "}")
    println(io, "  years:   ", first(sim.years), "–", last(sim.years), " (", length(sim.years), " runs)",
        finished(sim) ? ", finished" : ", next: $(sim.years[sim.clock.year_index]) step $(sim.clock.step)")
    println(io, "  spin-up: ", describe(sim.spinup), sim.io.spun_up ? " (done)" : "")
    _has_dir(sim.output) && println(io, "  run directory: ", sim.output.dir)
    outputs = sim.output.outputs
    print(io, "  outputs: ", isempty(outputs) ? "none" : join((basename(o.path) for o in outputs), ", "))
    return nothing
end

# ---------------------------------------------------------------------------
# Coupling and summaries
# ---------------------------------------------------------------------------

"""
    set_active_mask!(sim, mask; reset_newly_inactive = true)

[`set_active_mask!`](@ref) of the simulation's model.
"""
set_active_mask!(sim::Simulation, mask; kwargs...) = (set_active_mask!(sim.model, mask; kwargs...); sim)

"""
    set_input!(sim, name, values)

[`set_input!`](@ref) of the simulation's model, for the loaded forcing year.
"""
set_input!(sim::Simulation, name::Symbol, values) = (set_input!(sim.model, name, values); sim)

"""
    yearly_step!(sim::Simulation{<:BESSIModel}) -> NamedTuple

Run the next forcing year of `sim` (which must be at the start of a year) and return the coupling
fields of an ice-sheet model: `smb_ice_yearly`, the change of `smb_ice` over the year
(kg m⁻², zero in inactive columns), and `mean_T_srf_K`, the time-mean surface temperature (K, NaN
in inactive columns), as column vectors and, with a spatial grid, as y–x matrices (`*_grid`).
"""
function yearly_step!(sim::Simulation{<:BESSIModel})
    sim.clock.step == 1 || throw(ArgumentError("`yearly_step!` must start at the beginning of a forcing year."))
    finished(sim) && throw(ArgumentError("The simulation has finished."))
    model = sim.model
    _start!(sim)
    year = sim.years[sim.clock.year_index]
    update_forcing!(model, year)
    dt_days = model.forcing.data.dt_days
    smb_before = _host_vector(model.state.smb_ice; copy_array = true)
    Tsrf_sum = fill!(similar(model.state.Tsrf), zero(eltype(model.state.Tsrf)))
    for time_index in 1:forcing_steps(model)
        _advance!(sim, time_index)
        Tsrf_sum .+= model.state.Tsrf .* dt_days[time_index]
    end
    _finish_if_done!(sim)
    smb_ice_yearly = _host_vector(model.state.smb_ice; copy_array = true) .- smb_before
    mean_T_srf_K = _host_vector(Tsrf_sum; copy_array = true) ./ sum(dt_days)
    inactive = .!model.cache.active
    smb_ice_yearly[inactive] .= 0.0
    mean_T_srf_K[inactive] .= NaN
    on_grid(values) = has_spatial_coords(model.grid) ?
        _scatter_to_grid(values, model.grid) : Matrix{Float64}(undef, 0, 0)
    return (; year, smb_ice_yearly, mean_T_srf_K, smb_ice_yearly_grid = on_grid(smb_ice_yearly),
        mean_T_srf_K_grid = on_grid(mean_T_srf_K))
end

function _scatter_to_grid(values::Vector{Float64}, grid::Grid)
    out = fill(NaN, size(grid.mask))
    for (column, value) in enumerate(values)
        out[grid.js[column], grid.is[column]] = value
    end
    return out
end

"""
    smbstats(sim_or_model) -> NamedTuple

Column means (kg m⁻²) of the cumulative mass fluxes since the start of the model, and with column
areas in the forcing their totals (Gt): `smb_ice`, `melt`, `runoff`, `refreezing` (BESSI also
`sublimation`).
"""
smbstats(sim::Simulation) = smbstats(sim.model)

function smbstats(model::AbstractSnowModel)
    names = _smbstats_fields(model)
    area = forcing_area(model.forcing.source)
    active = model.cache.active
    means = NamedTuple{names}(map(name -> _active_mean(getproperty(model.state, name), active), names))
    area === nothing && return (; mean_kg_m2 = means)
    totals = NamedTuple{names}(map(name -> sum(_host_vector(getproperty(model.state, name))[active] .* area[active]) * 1e-6, names))
    return (; mean_kg_m2 = means, total_Gt = totals)
end

_smbstats_fields(::BESSIModel) = (:smb_ice, :melt, :runoff, :refreezing, :sublimation)
_smbstats_fields(::PDDModel) = (:smb_ice, :runoff)
_smbstats_fields(::ITMModel) = (:smb_ice, :melt_cum, :runoff_cum, :refreezing_cum)
_active_mean(values, active) = (host = _host_vector(values); sum(host[active]) / count(active))
