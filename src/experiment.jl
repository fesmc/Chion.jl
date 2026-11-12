"""
`Experiment`: a model, its forcing, a spin-up, and the outputs, run with `run!`.

    forcing = Forcing("/data/MAR/MARv3.14.3-10km-daily-ERA5-{year}.nc"; surface_fluxes=:parameterized)
    model   = BESSIModel(forcing)
    exp = Experiment(model, forcing; years=1940:2026, spinup=CycledSpinup(1940:1969; cycles=5),
        outputs=[Output("smb_monthly.nc"; frequency=:monthly),
                 Output("budget.csv"; variables=:budget, frequency=:yearly, aggregate=:area_sum)],
        backend=:gpu)
    result = run!(exp)          # result.state is the final state
"""

# ---------------------------------------------------------------------------
# Spin-up
# ---------------------------------------------------------------------------

abstract type AbstractSpinup end

"""No spin-up: start from the initial state of the model (or the `state` of the `Experiment`)."""
struct NoSpinup <: AbstractSpinup end

"""
    CycledSpinup(years; cycles=5)

Run the real forcing years `years` in order, `cycles` times, from the initial state.
"""
struct CycledSpinup <: AbstractSpinup
    years::Vector{Int}
    cycles::Int
end
CycledSpinup(years; cycles::Integer=5) = (cycles >= 1 || error("`cycles` must be at least 1.");
    CycledSpinup(collect(Int, years), Int(cycles)))

"""
    ClimatologySpinup(forcing::Forcing; years=200)

Repeat one year of `forcing` (e.g. a day-by-day climatology) `years` times.
"""
struct ClimatologySpinup <: AbstractSpinup
    forcing::Forcing
    years::Int
end
ClimatologySpinup(forcing::Forcing; years::Integer=200) = ClimatologySpinup(forcing, Int(years))

describe(::NoSpinup) = "none"
describe(s::CycledSpinup) = "$(s.cycles) cycles of the forcing years $(first(s.years))–$(last(s.years))"
describe(s::ClimatologySpinup) = "$(s.years) years of $(s.forcing.description)"

# ---------------------------------------------------------------------------
# Experiment
# ---------------------------------------------------------------------------

struct Experiment{M <: AbstractSnowModel}
    model::M
    forcing::Forcing
    years::Vector{Int}
    spinup::AbstractSpinup
    outputs::Vector{Output}
    backend::Symbol
    state::Any
    name::String
    verbose::Bool
end

"""
    Experiment(model, forcing; years=all, spinup=NoSpinup(), outputs=Output[], backend=:threads,
               state=nothing, name="chion", verbose=true)

`years` may repeat years (e.g. `repeat(1940:1969, 5)` to record every cycle). `state` continues
from a given state instead of the initial state.
"""
function Experiment(model::AbstractSnowModel, forcing::Forcing; years=nothing, spinup::AbstractSpinup=NoSpinup(),
        outputs=Output[], backend=:threads, state=nothing, name::AbstractString="chion", verbose::Bool=true)
    ys = years === nothing ? forcing_years(forcing) : collect(Int, years)
    available = Set(forcing_years(forcing))
    missing_years = setdiff(Set(ys), available)
    isempty(missing_years) || error("The forcing has no data for $(join(sort(collect(missing_years)), ", ")).")
    if spinup isa CycledSpinup
        m = setdiff(Set(spinup.years), available)
        isempty(m) || error("The forcing has no data for the spin-up years $(join(sort(collect(m)), ", ")).")
    end
    ncols(model.grid) == ncols(forcing_grid(forcing)) || error("The model and the forcing have different numbers of columns.")
    out = outputs isa Output ? [outputs] : collect(Output, outputs)
    paths = [o.path for o in out]
    length(unique(paths)) == length(paths) || error("Two outputs write to the same file.")
    return Experiment(model, forcing, ys, spinup, out, Symbol(backend), state, String(name), verbose)
end

"""`BESSIModel` on the columns of `forcing`."""
BESSIModel(forcing::Forcing; kwargs...) = BESSIModel(forcing_grid(forcing); kwargs...)
PDDModel(forcing::Forcing; kwargs...) = PDDModel(forcing_grid(forcing); kwargs...)
ITMModel(forcing::Forcing; kwargs...) = ITMModel(forcing_grid(forcing); kwargs...)

struct ExperimentResult
    state::Any
    outputs::Vector{String}
    years::Vector{Int}
    wall_seconds::Float64
end

Base.show(io::IO, r::ExperimentResult) =
    print(io, "ExperimentResult($(length(r.years)) years, $(round(r.wall_seconds; digits=1)) s, outputs: ",
        isempty(r.outputs) ? "none" : join(r.outputs, ", "), ")")

_log(exp::Experiment, msg) = exp.verbose && (println(stdout, "[$(exp.name)] ", msg); flush(stdout))

"""Run one year of forcing from `state` without output; returns the new host state."""
function _run_year(model, forcing::SnowpackForcing, state, backend, name)
    sim = Simulation(model; forcing, state, years=1, backend, write_netcdf=false, compute_year_metrics=false, name)
    integ = init_integrator(sim; io=devnull)
    while !finished(integ)
        _step_scheduled_range!(integ, integ.time_index:length(forcing.time_values))
    end
    finalize!(integ)
    return sim.now
end

function _spin_up(exp::Experiment, state)
    s = exp.spinup
    s isa NoSpinup && return state
    if s isa ClimatologySpinup
        # The whole climatology file (or its first year when it holds one file per year).
        f = s.forcing.source isa WholeSource ? _whole(s.forcing).forcing : load_forcing_year(s.forcing, first(forcing_years(s.forcing)))
        s.forcing.source isa WholeSource && s.forcing.transform !== nothing && (f = s.forcing.transform(f, year(first(f.time_values))))
        sim = Simulation(exp.model; forcing=f, state, years=s.years, backend=exp.backend, write_netcdf=false,
            compute_year_metrics=false, name="$(exp.name)_spinup")
        run!(sim; io=devnull).status == :complete || error("Spin-up failed.")
        return sim.now
    end
    for c in 1:s.cycles
        for y in s.years
            state = _run_year(exp.model, load_forcing_year(exp.forcing, y), state, exp.backend, "$(exp.name)_spinup_$(c)_$(y)")
            GC.gc(false)
        end
        _log(exp, "spin-up cycle $c of $(s.cycles) done")
    end
    return state
end

function _netcdf_metadata(exp::Experiment)
    return Dict{String, Any}(
        "title" => exp.name,
        "source_model" => "Chion",
        "model" => string(nameof(typeof(exp.model))),
        "forcing" => exp.forcing.description,
        "spinup" => describe(exp.spinup),
        "years" => "$(first(exp.years))–$(last(exp.years))",
        "created" => string(now()),
    )
end

"""
    run!(experiment) -> ExperimentResult

Spin up, then run every year of `experiment.years` and write the outputs.
"""
function run!(exp::Experiment)
    t0 = time()
    state = exp.state === nothing ? initial_state(exp.model) : exp.state
    _log(exp, "spin-up: $(describe(exp.spinup))")
    state = _spin_up(exp, state)
    meta = _netcdf_metadata(exp)
    area = forcing_area(exp.forcing)
    recorders = OutputRecorder[]
    status = :failed
    try
        for spec in exp.outputs
            push!(recorders, OutputRecorder(spec, exp.model, state, exp.model.grid, area, meta))
        end
        nyears = length(exp.years)
        for (iy, y) in enumerate(exp.years)
            forcing = load_forcing_year(exp.forcing, y)
            next_first = iy < nyears && exp.years[iy + 1] == y + 1 ? DateTime(y + 1) : nothing
            state = _run_year_with_outputs(exp, forcing, state, recorders, next_first, "$(exp.name)_$(y)")
            _log(exp, "year $y done ($iy of $nyears)")
            GC.gc(false)
        end
        status = :complete
    finally
        for rec in recorders
            _close_writer!(rec.writer, status)
        end
    end
    return ExperimentResult(state, [o.path for o in exp.outputs], exp.years, time() - t0)
end

"""
One forcing year with output recording. A period closes at a time step if the next step falls in
another calendar bucket; at the end of the year it closes unless the next year continues it.
"""
function _run_year_with_outputs(exp, forcing::SnowpackForcing, state, recorders, next_first, name)
    isempty(recorders) && return _run_year(exp.model, forcing, state, exp.backend, name)
    sim = Simulation(exp.model; forcing, state, years=1, backend=exp.backend, write_netcdf=false,
        compute_year_metrics=false, name)
    integ = init_integrator(sim; io=devnull)
    runtime = integ.model_runtime.data
    dev = runtime.state
    times = forcing.time_values
    nsteps = length(times)
    seconds_per_day = hasproperty(dev, :c) ? Float64(dev.c.seconds_per_day) : 86_400.0   # PDD and ITM states have no constants
    every_step = any(r -> r.needs_every_step, recorders)
    for rec in recorders
        rec.period_start === nothing && _start_period!(rec, dev, first(times))
    end
    closes_at(rec, k) = k < nsteps ? _bucket(rec.spec.frequency, times[k], k) != _bucket(rec.spec.frequency, times[k + 1], k + 1) :
        next_first === nothing || _bucket(rec.spec.frequency, times[k], k) != _bucket(rec.spec.frequency, next_first, k + 1)
    k = 1
    while k <= nsteps
        stop = every_step ? k : something(findfirst(j -> any(rec -> closes_at(rec, j), recorders), k:nsteps), nsteps - k + 1) + k - 1
        _step_scheduled_range!(integ, k:stop)
        for rec in recorders
            _accumulate_forcing!(rec, forcing, k:stop, seconds_per_day)
            if every_step
                _accumulate_step!(rec, dev, Float64(forcing.dt_days[k]))
            else
                rec.days += sum(Float64(forcing.dt_days[j]) for j in k:stop)
            end
        end
        closing = [rec for rec in recorders if closes_at(rec, stop)]
        if !isempty(closing)
            host = any(r -> r.needs_host_state, closing) ? (_sync_bessi_state!(sim.now, runtime, integ.timings); sim.now) : nothing
            for rec in closing
                _close_period!(rec, dev, host, times[stop])
                rec.period_start = nothing
                stop < nsteps && _start_period!(rec, dev, times[stop + 1])
            end
        end
        k = stop + 1
    end
    finalize!(integ)
    # Periods continuing into the next year restart from the synced host state's counters, which are
    # identical to the backend values.
    return sim.now
end
