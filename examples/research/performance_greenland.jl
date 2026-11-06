#!/usr/bin/env julia

# Wall time of one simulated year of the full Greenland setup (MAR v3.14.3 10 km daily forcing,
# 18,620 ice-sheet columns, parameterized SEB with dynamic albedo, up to 8 diurnal substeps).
#
# 1. Forcing I/O: load_mar + parameterized_forcing for one year (1989 first to compile the
#    readers, then the timed year, by default 1990).
# 2. Warm-up: one model year from the default initial state (includes compilation).
# 3. Timed: --repetitions one-year runs, each starting from a copy of the warmed-up state, split
#    into setup (host-to-device transfer and workspaces), integration and finalize.
#    --spinup-years=N (default 1) times from the state after N years instead (phase
#    one_year_spun_up), for comparison with a filled column.
#
#   CHION_BACKEND=gpu julia --project=. -O3 --threads 8 examples/research/performance_greenland.jl \
#       --output=output/performance/greenland_raw.csv
# CHION_NUMBER_TYPE=Float32 runs the model in single precision.
ENV["CHION_FORCING_MODE"] = "parameterized"
include(joinpath(@__DIR__, "run_gris_fully_prescribed_transient.jl"))
using Printf

arg(name, default) = (p = "--$name="; i = findfirst(a -> startswith(a, p), ARGS);
    isnothing(i) ? default : split(ARGS[i], "="; limit=2)[2])

synchronize(backend) = backend == :gpu ? CUDA.synchronize() : nothing

function load_year(dir, year)
    path = joinpath(dir, "MARv3.14.3-10km-daily-ERA5-$year.nc")
    seconds = @elapsed begin
        loaded = load_mar(path)
        forcing = parameterized_forcing(loaded.forcing; wind_speed=forcing_wind(path, loaded))
    end
    return loaded, forcing, seconds
end

function one_year(model, forcing, state, backend)
    sim = Simulation(model; forcing, state=deepcopy(state), years=1, backend, write_netcdf=false,
        compute_year_metrics=false)
    GC.gc()
    setup = @elapsed (integrator = init_integrator(sim; io=devnull); synchronize(backend))
    integration = @elapsed (run!(integrator); synchronize(backend))
    finalize = @elapsed result = finalize!(integrator)
    result.status == :complete || error("run ended with status $(result.status)")
    return (; setup, integration, finalize, state=sim.now)
end

function substep_fraction(model, forcing)
    p = model.parameters
    ncol, ntime = size(forcing.air_temperature)
    n = 0
    for t in 1:ntime, c in 1:ncol
        λ = Chion._solar_longitude_deg_from_calendar_day(Chion._calendar_day_of_year(forcing.time_values[t]))
        n += Chion._diurnal_shortwave_substep_count(forcing.dt_days[t], forcing.shortwave_down[c, t],
            forcing.air_temperature[c, t], p.diurnal_shortwave_min_air_temperature, forcing.latitude_deg[c, t],
            λ, p.diurnal_shortwave_threshold, p.diurnal_shortwave_max_substeps) > 1
    end
    return n / (ncol * ntime)
end

function greenland_main()
    backend = Symbol(get(ENV, "CHION_BACKEND", "gpu"))
    dir = get(ENV, "CHION_FORCING_DIR", DOMAIN_SETUP.forcing_dir)
    year = parse(Int, arg("year", "1990"))
    repetitions = parse(Int, arg("repetitions", "3"))
    spinup = parse(Int, arg("spinup-years", "1"))
    output = abspath(arg("output", "greenland_performance.csv"))
    backend == :gpu && !Chion.cuda_available() && error("CUDA is not available")
    device = backend == :gpu ? CUDA.name(CUDA.device()) : strip(Sys.cpu_info()[1].model)
    device = replace(device, "," => " ")

    _, _, compile_io = load_year(dir, year - 1)
    loaded, forcing, io_seconds = load_year(dir, year)
    number_type = Dict("Float64" => Float64, "Float32" => Float32)[get(ENV, "CHION_NUMBER_TYPE", "Float64")]
    model = build_parameterized_model(loaded.grid; number_type)
    ncol, nstep = length(loaded.grid.js), length(forcing.dt_days)
    fraction = substep_fraction(model, forcing)
    println("Greenland $year: $ncol columns, $nstep steps, substep fraction $(round(fraction; digits=3)), " *
        "forcing I/O $(round(io_seconds; digits=1)) s (first load incl. compilation $(round(compile_io; digits=1)) s)")

    warm = one_year(model, forcing, initial_state(model), backend)
    state0 = warm.state
    for _ in 2:spinup
        state0 = one_year(model, forcing, state0, backend).state
    end
    rows = [("warmup_from_initial_state", 0, warm)]
    for repetition in 1:repetitions
        r = one_year(model, forcing, state0, backend)
        push!(rows, (spinup == 1 ? "one_year_from_warm_state" : "one_year_spun_up", repetition, r))
        @printf("%s threads=%d rep=%d: setup %.3f s, integration %.3f s, finalize %.3f s\n",
            backend, Threads.nthreads(), repetition, r.setup, r.integration, r.finalize)
        flush(stdout)
    end
    mkpath(dirname(output))
    isfile(output) || open(io -> print(io, "backend,threads,device,year,columns,steps,diurnal_max_substeps," *
        "substep_fraction,forcing_io_s,phase,repetition,setup_s,integration_s,finalize_s\n"), output, "a")
    open(output, "a") do io
        for (phase, repetition, r) in rows
            println(io, join((backend, Threads.nthreads(), device, year, ncol, nstep,
                model.parameters.diurnal_shortwave_max_substeps, fraction, io_seconds, phase, repetition,
                r.setup, r.integration, r.finalize), ","))
        end
    end
    println("Results appended to $output")
end

greenland_main()
