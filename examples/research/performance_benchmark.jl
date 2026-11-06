#!/usr/bin/env julia

# Throughput benchmark of the default BESSIModel (fine near-surface layers, 5-layer ice substrate,
# SEMIX, dynamic albedo, up to 8 diurnal substeps) on synthetic daily forcing.
#
# The forcing is a seasonal Arctic cycle (365 daily steps per year) that exercises melt,
# percolation, refreezing and the diurnal substeps: per-column air temperature (full ncol x 365
# matrix; annual means from -19 to -7 degC, 13 K seasonal amplitude, 2 K synoptic variability with a
# column-dependent phase), per-column latitude 65-75 N, a daily-mean shortwave cycle for 70 N,
# 1.2 mm/day snowfall and 0.5 mm/day summer rain. Each case is first run for one year (probe);
# the number of years per repetition is then chosen so that one repetition takes about
# --target-seconds (fewer repetitions if they would exceed --case-budget-seconds; larger column
# counts are skipped once a probe year exceeds --max-probe-seconds). Compilation is excluded by a warm-up run. Times are split into setup
# (state/forcing transfer and workspaces), integration (the time loop) and finalize.
#
#   julia --project=. -O3 --threads 16 examples/research/performance_benchmark.jl \
#       --backend=threads --columns=100,1000,10000 --output=output/performance/scaling_raw.csv
# --number-type=Float32 runs the model in single precision (written to a separate CSV by the launchers).
using CUDA  # loads the GPU backend (ChionCUDAExt)
using Chion
using Dates: dayofyear
using Printf
using Statistics: mean

arg(name, default) = (p = "--$name="; i = findfirst(a -> startswith(a, p), ARGS);
    isnothing(i) ? default : split(ARGS[i], "="; limit=2)[2])
list(name, default) = parse.(Int, split(arg(name, default), ","))

const STEPS = 365
const DOY = collect(1:STEPS)

seasonal(t) = -cospi(2 * (t - 15) / 365)                  # -1 mid-January, +1 mid-July
declination(t) = 23.44 * sind(360 * (284 + t) / 365)
function toa_daily_mean(lat, t)
    δ = declination(t)
    c = -tand(lat) * tand(δ)
    h0 = c >= 1 ? 0.0 : c <= -1 ? π : acos(c)
    return 1361 / π * (h0 * sind(lat) * sind(δ) + cosd(lat) * cosd(δ) * sin(h0))
end
const SHORTWAVE = [0.55 * toa_daily_mean(70.0, t) for t in DOY]   # W m-2, daily mean at the surface
const SNOWFALL = fill(1.2 / 86_400, STEPS)                       # kg m-2 s-1
const RAINFALL = [152 <= t <= 243 ? 0.5 / 86_400 : 0.0 for t in DOY]

annual_mean_c(f) = -19 + 12 * f            # f in [0, 1] across columns
column_fraction(c, ncol) = ncol == 1 ? 0.5 : (c - 1) / (ncol - 1)

function synthetic_forcing(ncol)
    air = Matrix{Float64}(undef, ncol, STEPS)
    Threads.@threads for c in 1:ncol
        f = column_fraction(c, ncol)
        phase = 2π * ((c * 0.618034) % 1)
        for t in 1:STEPS
            air[c, t] = 273.15 + annual_mean_c(f) + 13 * seasonal(t) + 2 * sin(2π * t / 9.3 + phase)
        end
    end
    latitude = [65 + 10 * ((c * 0.381966) % 1) for c in 1:ncol]
    forcing = SnowpackForcing(dt_days=ones(STEPS), ncol=ncol, air_temperature=air, snowfall_rate=SNOWFALL,
        rainfall_rate=RAINFALL, shortwave_down=SHORTWAVE, latitude_deg=latitude)
    return forcing
end

"""Fraction of column-days that the default model splits into diurnal substeps."""
function substep_fraction(model, forcing)
    p = model.parameters
    ncol, ntime = size(forcing.air_temperature)
    sample = unique(round.(Int, range(1, ncol; length=min(ncol, 2000))))
    n = 0
    for t in 1:ntime, c in sample
        λ = Chion._solar_longitude_deg_from_calendar_day(dayofyear(forcing.time_values[t]))
        n += Chion._diurnal_shortwave_substep_count(forcing.dt_days[t], forcing.shortwave_down[c, t],
            forcing.air_temperature[c, t], p.diurnal_shortwave_min_air_temperature, forcing.latitude_deg[c, t],
            λ, p.diurnal_shortwave_threshold, p.diurnal_shortwave_max_substeps) > 1
    end
    return n / (length(sample) * ntime)
end

synchronize(backend) = backend == :gpu ? CUDA.synchronize() : nothing

function timed_run(ncol, ntot, years, backend, forcing; number_type=Float64)
    model = BESSIModel(SnowpackGrid(ncol); Ntot=ntot, number_type)
    sim = Simulation(model; forcing, years, backend, write_netcdf=false, history_year_stride=0,
        compute_year_metrics=false)
    GC.gc()
    setup = @elapsed begin
        integrator = init_integrator(sim; io=devnull)
        synchronize(backend)
    end
    integration = @elapsed begin
        run!(integrator)
        synchronize(backend)
    end
    finalize = @elapsed result = finalize!(integrator)
    result.status == :complete || error("run ended with status $(result.status)")
    state = sim.now
    perm(x) = mean(Array(x)) / years
    diagnostics = (melt=perm(state.melt), refreezing=perm(state.refreezing), runoff=perm(state.runoff),
        ice_exposed=mean(Array(state.thickness) .< 0.01))
    return (; setup, integration, finalize, diagnostics, model)
end

const HEADER = "backend,threads,device,time_block_steps,columns,ntot,years,repetition,setup_s,integration_s," *
    "finalize_s,seconds_per_year,column_years_per_s,mean_melt_kg_m2_yr,mean_refreezing_kg_m2_yr," *
    "mean_runoff_kg_m2_yr,ice_exposed_fraction,substep_fraction"

function main()
    backend = Symbol(arg("backend", "threads"))
    columns = list("columns", "100,1000,10000")
    ntots = list("ntots", string(Chion.DEFAULT_NTOT))
    repetitions = parse(Int, arg("repetitions", "3"))
    target = parse(Float64, arg("target-seconds", "10"))
    max_years = parse(Int, arg("max-years", "100"))
    max_probe = parse(Float64, arg("max-probe-seconds", "60"))
    budget = parse(Float64, arg("case-budget-seconds", "240"))
    output = abspath(arg("output", "performance_scaling.csv"))
    number_type = Dict("Float64" => Float64, "Float32" => Float32)[arg("number-type", "Float64")]
    backend == :gpu && !Chion.cuda_available() && error("CUDA is not available")
    device = backend == :gpu ? replace(CUDA.name(CUDA.device()), "," => " ") :
        replace(strip(Sys.cpu_info()[1].model), "," => " ")
    block = get(ENV, backend == :gpu ? "CHION_GPU_TIME_BLOCK_STEPS" : "CHION_CPU_TIME_BLOCK_STEPS", "1")
    mkpath(dirname(output))
    isfile(output) || open(io -> println(io, HEADER), output, "a")

    warm = synthetic_forcing(64)
    for ntot in ntots
        timed_run(64, ntot, 2, backend, warm; number_type)
    end
    println("warm-up done: backend=$backend threads=$(Threads.nthreads()) device=$device block=$block number_type=$number_type")
    fraction_cache = Dict{Int, Float64}()

    for ntot in ntots
        for ncol in columns
            forcing = synthetic_forcing(ncol)
            fraction = get!(fraction_cache, ncol) do
                substep_fraction(BESSIModel(SnowpackGrid(ncol)), forcing)
            end
            probe = timed_run(ncol, ntot, 1, backend, forcing; number_type)
            years = clamp(ceil(Int, target / probe.integration), 1, max_years)
            reps = clamp(floor(Int, budget / (years * probe.integration)), 1, repetitions)
            for repetition in 0:reps
                r = repetition == 0 ? probe : timed_run(ncol, ntot, years, backend, forcing; number_type)
                y = repetition == 0 ? 1 : years
                d = r.diagnostics
                row = (backend, Threads.nthreads(), device, block, ncol, ntot, y, repetition, r.setup, r.integration,
                    r.finalize, r.integration / y, ncol * y / r.integration, d.melt, d.refreezing, d.runoff,
                    d.ice_exposed, fraction)
                open(io -> println(io, join(row, ",")), output, "a")
                @printf("%s threads=%d block=%s ncol=%d ntot=%d rep=%d years=%d: %.3f s/yr, %.3g col-yr/s (melt %.0f, refr %.0f kg/m2/yr, substeps %.2f)\n",
                    backend, Threads.nthreads(), block, ncol, ntot, repetition, y, r.integration / y,
                    ncol * y / r.integration, d.melt, d.refreezing, fraction)
                flush(stdout)
            end
            forcing = nothing
            GC.gc()
            if probe.integration > max_probe
                println("probe year took $(round(probe.integration; digits=1)) s > $max_probe s; skipping larger column counts")
                break
            end
        end
    end
    println("Results appended to $output")
end

main()
