#!/usr/bin/env julia

using Chion

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SCRIPT = abspath(@__FILE__)

function arg(name, default)
    prefix = "--$name="
    for value in ARGS
        startswith(value, prefix) && return split(value, "="; limit=2)[2]
    end
    return default
end
list(name, default) = parse.(Int, split(arg(name, default), ","))

function simulation(columns, ntot, years, steps, backend)
    forcing = SnowpackForcing(
        dt_days=fill(365.0 / steps, steps),
        ncol=columns,
        air_temperature_c=-5.0,
        snowfall_mm_day=1.0,
        rainfall_mm_day=0.0,
        shortwave_down=100.0,
    )
    return Simulation(
        BESSIModel(SnowpackGrid(columns); Ntot=ntot);
        forcing,
        years,
        backend,
        write_netcdf=false,
        history_year_stride=0,
        compute_year_metrics=false,
    )
end

function worker()
    backend = Symbol(arg("backend", "threads"))
    columns = list("columns", "128,1024,10000")
    ntots = list("ntots", "4,8,12")
    years = parse(Int, arg("years", "1"))
    steps = parse(Int, arg("steps", "365"))
    repetitions = parse(Int, arg("repetitions", "3"))
    output = arg("output", "snowpack_columns.csv")

    backend == :gpu && !Chion.cuda_available() && error("CUDA is not available")
    run!(simulation(2, first(ntots), 1, 2, backend); io=devnull)

    for ntot in ntots, ncol in columns, repetition in 1:repetitions
        sim = simulation(ncol, ntot, years, steps, backend)
        GC.gc()
        result = nothing
        seconds = @elapsed begin
            result = run!(sim; io=devnull)
            backend == :gpu && Chion.CUDA.synchronize()
        end
        open(output, "a") do io
            println(io, join((backend, Threads.nthreads(), ncol, ntot, repetition, seconds, result.status), ","))
        end
        println("$backend threads=$(Threads.nthreads()) columns=$ncol ntot=$ntot: $(round(seconds; digits=3)) s")
    end
end

function main()
    output = abspath(arg("output", "snowpack_columns.csv"))
    mkpath(dirname(output))
    write(output, "backend,threads,columns,ntot,repetition,seconds,status\n")

    common = [
        "--worker",
        "--columns=$(arg("columns", "128,1024,16384,65536,524288,1048576"))",
        "--ntots=$(arg("ntots", "2,4,8,12,20"))",
        "--years=$(arg("years", "100"))",
        "--steps=$(arg("steps", "365"))",
        "--repetitions=$(arg("repetitions", "3"))",
        "--output=$output",
    ]
    for backend in split(arg("backends", "threads"), ",")
        thread_counts = backend == "threads" ? list("threads", "32") : [1]
        for threads in thread_counts
            run(`$(Base.julia_cmd()) --threads=$threads --project=$ROOT $SCRIPT --backend=$backend $(common)`)
        end
    end
    println("Results written to $output")
end

"--worker" in ARGS ? worker() : main()
