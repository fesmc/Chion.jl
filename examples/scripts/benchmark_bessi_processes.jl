#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using BenchmarkTools
using Chion

function print_help()
    println("Usage:")
    println("  julia --project=. examples/scripts/benchmark_bessi_processes.jl [options]")
    println()
    println("Options:")
    println("  --seconds=N      Approximate seconds per benchmark (default: 3)")
    println("  --samples=N      Maximum samples per benchmark (default: 10_000)")
    println("  --evals=N        Evaluations per sample (default: 1)")
end

arg_value(args, name, default) = begin
    prefix = "--$(name)="
    for arg in args
        startswith(arg, prefix) && return arg[length(prefix)+1:end]
    end
    return default
end

has_flag(args, name) = any(==("--$(name)"), args)

function seed_domain!(domain)
    c = domain.c
    fill!(domain.N, 4)
    fill!(domain.mass, 0.0)
    fill!(domain.mass_w, 0.0)
    fill!(domain.density, 300.0)
    fill!(domain.temperature, c.T0 - 12.0)
    fill!(domain.mass_base, 0.0)
    fill!(domain.smb_ice, 0.0)
    fill!(domain.runoff, 0.0)
    fill!(domain.Tsrf, c.T0 - 12.0)
    fill!(domain.snow_cover, 1.0)
    fill!(domain.albedo_dynamic, c.alpha_dry)

    domain.mass[1, 1] = 85.0
    domain.mass[2, 1] = 120.0
    domain.mass[3, 1] = 160.0
    domain.mass[4, 1] = 210.0
    domain.mass_w[1, 1] = 1.5
    domain.mass_w[2, 1] = 0.5
    domain.density[1, 1] = 290.0
    domain.density[2, 1] = 330.0
    domain.density[3, 1] = 380.0
    domain.density[4, 1] = 450.0
    domain.temperature[1, 1] = c.T0 - 6.0
    domain.temperature[2, 1] = c.T0 - 10.0
    domain.temperature[3, 1] = c.T0 - 14.0
    domain.temperature[4, 1] = c.T0 - 18.0
    domain.Tsrf[1] = domain.temperature[1, 1]
    return domain
end

function fresh_domain(template)
    domain = Chion.SnowpackDomain(;
        c=template.c,
        Ntot=template.Ntot,
        ncol=template.ncol,
        mass_max=template.mass_max,
        mass_split=template.mass_split,
        mass_min=template.mass_min,
        rho_max=template.rho_max,
    )
    return Chion._copy_domain_state!(domain, template)
end

function benchmark_suite(; seconds=3.0, samples=10_000, evals=1)
    grid = SnowpackGrid(1)
    model = BESSIModel(grid; Ntot=8)
    template = seed_domain!(initial_state(model).domain)
    c = template.c

    dt_days = 1.0
    dt_seconds = dt_days * c.seconds_per_day
    air_temperature = c.T0 - 8.0
    snowfall_rate = 2.0 / c.seconds_per_day
    rainfall_rate = 0.0
    shortwave_down = 180.0
    wind_speed = 5.0
    forcing = Chion.SnowpackStepForcing(
        air_temperature,
        snowfall_rate + rainfall_rate,
        dt_days,
        snowfall_rate,
        rainfall_rate,
        shortwave_down,
        wind_speed,
        0.0,
        0.0,
        0.0,
        0.0,
        false,
        false,
        false,
        false,
        false,
        0.0,
        0.0,
    )
    latent_linear, latent_constant = Chion._diagnose_latent_heat_flux_coefficients(
        true,
        c,
        air_temperature,
        snowfall_rate,
        rainfall_rate,
    )

    benchmarks = Pair{String,Any}[]

    push!(benchmarks, "accumulation" => (@benchmarkable Chion._apply_accumulation!(
            d.N, d.mass, d.mass_w, d.density, d.temperature, d.mass_base, d.smb_ice,
            d.runoff, d.Tsrf, d.snow_cover, d.albedo_dynamic, 1, d.c, d.Ntot,
            d.mass_max, d.mass_split, d.mass_min, $snowfall_rate, $rainfall_rate,
            $dt_seconds; air_temperature=$air_temperature, wind_speed=$wind_speed,
        ) setup=(d = fresh_domain($template)) seconds=seconds samples=samples evals=evals))

    push!(benchmarks, "surface_albedo" => (@benchmarkable Chion._update_surface_albedo_arrays!(
            d.N, d.mass, d.mass_w, d.density, d.temperature, d.albedo_dynamic, 1, d.c,
        ) setup=(d = fresh_domain($template)) seconds=seconds samples=samples evals=evals))

    push!(benchmarks, "densification" => (@benchmarkable Chion._go_densification!(
            d.N, d.mass, d.density, d.temperature, 1, d.c, $snowfall_rate, $dt_seconds,
        ) setup=(d = fresh_domain($template)) seconds=seconds samples=samples evals=evals))

    push!(benchmarks, "energy_flux" => (@benchmarkable Chion._go_energy_flux_resolved!(
            d.N, d.mass, d.mass_w, d.density, d.temperature, d.Tsrf, d.albedo_dynamic,
            1, d.c, w.energy, $air_temperature, $shortwave_down, $latent_linear,
            $latent_constant, $dt_seconds, false, 0.0, false, 0.0, false, 0.0,
            false, 0.0,
        ) setup=(d = fresh_domain($template); w = Chion.ColumnarStepWorkspace(d)) seconds=seconds samples=samples evals=evals))

    push!(benchmarks, "melt" => (@benchmarkable Chion._apply_melt!(
            d.N, d.mass, d.mass_w, d.density, d.temperature, d.runoff, d.Tsrf,
            d.albedo_dynamic, 1, d.mass_split, d.mass_min, 0.25, d.c,
        ) setup=(d = fresh_domain($template)) seconds=seconds samples=samples evals=evals))

    push!(benchmarks, "percolation" => (@benchmarkable Chion._go_percolation!(
            d.N, d.mass, d.mass_w, d.density, 1, d.c.rho_i, d.c.rho_w,
        ) setup=(d = fresh_domain($template)) seconds=seconds samples=samples evals=evals))

    push!(benchmarks, "refreezing" => (@benchmarkable Chion._go_refreezing!(
            d.N, d.mass_w, d.mass, d.density, d.temperature, 1, d.c.T0, d.c.ci,
            d.c.Lm, d.c.rho_i,
        ) setup=(d = fresh_domain($template)) seconds=seconds samples=samples evals=evals))

    push!(benchmarks, "full_column_step" => (@benchmarkable Chion._step_state_resolved!(
            d.N, d.mass, d.mass_w, d.density, d.temperature, d.mass_base, d.smb_ice,
            d.runoff, d.melt, d.refreezing, d.Tsrf, d.snow_cover, d.albedo_dynamic, 1, d.c, d.Ntot,
            d.mass_max, d.mass_split, d.mass_min, $forcing, w, true,
        ) setup=(d = fresh_domain($template); w = Chion.ColumnarStepWorkspace(d)) seconds=seconds samples=samples evals=evals))

    return benchmarks
end

function main(args::Vector{String})
    if has_flag(args, "help")
        print_help()
        return
    end

    seconds = parse(Float64, arg_value(args, "seconds", "3"))
    samples = parse(Int, replace(arg_value(args, "samples", "10000"), "_" => ""))
    evals = parse(Int, arg_value(args, "evals", "1"))

    println("Benchmarking BESSI CPU process kernels")
    println("seconds=$(seconds), samples=$(samples), evals=$(evals)")
    println()

    for (name, benchmark) in benchmark_suite(; seconds=seconds, samples=samples, evals=evals)
        println("== ", name, " ==")
        trial = run(benchmark)
        show(stdout, MIME("text/plain"), trial)
        println()
        println()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
