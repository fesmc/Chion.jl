#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

include("gris_equilibrium_backend.jl")

function print_compare_help()
    println("Usage:")
    println("  julia --project=. examples/scripts/compare_gris_equilibrium_api.jl [options]")
    println()
    println("This compares the existing examples/scripts/run_gris_equilibrium.jl runner")
    println("against the new API-backed examples/scripts/run_gris_equilibrium_api.jl runner")
    println("in an old-compatible configuration.")
    println()
    println("Options:")
    println("  --nc=PATH                    MAR NetCDF/HDF5 file")
    println("  --backend=threads|cpu|gpu    Backend to compare (default: threads)")
    println("  --max-cycles=N               Max cycles for both runners (default: 1)")
    println("  --mask-threshold=VALUE       Mask threshold passed to both runners (default: 50)")
    println("  --write-output               Enable summary/plot/CSV outputs in both runners")
    println("  --write-nc                   Enable NetCDF output in both runners")
    println("  --flip-turbulent-fluxes      Multiply SHF and LHF by -1 in both runners")
    println("  --help                       Show this message")
end

function capture_command(cmd::Cmd)
    out = PipeBuffer()
    err = PipeBuffer()
    proc = run(pipeline(ignorestatus(cmd), stdout=out, stderr=err))
    return (
        exitcode = proc.exitcode,
        stdout = String(take!(out)),
        stderr = String(take!(err)),
        command = sprint(show, cmd),
    )
end

function extract_metric(name::AbstractString, output::AbstractString)
    match_obj = match(Regex("^" * name * "\\s*:\\s*([^\\n]+)\$", "m"), output)
    return isnothing(match_obj) ? nothing : strip(match_obj.captures[1])
end

function extract_float_metric(name::AbstractString, output::AbstractString)
    value = extract_metric(name, output)
    return isnothing(value) ? nothing : parse(Float64, replace(value, " s" => ""))
end

function extract_cycle_lines(output::AbstractString)
    return [line for line in split(output, '\n') if startswith(line, "cycle=")]
end

function first_difference(old_lines::AbstractVector{<:AbstractString}, new_lines::AbstractVector{<:AbstractString})
    n = min(length(old_lines), length(new_lines))
    for idx in 1:n
        old_lines[idx] == new_lines[idx] || return (idx=idx, old=old_lines[idx], new=new_lines[idx])
    end
    if length(old_lines) != length(new_lines)
        idx = n + 1
        old_line = idx <= length(old_lines) ? old_lines[idx] : "<missing>"
        new_line = idx <= length(new_lines) ? new_lines[idx] : "<missing>"
        return (idx=idx, old=old_line, new=new_line)
    end
    return nothing
end

function build_runner_args(args::Vector{String})
    backend_name = lowercase(strip(arg_value(args, "backend", "threads")))
    backend = backend_name == "cpu" ? :threads : Symbol(backend_name)
    backend in (:threads, :gpu) || error("Unsupported backend '$backend_name'. Use `threads`, `cpu`, or `gpu`.")
    write_output = has_flag(args, "write-output")
    write_netcdf = has_flag(args, "write-nc")
    !write_output && write_netcdf && error("The old runner does not support NetCDF-only mode. Use --write-output together with --write-nc, or omit --write-nc.")

    nc_path = arg_value(args, "nc", DEFAULT_GRIS_API_NC_PATH)
    isempty(nc_path) && error("Pass --nc=PATH or place the MAR file at $(DEFAULT_GRIS_API_NC_PATH).")

    old_args = String[
        "--nc=$(nc_path)",
        "--backend=$(backend == :threads ? "threads" : "gpu")",
        "--max-cycles=$(arg_value(args, "max-cycles", "1"))",
        "--mask-threshold=$(arg_value(args, "mask-threshold", "50.0"))",
    ]
    new_args = copy(old_args)

    if !write_output
        push!(old_args, "--no-output")
        push!(new_args, "--no-output")
    end
    if !write_netcdf
        push!(new_args, "--no-nc")
        if write_output
            push!(old_args, "--no-nc")
        end
    end
    if has_flag(args, "flip-turbulent-fluxes")
        push!(old_args, "--flip-turbulent-fluxes")
        push!(new_args, "--flip-turbulent-fluxes")
    end
    return old_args, new_args
end

function summarize_run(label::AbstractString, output::AbstractString)
    return (
        label = label,
        status = extract_metric("Status", output),
        cycles = extract_metric("Cycles", output),
        simulation_wall = extract_float_metric("Simulation wall", output),
        run_wall_total = extract_float_metric("Run wall total", output),
        model_step_wall = extract_float_metric("Model step wall", output),
        cycle_lines = extract_cycle_lines(output),
    )
end

function main(args::Vector{String})
    if has_flag(args, "help")
        print_compare_help()
        return
    end

    ensure_tools!()

    old_args, new_args = build_runner_args(args)
    project_dir = joinpath(@__DIR__, "..", "..")
    old_script = joinpath(@__DIR__, "run_gris_equilibrium.jl")
    new_script = joinpath(@__DIR__, "run_gris_equilibrium_api.jl")

    old_run = capture_command(`julia --project=$(project_dir) $(old_script) $(old_args)`)
    old_run.exitcode == 0 || error("Old runner failed.\nCommand: $(old_run.command)\nSTDERR:\n$(old_run.stderr)\nSTDOUT:\n$(old_run.stdout)")

    new_run = capture_command(`julia --project=$(project_dir) $(new_script) $(new_args)`)
    new_run.exitcode == 0 || error("API runner failed.\nCommand: $(new_run.command)\nSTDERR:\n$(new_run.stderr)\nSTDOUT:\n$(new_run.stdout)")

    old_summary = summarize_run("old", old_run.stdout)
    new_summary = summarize_run("api", new_run.stdout)
    diff = first_difference(old_summary.cycle_lines, new_summary.cycle_lines)
    same_status = old_summary.status == new_summary.status
    same_cycles = old_summary.cycles == new_summary.cycles
    same_logs = isnothing(diff)

    println("Comparison complete.")
    println("Old runner : run_wall_total=$(old_summary.run_wall_total) s simulation_wall=$(old_summary.simulation_wall) s status=$(old_summary.status) cycles=$(old_summary.cycles)")
    println("API runner : run_wall_total=$(new_summary.run_wall_total) s simulation_wall=$(new_summary.simulation_wall) s status=$(new_summary.status) cycles=$(new_summary.cycles)")
    if !isnothing(old_summary.run_wall_total) && !isnothing(new_summary.run_wall_total) && new_summary.run_wall_total > 0
        println(@sprintf("Wall ratio  : %.3fx (old/api)", old_summary.run_wall_total / new_summary.run_wall_total))
        println(@sprintf("Wall delta  : %.3f s (api - old)", new_summary.run_wall_total - old_summary.run_wall_total))
    end

    if same_status && same_cycles && same_logs
        println("Results     : cycle logs, status, and completed cycle count match exactly.")
        return
    end

    println("Results     : mismatch detected.")
    println("Status match: $(same_status)")
    println("Cycles match: $(same_cycles)")
    println("Logs match  : $(same_logs)")
    if !isnothing(diff)
        println("First diff  : cycle line $(diff.idx)")
        println("Old         : $(diff.old)")
        println("API         : $(diff.new)")
    end
    exit(1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
