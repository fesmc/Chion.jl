"""Checkpoint and restart helpers for initialized `Simulation` runs."""

const CHECKPOINT_FORMAT_VERSION = 2

function _checkpoint_run_options(options::RunOptions)
    return (
        name=options.name,
        input_label=options.input_label,
        output_dir=options.output_dir,
        netcdf_path=options.netcdf_path,
        write_outputs=options.write_outputs,
        write_netcdf=options.write_netcdf,
        netcdf_variables=copy(options.netcdf_variables),
        years=options.years,
        backend=options.backend,
        history_year_stride=options.history_year_stride,
    )
end

function _run_options_from_checkpoint(options)
    return RunOptions(
        name=options.name,
        input_label=options.input_label,
        output_dir=options.output_dir,
        netcdf_path=options.netcdf_path,
        write_outputs=options.write_outputs,
        write_netcdf=options.write_netcdf,
        netcdf_variables=options.netcdf_variables,
        years=options.years,
        backend=options.backend,
        history_year_stride=options.history_year_stride,
    )
end

function _checkpoint_diagnostics(diagnostics::DiagnosticsRuntime)
    return (
        prev=deepcopy(diagnostics.prev),
        final=deepcopy(diagnostics.final),
        step_summary=deepcopy(diagnostics.step_summary),
        deltas=deepcopy(diagnostics.deltas),
        previous=deepcopy(diagnostics.previous),
        previous_year_smb_ice=copy(diagnostics.previous_year_smb_ice),
        history=deepcopy(diagnostics.history),
    )
end

function _restore_diagnostics!(diagnostics::DiagnosticsRuntime, checkpoint)
    diagnostics.prev = deepcopy(checkpoint.prev)
    diagnostics.final = deepcopy(checkpoint.final)
    diagnostics.step_summary = deepcopy(checkpoint.step_summary)
    diagnostics.deltas = deepcopy(checkpoint.deltas)
    diagnostics.previous = deepcopy(checkpoint.previous)
    diagnostics.previous_year_smb_ice = copy(checkpoint.previous_year_smb_ice)
    diagnostics.history = deepcopy(checkpoint.history)
    return diagnostics
end

function _checkpoint_output(output::OutputRuntime)
    return (
        nc_path=output.nc_path,
        monthly_sums=deepcopy(output.monthly_sums),
        monthly_count=copy(output.monthly_count),
        step_vectors=deepcopy(output.step_vectors),
        steps_written=output.steps_written,
    )
end

function _netcdf_var_name(key::Symbol)
    for spec in NC_SPECS
        spec.key == key && return spec.name
    end
    return ""
end

function _netcdf_dim_capacity(ds, name::AbstractString)
    haskey(ds, name) && return length(ds[name])
    return typemax(Int)
end

function _netcdf_can_store_restart(ds, options::RunOptions, schedule)
    _netcdf_dim_capacity(ds, "year") >= options.years || return false
    schedule === nothing && return true
    _netcdf_dim_capacity(ds, "month") >= schedule.nmonth_total || return false
    expected_steps = options.years * length(schedule.annual_output.source_indices)
    _netcdf_dim_capacity(ds, "step") >= expected_steps || return false
    return true
end

function _continued_netcdf_path(nc_path::AbstractString, years::Int)
    root, ext = splitext(nc_path)
    isempty(ext) && (ext = ".nc")
    for suffix in Iterators.flatten((("",), ("_$(i)" for i in 2:1000)))
        candidate = "$(root)_to$(years)y$(suffix)$(ext)"
        isfile(candidate) || return candidate
    end
    error("Could not find an unused continued NetCDF path for '$(nc_path)'.")
end

function _reopen_netcdf_writer(nc_path::AbstractString, options::RunOptions, schedule)
    options.write_netcdf || return nothing
    isempty(nc_path) && error("Checkpoint requested NetCDF restart, but no NetCDF path was stored.")
    isfile(nc_path) || error("Cannot restart NetCDF output because '$(abspath(nc_path))' does not exist.")

    ds = NCDataset(nc_path, "a")
    vars = Dict{Symbol, Any}()
    try
        if !_netcdf_can_store_restart(ds, options, schedule)
            error(
                "Checkpoint NetCDF '$(abspath(nc_path))' is too small for restart target years=$(options.years). " *
                "Pass --netcdf-path=PATH for a new extended output file, or disable NetCDF with --no-nc.",
            )
        end
        for key in options.netcdf_variables
            name = _netcdf_var_name(key)
            !isempty(name) && haskey(ds, name) && (vars[key] = ds[name])
        end
        haskey(ds, "step_valid") && (vars[:step_valid] = ds["step_valid"])
        max_steps = haskey(ds, "step") ? length(ds["step"]) : options.years
        return NetCDFWriter(ds, vars, max_steps, options.years)
    catch
        close(ds)
        rethrow()
    end
end

function _extend_monthly_sums(monthly_sums, active::Bool, nmonth_total::Int, ncol::Int)
    active || return _allocate_monthly_sums(false, 0, ncol)
    return NamedTuple{MONTHLY_GRID_KEYS}(ntuple(i -> begin
        key = MONTHLY_GRID_KEYS[i]
        previous = getfield(monthly_sums, key)
        if isempty(previous)
            error("Restart requested monthly NetCDF outputs, but the checkpoint has no saved monthly accumulator for $(key).")
        end
        size(previous, 2) == ncol || error("Checkpoint monthly accumulator $(key) has $(size(previous, 2)) columns; expected $(ncol).")
        size(previous, 1) <= nmonth_total || error("Checkpoint monthly accumulator $(key) has more months than restart target years allow.")
        extended = zeros(Float64, nmonth_total, ncol)
        extended[1:size(previous, 1), :] .= previous
        extended
    end, length(MONTHLY_GRID_KEYS)))
end

function _extend_monthly_count(monthly_count::Vector{Int32}, active::Bool, nmonth_total::Int)
    active || return Int32[]
    isempty(monthly_count) && error("Restart requested monthly NetCDF outputs, but the checkpoint has no saved monthly counts.")
    length(monthly_count) <= nmonth_total || error("Checkpoint monthly counts have more months than restart target years allow.")
    extended = zeros(Int32, nmonth_total)
    extended[1:length(monthly_count)] .= monthly_count
    return extended
end

function _new_restart_netcdf_writer(
    nc_path::AbstractString,
    sim,
    model_runtime::ModelRuntime,
    diagnostics::DiagnosticsRuntime,
    options::RunOptions,
    schedule,
)
    initial_thickness = scatter_to_grid(diagnostics.prev.thickness, model_runtime.grid.js, model_runtime.grid.is, _grid_shape(model_runtime.grid))
    return init_netcdf(
        nc_path,
        options,
        sim.forcing.time_values,
        model_layer_count(sim.model, sim.now, model_runtime.backend),
        model_runtime.grid,
        initial_thickness,
        schedule.month_of_year,
        schedule.source_month_code,
        schedule.annual_output.source_indices,
        schedule.annual_output.source_codes,
    )
end

function _restart_output_runtime(sim, model_runtime::ModelRuntime, diagnostics::DiagnosticsRuntime, options::RunOptions, checkpoint)
    schedule = options.write_netcdf ? _prepare_output_schedule(sim.forcing.time_values, options.years) : nothing
    writer = nothing
    nc_path = checkpoint.nc_path
    if options.write_netcdf
        requested_path = resolve_netcdf_path(options)
        explicit_new_path = !isempty(options.netcdf_path) && abspath(requested_path) != abspath(checkpoint.nc_path)
        if explicit_new_path
            nc_path = requested_path
            if diagnostics.need_step_outputs && checkpoint.steps_written > 0
                error("Cannot move checkpointed step NetCDF output to a new file because earlier step grids are not stored in the checkpoint.")
            end
            writer = _new_restart_netcdf_writer(nc_path, sim, model_runtime, diagnostics, options, schedule)
        else
            try
                writer = _reopen_netcdf_writer(checkpoint.nc_path, options, schedule)
            catch err
                if !isa(err, ErrorException) || occursin("--netcdf-path", sprint(showerror, err)) == false
                    rethrow()
                end
                nc_path = _continued_netcdf_path(checkpoint.nc_path, options.years)
                if diagnostics.need_step_outputs && checkpoint.steps_written > 0
                    error("Cannot auto-create continued NetCDF with step outputs because earlier step grids are not stored in the checkpoint. Pass --no-nc or omit step outputs.")
                end
                println("Restart NetCDF extension: writing continued output to $(nc_path)")
                writer = _new_restart_netcdf_writer(nc_path, sim, model_runtime, diagnostics, options, schedule)
            end
        end
    end
    nmonth_total = schedule === nothing ? 0 : schedule.nmonth_total
    return OutputRuntime(
        schedule,
        writer,
        nc_path,
        _extend_monthly_sums(checkpoint.monthly_sums, diagnostics.need_monthly_outputs, nmonth_total, model_runtime.ncol),
        _extend_monthly_count(checkpoint.monthly_count, diagnostics.need_monthly_outputs, nmonth_total),
        deepcopy(checkpoint.step_vectors),
        checkpoint.steps_written,
    )
end

function _sync_checkpoint_state!(integrator::SimulationIntegrator)
    finalize_state!(
        integrator.sim.model,
        integrator.sim.now,
        integrator.model_runtime.backend,
        integrator.options,
        integrator.timings,
    )
    return nothing
end

function _sync_checkpoint_output!(integrator::SimulationIntegrator)
    integrator.output.writer === nothing && return nothing
    try
        sync(integrator.output.writer.dataset)
    catch
        nothing
    end
    return nothing
end

function _checkpoint_payload(integrator::SimulationIntegrator)
    return (
        format=:ChionSimulationCheckpoint,
        version=CHECKPOINT_FORMAT_VERSION,
        sim=deepcopy(integrator.sim),
        options=_checkpoint_run_options(integrator.options),
        time_index=integrator.time_index,
        completed_years=integrator.completed_years,
        current_forcing=deepcopy(integrator.current_forcing),
        diagnostics=_checkpoint_diagnostics(integrator.diagnostics),
        output=_checkpoint_output(integrator.output),
    )
end

function _checkpoint_parent(path::AbstractString)
    parent = dirname(path)
    return isempty(parent) ? "." : parent
end

"""
    checkpoint!(integrator, path)

Write a restart checkpoint for an initialized `SimulationIntegrator`.

The checkpoint stores serializable model state, run counters, diagnostics, and
output aggregation buffers. Live runtime resources such as GPU arrays, progress
meters, and NetCDF file handles are rebuilt by `restart_integrator`.
"""
function checkpoint!(integrator::SimulationIntegrator, path::AbstractString)
    integrator.finalized && error("Cannot checkpoint a finalized integrator.")
    _sync_checkpoint_state!(integrator)
    _sync_checkpoint_output!(integrator)
    mkpath(_checkpoint_parent(path))
    tmp_path = string(path, ".tmp")
    open(tmp_path, "w") do io
        serialize(io, _checkpoint_payload(integrator))
    end
    mv(tmp_path, path; force=true)
    return String(path)
end

function _read_checkpoint(path::AbstractString)
    checkpoint = open(deserialize, path)
    checkpoint.format == :ChionSimulationCheckpoint || error("Unsupported checkpoint format in '$(abspath(path))'.")
    checkpoint.version == CHECKPOINT_FORMAT_VERSION || error("Unsupported checkpoint version $(checkpoint.version); expected $(CHECKPOINT_FORMAT_VERSION).")
    return checkpoint
end

"""
    restart_integrator(path; io=stdout)

Load a checkpoint written by `checkpoint!` and return an initialized
`SimulationIntegrator` that can continue with `run!`, `step!`, or `finalize!`.
"""
function _restart_integrator_from_checkpoint(checkpoint; io::IO=stdout, sim_transform=identity, options_transform=(options, checkpoint) -> options)
    sim = sim_transform(deepcopy(checkpoint.sim))
    options = options_transform(_run_options_from_checkpoint(checkpoint.options), checkpoint)
    timings = StepTimingStats()
    init_problem!(sim, options)
    model_runtime = init_model_runtime!(sim, options, timings)
    diagnostics = init_diagnostics!(sim, model_runtime, options, timings)
    _restore_diagnostics!(diagnostics, checkpoint.diagnostics)
    output = _restart_output_runtime(sim, model_runtime, diagnostics, options, checkpoint.output)
    return _new_integrator(
        sim,
        options,
        io,
        timings,
        model_runtime,
        diagnostics,
        output;
        time_index=checkpoint.time_index,
        completed_years=checkpoint.completed_years,
        current_forcing=deepcopy(checkpoint.current_forcing),
    )
end

function restart_integrator(path::AbstractString; io::IO=stdout, sim_transform=identity, options_transform=(options, checkpoint) -> options)
    checkpoint = _read_checkpoint(path)
    return _restart_integrator_from_checkpoint(checkpoint; io=io, sim_transform=sim_transform, options_transform=options_transform)
end

load_checkpoint(path::AbstractString; io::IO=stdout) = restart_integrator(path; io=io)
