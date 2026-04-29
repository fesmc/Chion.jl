"""Initialized output runtime setup, step writes, and final output writes."""

function _generic_final_output_vector(model::AbstractSnowModel, state::AbstractSnowModelState, runtime, spec, final_state, deltas)
    source = spec.source
    values = source === :final_state ? getfield(final_state, spec.field) :
        source === :deltas ? getfield(deltas, spec.field) :
        source === :domain && spec.field === :smb_ice ? _model_smb_ice_vector(model, state, runtime) :
        source === :domain && spec.field === :runoff ? _model_runoff_vector(model, state, runtime) :
        error("Unsupported final output source `$(source)`.")
    return _host_vector(values; copy_array=true)
end

model_final_output_vector(model::AbstractSnowModel, state::AbstractSnowModelState, runtime, spec, final_state, deltas) =
    _generic_final_output_vector(model, state, runtime, spec, final_state, deltas)

function model_final_output_vector(model::PDDModel, state::PDDState, runtime, spec, final_state, deltas)
    spec.key === :final_base_mass && return zeros(Float64, ncols(model.grid))
    return _generic_final_output_vector(model, state, runtime, spec, final_state, deltas)
end

function _scatter_model_final_grids(model::AbstractSnowModel, state::AbstractSnowModelState, runtime, final_state, deltas, layout)
    return NamedTuple{FINAL_GRID_KEYS}(ntuple(i -> begin
        values = model_final_output_vector(model, state, runtime, FINAL_OUTPUT_SPECS[i], final_state, deltas)
        scatter_to_grid(values, layout.js, layout.is, _grid_shape(layout))
    end, length(FINAL_OUTPUT_SPECS)))
end

model_layer_grids(::AbstractSnowModel, ::AbstractSnowModelState, runtime, layout, need_layer_outputs::Bool, timings::StepTimingStats) =
    _empty_layer_grids()

function model_layer_grids(::BESSIModel, ::BESSIState, runtime, layout, need_layer_outputs::Bool, timings::StepTimingStats)
    need_layer_outputs || return _empty_layer_grids()
    final_domain = runtime.is_gpu ? time_block!(timings, :gpu_transfer) do
        cpu_domain(runtime.domain)
    end : runtime.domain
    return time_block!(timings, :collect_final_layer_grids) do
        collect_final_layer_grids(final_domain, layout.js, layout.is, _grid_shape(layout), final_domain.Ntot)
    end
end

function init_io!(context, model_runtime::IntegratorModelRuntime, diagnostics::DiagnosticsRuntime, options::RunOptions, timings::StepTimingStats)
    schedule = nothing
    writer = nothing
    nc_path = ""
    if options.write_netcdf
        schedule = time_block!(timings, :prepare_output_schedule) do
            _prepare_output_schedule(context.forcing.time_values, options.years)
        end
        initial_thickness_vec = copy(diagnostics.prev.thickness)
        initial_thickness = time_block!(timings, :prepare_initial_output_fields_grid) do
            scatter_to_grid(initial_thickness_vec, model_runtime.grid.js, model_runtime.grid.is, _grid_shape(model_runtime.grid))
        end
        nc_path = resolve_netcdf_path(options)
        writer = time_block!(timings, :init_netcdf) do
            init_netcdf(
                nc_path,
                options,
                context.forcing.time_values,
                model_layer_count(context.model, context.state, model_runtime.runtime),
                model_runtime.grid,
                initial_thickness,
                schedule.month_of_year,
                schedule.source_month_code,
                schedule.annual_output.source_indices,
                schedule.annual_output.source_codes,
            )
        end
    end

    monthly_total = diagnostics.need_monthly_outputs && schedule !== nothing ? schedule.nmonth_total : 0
    monthly_sums = _allocate_monthly_sums(diagnostics.need_monthly_outputs, monthly_total, model_runtime.ncol)
    monthly_count = diagnostics.need_monthly_outputs ? zeros(Int32, monthly_total) : Int32[]
    step_vectors = _allocate_step_vectors(diagnostics.need_step_outputs, model_runtime.ncol)
    return OutputRuntime(schedule, writer, nc_path, monthly_sums, monthly_count, step_vectors, 0)
end

function maybe_write_step_outputs!(integrator::SimulationIntegrator)
    diagnostics = integrator.diagnostics
    output = integrator.output
    stepper = integrator.stepper
    diagnostics.need_step_outputs || return nothing
    output.schedule.annual_output.write_output[stepper.time_index] || return nothing

    output.steps_written += 1
    if output.writer !== nothing
        step_grids = time_block!(integrator.timings, :step_output_prepare) do
            _step_output_grids(output.step_vectors, integrator.model_runtime.grid)
        end
        time_block!(integrator.timings, :step_output_write) do
            for key in OUTPUT_GROUPS.step
                maybe_write_step_output!(output.writer, output.steps_written, key, getfield(step_grids, key))
            end
        end
        _reset_step_vectors!(output.step_vectors)
    end
    return nothing
end

function finalize_output_runtime!(integrator::SimulationIntegrator, status::Symbol, years_completed::Int)
    output = integrator.output
    output.writer === nothing && return nothing

    model_runtime = integrator.model_runtime
    diagnostics = integrator.diagnostics
    history = diagnostics.history
    final_state = _current_year_summary!(integrator)
    final_grids = time_block!(integrator.timings, :scatter_final_outputs) do
        _scatter_model_final_grids(
            integrator.sim.model,
            integrator.sim.now,
            model_runtime.runtime,
            final_state,
            diagnostics.deltas,
            model_runtime.grid,
        )
    end
    layer_grids = model_layer_grids(
        integrator.sim.model,
        integrator.sim.now,
        model_runtime.runtime,
        model_runtime.grid,
        diagnostics.need_layer_outputs,
        integrator.timings,
    )
    monthly_grids = diagnostics.need_monthly_outputs ? time_block!(integrator.timings, :aggregate_monthly_outputs) do
        _finalize_monthly_grids(output.monthly_sums, output.monthly_count, model_runtime.grid)
    end : empty_monthly_grids()
    time_block!(integrator.timings, :write_netcdf) do
        finalize_netcdf!(output.writer, final_grids, layer_grids, history, monthly_grids, status, years_completed, output.steps_written)
    end
    return nothing
end
