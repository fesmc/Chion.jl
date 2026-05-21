"""Runtime summaries, step diagnostics, history, and yearly deltas."""

const SUMMARY_BUFFER_NAMES = (:thickness, :wet_mass, :bulk_density, :base_mass, :smb_ice, :liquid_water, :runoff, :pdd, :melt, :refreezing, :albedo)
const YEAR_BUFFER_NAMES = (:thickness, :wet_mass, :bulk_density, :base_mass)

_named_buffers(names::NTuple{N, Symbol}, build::F) where {N, F <: Function} = NamedTuple{names}(ntuple(_ -> build(), N))

allocate_summary_buffers(n::Int) = _named_buffers(SUMMARY_BUFFER_NAMES, () -> Vector{Float64}(undef, n))
allocate_summary_buffers(domain::AbstractSnowpackDomain, n::Int) = _named_buffers(SUMMARY_BUFFER_NAMES, () -> similar(domain.mass, Float64, n))
allocate_year_summary_buffers(n::Int) = _named_buffers(YEAR_BUFFER_NAMES, () -> Vector{Float64}(undef, n))
allocate_year_summary_buffers(domain::AbstractSnowpackDomain, n::Int) = _named_buffers(YEAR_BUFFER_NAMES, () -> similar(domain.mass, Float64, n))

_allocate_year_backend_buffers(::AbstractSnowModel, ::AbstractSnowModelState, runtime, ncol::Int) = allocate_year_summary_buffers(ncol)
_allocate_year_backend_buffers(::BESSIModel, ::BESSIState, runtime, ncol::Int) = allocate_year_summary_buffers(runtime.domain, ncol)
_allocate_year_backend_buffers(::PDDModel, ::PDDState, runtime, ncol::Int) =
    _named_buffers(YEAR_BUFFER_NAMES, () -> similar(runtime.snowpack_swe, Float64, ncol))

_allocate_step_backend_buffers(::AbstractSnowModel, ::AbstractSnowModelState, runtime, ncol::Int) = allocate_summary_buffers(ncol)
_allocate_step_backend_buffers(::BESSIModel, ::BESSIState, runtime, ncol::Int) = allocate_summary_buffers(runtime.domain, ncol)
_allocate_step_backend_buffers(::PDDModel, ::PDDState, runtime, ncol::Int) =
    _named_buffers(SUMMARY_BUFFER_NAMES, () -> similar(runtime.snowpack_swe, Float64, ncol))

function summarize_year_state!(summary, ::BESSIModel, ::BESSIState, runtime)
    summarize_year_state!(
        summary.thickness,
        summary.wet_mass,
        summary.bulk_density,
        summary.base_mass,
        runtime.domain,
    )
    return summary
end

function summarize_year_state!(summary, ::PDDModel, ::PDDState, runtime)
    summary.thickness .= runtime.snowpack_swe ./ 1000.0
    summary.wet_mass .= runtime.snowpack_swe
    fill!(summary.bulk_density, NaN)
    summary.base_mass .= runtime.smb_ice
    return summary
end

function summarize_step_state!(summary, ::BESSIModel, ::BESSIState, runtime)
    summarize_domain_state!(
        summary.thickness,
        summary.wet_mass,
        summary.bulk_density,
        summary.base_mass,
        summary.smb_ice,
        summary.liquid_water,
        summary.runoff,
        summary.melt,
        summary.refreezing,
        summary.albedo,
        runtime.domain,
    )
    fill!(summary.pdd, 0.0)
    return summary
end

function summarize_step_state!(summary, ::PDDModel, ::PDDState, runtime)
    summary.thickness .= runtime.snowpack_swe ./ 1000.0
    summary.wet_mass .= runtime.snowpack_swe
    fill!(summary.bulk_density, NaN)
    summary.base_mass .= runtime.smb_ice
    summary.smb_ice .= runtime.smb_ice
    fill!(summary.liquid_water, 0.0)
    summary.runoff .= runtime.runoff
    summary.pdd .= runtime.pdd_sum
    fill!(summary.melt, NaN)
    fill!(summary.refreezing, NaN)
    fill!(summary.albedo, NaN)
    return summary
end

_model_smb_ice_vector(::BESSIModel, ::BESSIState, runtime) = _host_vector(runtime.domain.smb_ice; copy_array=true)
_model_smb_ice_vector(::PDDModel, ::PDDState, runtime) = _host_vector(runtime.smb_ice; copy_array=true)
_model_runoff_vector(::BESSIModel, ::BESSIState, runtime) = _host_vector(runtime.domain.runoff; copy_array=true)
_model_runoff_vector(::PDDModel, ::PDDState, runtime) = _host_vector(runtime.runoff; copy_array=true)
_model_pdd_vector(::AbstractSnowModel, ::AbstractSnowModelState, runtime, ncol::Int) = zeros(Float64, ncol)
_model_pdd_vector(::PDDModel, ::PDDState, runtime, ncol::Int) = _host_vector(runtime.pdd_sum; copy_array=true)
_model_melt_vector(::AbstractSnowModel, ::AbstractSnowModelState, runtime, ncol::Int) = zeros(Float64, ncol)
_model_melt_vector(::BESSIModel, ::BESSIState, runtime, ncol::Int) = _host_vector(runtime.domain.melt; copy_array=true)
_model_refreezing_vector(::AbstractSnowModel, ::AbstractSnowModelState, runtime, ncol::Int) = zeros(Float64, ncol)
_model_refreezing_vector(::BESSIModel, ::BESSIState, runtime, ncol::Int) = _host_vector(runtime.domain.refreezing; copy_array=true)

function _update_year_smb_delta!(last_delta::Vector{Float64}, previous_year_smb_ice::Vector{Float64}, model::AbstractSnowModel, state::AbstractSnowModelState, runtime)
    current = _model_smb_ice_vector(model, state, runtime)
    last_delta .= current .- previous_year_smb_ice
    previous_year_smb_ice .= current
    return nothing
end

function init_diagnostics!(sim, model_runtime::ModelRuntime, options::RunOptions, timings::StepTimingStats)
    model = sim.model
    state = sim.now
    runtime = model_runtime.backend
    ncol = model_runtime.ncol
    prev = allocate_year_summary_buffers(ncol)
    final = allocate_year_summary_buffers(ncol)
    backend_year_summary = _allocate_year_backend_buffers(model, state, runtime, ncol)
    time_block!(timings, :summarize_columns_initial) do
        summarize_year_state!(backend_year_summary, model, state, runtime)
        _copy_summary_fields!(prev, backend_year_summary, YEAR_BUFFER_NAMES)
    end

    selected = Set(options.netcdf_variables)
    need_step_outputs = options.write_netcdf && any(var -> var in selected, OUTPUT_GROUPS.step)
    need_monthly_outputs = options.write_netcdf && any(var -> var in selected, OUTPUT_GROUPS.monthly)
    need_layer_outputs = options.write_netcdf && any(var -> var in selected, OUTPUT_GROUPS.layers)
    need_last_year_smb_delta = options.write_netcdf && (:last_year_delta_ice_sheet_smb in selected)
    need_step_diagnostics = need_step_outputs || need_monthly_outputs

    step_summary = need_step_diagnostics ? allocate_summary_buffers(ncol) : nothing
    backend_step_summary = need_step_diagnostics ? _allocate_step_backend_buffers(model, state, runtime, ncol) : nothing
    deltas = (
        thickness=fill(NaN, ncol),
        wet_mass=fill(NaN, ncol),
        base_mass=fill(NaN, ncol),
        ice_sheet_smb=need_last_year_smb_delta ? fill(NaN, ncol) : Float64[],
    )
    previous = (
        base_mass=_host_vector(prev.base_mass; copy_array=true),
        wet_mass=_host_vector(prev.wet_mass; copy_array=true),
        smb_ice=_model_smb_ice_vector(model, state, runtime),
        runoff=_model_runoff_vector(model, state, runtime),
        pdd=_model_pdd_vector(model, state, runtime, ncol),
        melt=_model_melt_vector(model, state, runtime, ncol),
        refreezing=_model_refreezing_vector(model, state, runtime, ncol),
    )
    previous_year_smb_ice = need_last_year_smb_delta ? _model_smb_ice_vector(model, state, runtime) : Float64[]

    return DiagnosticsRuntime(
        need_step_outputs,
        need_monthly_outputs,
        need_layer_outputs,
        need_last_year_smb_delta,
        need_step_diagnostics,
        prev,
        final,
        backend_year_summary,
        step_summary,
        backend_step_summary,
        deltas,
        previous,
        previous_year_smb_ice,
        NamedTuple[],
    )
end

function accumulate_step_diagnostics!(integrator::SimulationIntegrator)
    diagnostics = integrator.diagnostics
    diagnostics.need_step_diagnostics || return nothing

    model = integrator.sim.model
    state = integrator.sim.now
    model_runtime = integrator.model_runtime
    output = integrator.output
    runtime = model_runtime.backend
    month_idx = diagnostics.need_monthly_outputs ?
        integrator.completed_years * output.schedule.nmonth_per_year + output.schedule.step_month[integrator.time_index] :
        0
    time_counted_block!(integrator.timings, :step_diagnostics, model_runtime.ncol) do
        summarize_step_state!(diagnostics.backend_step_summary, model, state, runtime)
        _copy_summary_fields!(diagnostics.step_summary, diagnostics.backend_step_summary, SUMMARY_BUFFER_NAMES)
        _accumulate_step_diagnostics!(
            diagnostics.step_summary,
            diagnostics.previous,
            output.monthly_sums,
            output.step_vectors,
            month_idx,
            diagnostics.need_monthly_outputs,
            diagnostics.need_step_outputs,
        )
    end
    diagnostics.need_monthly_outputs && (output.monthly_count[month_idx] += 1)
    return nothing
end

function _complete_year!(integrator::SimulationIntegrator)
    model = integrator.sim.model
    state = integrator.sim.now
    model_runtime = integrator.model_runtime
    diagnostics = integrator.diagnostics
    runtime = model_runtime.backend
    year = integrator.completed_years + 1

    time_block!(integrator.timings, :summarize_columns_year) do
        summarize_year_state!(diagnostics.backend_year_summary, model, state, runtime)
        _copy_summary_fields!(diagnostics.final, diagnostics.backend_year_summary, YEAR_BUFFER_NAMES)
    end

    if should_record_year_metrics(year, integrator.options.years, integrator.options.history_year_stride)
        record = time_block!(integrator.timings, :year_metrics) do
            make_year_record_and_deltas!(
                year,
                diagnostics.deltas.thickness,
                diagnostics.deltas.wet_mass,
                diagnostics.deltas.base_mass,
                diagnostics.final.thickness,
                diagnostics.final.wet_mass,
                diagnostics.final.bulk_density,
                diagnostics.final.base_mass,
                diagnostics.prev.thickness,
                diagnostics.prev.wet_mass,
                diagnostics.prev.base_mass,
            )
        end
        push!(diagnostics.history, record)
        time_block!(integrator.timings, :year_logging) do
            println(integrator.io, year_log_line(record))
        end
    end

    if diagnostics.need_last_year_smb_delta
        time_block!(integrator.timings, :year_state_deltas) do
            _update_year_smb_delta!(
                diagnostics.deltas.ice_sheet_smb,
                diagnostics.previous_year_smb_ice,
                model,
                state,
                runtime,
            )
        end
    end

    diagnostics.prev, diagnostics.final = diagnostics.final, diagnostics.prev
    integrator.completed_years = year
    next!(integrator.progress)
    return nothing
end

function _current_year_summary!(integrator::SimulationIntegrator)
    model = integrator.sim.model
    state = integrator.sim.now
    runtime = integrator.model_runtime.backend
    diagnostics = integrator.diagnostics
    time_block!(integrator.timings, :summarize_columns_final) do
        summarize_year_state!(diagnostics.backend_year_summary, model, state, runtime)
        _copy_summary_fields!(diagnostics.final, diagnostics.backend_year_summary, YEAR_BUFFER_NAMES)
    end
    return diagnostics.final
end
