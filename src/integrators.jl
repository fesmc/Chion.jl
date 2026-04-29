"""Initialized runtime and stepping internals for `Simulation`."""

_model_grid(model::AbstractSnowModel) = model.grid
_model_column_count(model::AbstractSnowModel) = ncols(_model_grid(model))
model_layer_count(::AbstractSnowModel, ::AbstractSnowModelState, runtime) = 0
model_layer_count(::BESSIModel, ::BESSIState, runtime) = runtime.domain.Ntot

model_output_groups(::BESSIModel) = (:final, :layers, :history, :monthly, :step)
model_output_groups(::PDDModel) = (:final, :history, :step)
model_output_groups(::ITMModel) = ()

function _validate_model_outputs!(model::AbstractSnowModel, options::RunOptions)
    options.write_netcdf || return nothing
    supported = model_output_groups(model)
    allowed = Symbol[]
    for group in supported
        append!(allowed, getproperty(OUTPUT_GROUPS, group))
    end
    unsupported = setdiff(options.netcdf_variables, allowed)
    isempty(unsupported) || error("$(typeof(model)) output supports only $(join(string.(supported), ", ")) NetCDF groups; unsupported: $(join(string.(unsupported), ", ")).")
    return nothing
end

function _validate_integrator_setup!(sim, options::RunOptions)
    model = sim.model
    state = sim.now
    forcing = sim.forcing
    grid = _model_grid(model)
    ncol = _model_column_count(model)
    size(forcing.air_temperature, 1) == ncol || error("Forcing column count must match the model column count.")
    spatial_grid = has_spatial_coords(grid)
    options.write_netcdf && !spatial_grid && error("NetCDF output requires a grid with spatial coordinates.")
    spatial_grid && length(grid.js) != ncol && error("Grid point count must match the domain column count.")
    _validate_model_outputs!(model, options)
    return (model=model, state=state, forcing=forcing, grid=grid, ncol=ncol)
end

function _prepare_backend!(timings::StepTimingStats, domain::SnowpackDomain, forcing::SnowpackForcing; is_gpu::Bool)
    step_fields = forcing
    if is_gpu
        cuda_available() || error("`backend=gpu` requested, but CUDA is not functional in the current environment.")
        domain = time_block!(timings, :gpu_transfer) do
            gpu_domain(domain)
        end
        step_fields = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), forcing)
        end
        workspace = time_block!(timings, :gpu_transfer) do
            ColumnarStepWorkspace(domain)
        end
        return (domain=domain, step_fields=step_fields, workspace=workspace, is_gpu=true)
    end
    workspace = time_block!(timings, :create_workspaces) do
        ColumnarStepWorkspace(domain)
    end
    return (domain=domain, step_fields=step_fields, workspace=workspace, is_gpu=false)
end

function prepare_runtime!(model::BESSIModel, state::BESSIState, forcing::SnowpackForcing, options::RunOptions, timings::StepTimingStats)
    return _prepare_backend!(timings, state.domain, forcing; is_gpu=options.backend == :gpu)
end

function prepare_runtime!(model::PDDModel, state::PDDState, forcing::SnowpackForcing, options::RunOptions, timings::StepTimingStats)
    is_gpu = options.backend == :gpu
    if is_gpu
        cuda_available() || error("`backend=gpu` requested, but CUDA is not functional in the current environment.")
        snowpack_swe = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), state.snowpack_swe)
        end
        smb_ice = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), state.smb_ice)
        end
        runoff = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), state.runoff)
        end
        pdd_sum = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), state.pdd_sum)
        end
        step_fields = time_block!(timings, :gpu_transfer) do
            adapt(gpu_storage_type(), forcing)
        end
        return (snowpack_swe=snowpack_swe, smb_ice=smb_ice, runoff=runoff, pdd_sum=pdd_sum, step_fields=step_fields, is_gpu=true)
    end
    ncol = ncols(model.grid)
    scratch = (
        a=Vector{Float64}(undef, ncol),
        b=Vector{Float64}(undef, ncol),
        c=Vector{Float64}(undef, ncol),
        d=Vector{Float64}(undef, ncol),
        e=Vector{Float64}(undef, ncol),
        f=Vector{Float64}(undef, ncol),
    )
    return (snowpack_swe=state.snowpack_swe, smb_ice=state.smb_ice, runoff=state.runoff, pdd_sum=state.pdd_sum, step_fields=forcing, scratch=scratch, is_gpu=false)
end

function prepare_runtime!(::ITMModel, ::AbstractSnowModelState, ::SnowpackForcing, ::RunOptions, ::StepTimingStats)
    error("ITMModel is not yet implemented. Physics coming soon.")
end

const SUMMARY_BUFFER_NAMES = (:thickness, :wet_mass, :bulk_density, :base_mass, :smb_ice, :liquid_water, :runoff, :pdd)
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
    return summary
end

function step_model!(::BESSIModel, ::BESSIState, runtime, forcing::SnowpackForcing, time_index::Int)
    step!(runtime.domain, forcing, time_index, runtime.workspace)
    return nothing
end

function step_model!(model::PDDModel, ::PDDState, runtime, forcing::SnowpackForcing, time_index::Int)
    if runtime.snowpack_swe isa Vector{Float64}
        pdd_step!(
            runtime.snowpack_swe,
            runtime.smb_ice,
            runtime.runoff,
            runtime.pdd_sum,
            forcing,
            time_index,
            model.ddf_snow,
            model.ddf_ice,
            model.refreezing_fraction,
            runtime.scratch,
        )
    else
        pdd_step!(
            runtime.snowpack_swe,
            runtime.smb_ice,
            runtime.runoff,
            runtime.pdd_sum,
            forcing,
            time_index,
            model.ddf_snow,
            model.ddf_ice,
            model.refreezing_fraction,
        )
    end
    return nothing
end

_model_smb_ice_vector(::BESSIModel, ::BESSIState, runtime) = _host_vector(runtime.domain.smb_ice; copy_array=true)
_model_smb_ice_vector(::PDDModel, ::PDDState, runtime) = _host_vector(runtime.smb_ice; copy_array=true)
_model_runoff_vector(::BESSIModel, ::BESSIState, runtime) = _host_vector(runtime.domain.runoff; copy_array=true)
_model_runoff_vector(::PDDModel, ::PDDState, runtime) = _host_vector(runtime.runoff; copy_array=true)
_model_pdd_vector(::AbstractSnowModel, ::AbstractSnowModelState, runtime, ncol::Int) = zeros(Float64, ncol)
_model_pdd_vector(::PDDModel, ::PDDState, runtime, ncol::Int) = _host_vector(runtime.pdd_sum; copy_array=true)

function _update_year_smb_delta!(last_delta::Vector{Float64}, previous_year_smb_ice::Vector{Float64}, model::AbstractSnowModel, state::AbstractSnowModelState, runtime)
    current = _model_smb_ice_vector(model, state, runtime)
    last_delta .= current .- previous_year_smb_ice
    previous_year_smb_ice .= current
    return nothing
end

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

finalize_state!(::AbstractSnowModel, ::AbstractSnowModelState, runtime, ::RunOptions, ::StepTimingStats) = nothing

function finalize_state!(::BESSIModel, state::BESSIState, runtime, options::RunOptions, timings::StepTimingStats)
    runtime.is_gpu || return nothing
    time_block!(timings, :gpu_transfer) do
        _copy_domain_state!(state.domain, cpu_domain(runtime.domain))
    end
    return nothing
end

function finalize_state!(::PDDModel, state::PDDState, runtime, options::RunOptions, timings::StepTimingStats)
    runtime.is_gpu || return nothing
    time_block!(timings, :gpu_transfer) do
        copyto!(state.snowpack_swe, Array(runtime.snowpack_swe))
        copyto!(state.smb_ice, Array(runtime.smb_ice))
        copyto!(state.runoff, Array(runtime.runoff))
        copyto!(state.pdd_sum, Array(runtime.pdd_sum))
    end
    return nothing
end

struct IntegratorClocks
    run_wall_t0::Int
    simulation_wall_t0::Int
end

struct IntegratorModelRuntime
    runtime
    ncol::Int
    grid
end

mutable struct DiagnosticsRuntime
    selected::Set{Symbol}
    need_step_outputs::Bool
    need_monthly_outputs::Bool
    need_layer_outputs::Bool
    need_last_year_smb_delta::Bool
    need_step_diagnostics::Bool
    prev
    final
    backend_year_summary
    step_summary
    backend_step_summary
    deltas
    previous
    previous_year_smb_ice::Vector{Float64}
    history::Vector{NamedTuple}
end

mutable struct OutputRuntime
    schedule
    writer
    nc_path::String
    monthly_sums
    monthly_count::Vector{Int32}
    step_vectors
    steps_written::Int
end

mutable struct StepperRuntime
    time_index::Int
    completed_years::Int
    current_forcing::SnowpackForcing
    progress
end

mutable struct SimulationIntegrator
    sim
    options::RunOptions
    io::IO
    timings::StepTimingStats
    clocks::IntegratorClocks
    model_runtime::IntegratorModelRuntime
    diagnostics::DiagnosticsRuntime
    output::OutputRuntime
    stepper::StepperRuntime
    finalized::Bool
    result::Union{Nothing, SimulationResult}
end

function init_model_runtime!(context, options::RunOptions, timings::StepTimingStats)
    runtime = prepare_runtime!(context.model, context.state, context.forcing, options, timings)
    return IntegratorModelRuntime(runtime, context.ncol, context.grid)
end

function init_diagnostics!(context, model_runtime::IntegratorModelRuntime, options::RunOptions, timings::StepTimingStats)
    model = context.model
    state = context.state
    runtime = model_runtime.runtime
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
        smb_ice=_model_smb_ice_vector(model, state, runtime),
        runoff=_model_runoff_vector(model, state, runtime),
        pdd=_model_pdd_vector(model, state, runtime, ncol),
    )
    previous_year_smb_ice = need_last_year_smb_delta ? _model_smb_ice_vector(model, state, runtime) : Float64[]

    return DiagnosticsRuntime(
        selected,
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

function _single_step_forcing_template(forcing::SnowpackForcing)
    return SnowpackForcing(
        dt_days=[first(forcing.dt_days)],
        ncol=size(forcing.air_temperature, 1),
        air_temperature=forcing.air_temperature[:, 1:1],
        snowfall_rate=forcing.snowfall_rate[:, 1:1],
        rainfall_rate=forcing.rainfall_rate[:, 1:1],
        shortwave_down=forcing.shortwave_down[:, 1:1],
        wind_speed=forcing.wind_speed[:, 1:1],
        q_lw_down=forcing.q_lw_down[:, 1:1],
        has_q_lw_down=forcing.has_q_lw_down[:, 1:1],
        q_sh=forcing.q_sh[:, 1:1],
        has_q_sh=forcing.has_q_sh[:, 1:1],
        q_lh=forcing.q_lh[:, 1:1],
        has_q_lh=forcing.has_q_lh[:, 1:1],
        time_values=[first(forcing.time_values)],
    )
end

function init_stepper_state!(context, options::RunOptions, io::IO)
    progress = Progress(options.years; desc="Running years: ", output=io, showspeed=true)
    return StepperRuntime(1, 0, _single_step_forcing_template(context.forcing), progress)
end

function _new_integrator(sim, options::RunOptions, io::IO, timings::StepTimingStats, model_runtime, diagnostics, output, stepper)
    return SimulationIntegrator(
        sim,
        options,
        io,
        timings,
        IntegratorClocks(time_ns(), time_ns()),
        model_runtime,
        diagnostics,
        output,
        stepper,
        false,
        nothing,
    )
end

_finished(integrator::SimulationIntegrator) =
    integrator.stepper.completed_years >= integrator.options.years

function _external_forcing_for_runtime(integrator::SimulationIntegrator)
    forcing = integrator.stepper.current_forcing
    integrator.options.backend == :gpu || return forcing
    return time_block!(integrator.timings, :gpu_transfer) do
        adapt(gpu_storage_type(), forcing)
    end
end

_scheduled_forcing_for_runtime(integrator::SimulationIntegrator) =
    integrator.model_runtime.runtime.step_fields

function _advance_with_forcing!(integrator::SimulationIntegrator, forcing::SnowpackForcing, time_index::Int)
    integrator.finalized && error("Cannot step a finalized integrator.")
    _finished(integrator) && error("Cannot step an integrator that has already completed all years.")

    model = integrator.sim.model
    state = integrator.sim.now
    model_runtime = integrator.model_runtime
    diagnostics = integrator.diagnostics
    output = integrator.output
    stepper = integrator.stepper
    runtime = model_runtime.runtime

    time_counted_block!(integrator.timings, :model_step_wall, model_runtime.ncol) do
        step_model!(model, state, runtime, forcing, time_index)
    end

    if diagnostics.need_step_diagnostics
        month_idx = diagnostics.need_monthly_outputs ?
            stepper.completed_years * output.schedule.nmonth_per_year + output.schedule.step_month[stepper.time_index] :
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
    end

    if diagnostics.need_step_outputs && output.schedule.annual_output.write_output[stepper.time_index]
        output.steps_written += 1
        if output.writer !== nothing
            step_grids = time_block!(integrator.timings, :step_output_prepare) do
                _step_output_grids(output.step_vectors, model_runtime.grid)
            end
            time_block!(integrator.timings, :step_output_write) do
                for key in OUTPUT_GROUPS.step
                    maybe_write_step_output!(output.writer, output.steps_written, key, getfield(step_grids, key))
                end
            end
            _reset_step_vectors!(output.step_vectors)
        end
    end

    if stepper.time_index == length(integrator.sim.forcing.time_values)
        _complete_year!(integrator)
        stepper.time_index = 1
    else
        stepper.time_index += 1
    end
    return nothing
end

function _complete_year!(integrator::SimulationIntegrator)
    model = integrator.sim.model
    state = integrator.sim.now
    model_runtime = integrator.model_runtime
    diagnostics = integrator.diagnostics
    stepper = integrator.stepper
    runtime = model_runtime.runtime
    year = stepper.completed_years + 1

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
    stepper.completed_years = year
    next!(stepper.progress)
    return nothing
end

function _step_scheduled!(integrator::SimulationIntegrator)
    forcing = _scheduled_forcing_for_runtime(integrator)
    return _advance_with_forcing!(integrator, forcing, integrator.stepper.time_index)
end

function _step_n!(integrator::SimulationIntegrator, n::Integer)
    n >= 0 || error("Step count must be non-negative.")
    for _ in 1:Int(n)
        _step_scheduled!(integrator)
    end
    return nothing
end

function _step_external!(integrator::SimulationIntegrator, Δt_days::Real, force_dt::Bool=true)
    force_dt || error("Chion's initialized stepper requires `force_dt=true`, matching the FastIsostasy coupling pattern.")
    Δt_days > 0 || error("`Δt_days` must be positive.")
    integrator.stepper.current_forcing.dt_days[1] = Float64(Δt_days)
    forcing = _external_forcing_for_runtime(integrator)
    return _advance_with_forcing!(integrator, forcing, 1)
end

function _run_integrator!(integrator::SimulationIntegrator)
    while !_finished(integrator)
        _step_scheduled!(integrator)
    end
    return nothing
end

function _assign_numeric_step_field!(dest::AbstractMatrix, value, ncol::Int, name::AbstractString; transform=identity)
    isnothing(value) && return false
    dest[:, :] .= transform.(_forcing_numeric_matrix(value, ncol, 1, name))
    return true
end

function _assign_bool_step_field!(dest::AbstractMatrix{Bool}, value, ncol::Int, name::AbstractString)
    isnothing(value) && return false
    dest[:, :] .= _forcing_bool_matrix(value, ncol, 1, name)
    return true
end

function _set_forcing!(
    integrator::SimulationIntegrator;
    air_temperature=nothing,
    snowfall_rate=nothing,
    rainfall_rate=nothing,
    air_temperature_c=nothing,
    snowfall_mm_day=nothing,
    rainfall_mm_day=nothing,
    shortwave_down=nothing,
    wind_speed=nothing,
    q_lw_down=nothing,
    has_q_lw_down=nothing,
    q_sh=nothing,
    has_q_sh=nothing,
    q_lh=nothing,
    has_q_lh=nothing,
    time_value=nothing,
)
    f = integrator.stepper.current_forcing
    ncol = integrator.model_runtime.ncol
    has_native = !isnothing(air_temperature) || !isnothing(snowfall_rate) || !isnothing(rainfall_rate)
    has_user = !isnothing(air_temperature_c) || !isnothing(snowfall_mm_day) || !isnothing(rainfall_mm_day)
    has_native && has_user && error("Pass either model-native forcing fields or user-facing fields, not both.")
    if has_user
        _assign_numeric_step_field!(f.air_temperature, air_temperature_c, ncol, "air_temperature_c"; transform=x -> x + 273.15)
        _assign_numeric_step_field!(f.snowfall_rate, snowfall_mm_day, ncol, "snowfall_mm_day"; transform=x -> x / 86_400.0)
        _assign_numeric_step_field!(f.rainfall_rate, rainfall_mm_day, ncol, "rainfall_mm_day"; transform=x -> x / 86_400.0)
    else
        _assign_numeric_step_field!(f.air_temperature, air_temperature, ncol, "air_temperature")
        _assign_numeric_step_field!(f.snowfall_rate, snowfall_rate, ncol, "snowfall_rate")
        _assign_numeric_step_field!(f.rainfall_rate, rainfall_rate, ncol, "rainfall_rate")
    end
    _assign_numeric_step_field!(f.shortwave_down, shortwave_down, ncol, "shortwave_down")
    _assign_numeric_step_field!(f.wind_speed, wind_speed, ncol, "wind_speed")
    if _assign_numeric_step_field!(f.q_lw_down, q_lw_down, ncol, "q_lw_down") && isnothing(has_q_lw_down)
        fill!(f.has_q_lw_down, true)
    end
    _assign_bool_step_field!(f.has_q_lw_down, has_q_lw_down, ncol, "has_q_lw_down")
    if _assign_numeric_step_field!(f.q_sh, q_sh, ncol, "q_sh") && isnothing(has_q_sh)
        fill!(f.has_q_sh, true)
    end
    _assign_bool_step_field!(f.has_q_sh, has_q_sh, ncol, "has_q_sh")
    if _assign_numeric_step_field!(f.q_lh, q_lh, ncol, "q_lh") && isnothing(has_q_lh)
        fill!(f.has_q_lh, true)
    end
    _assign_bool_step_field!(f.has_q_lh, has_q_lh, ncol, "has_q_lh")
    isnothing(time_value) || (f.time_values[1] = DateTime(time_value))
    return integrator
end

function _current_year_summary!(integrator::SimulationIntegrator)
    model = integrator.sim.model
    state = integrator.sim.now
    runtime = integrator.model_runtime.runtime
    diagnostics = integrator.diagnostics
    time_block!(integrator.timings, :summarize_columns_final) do
        summarize_year_state!(diagnostics.backend_year_summary, model, state, runtime)
        _copy_summary_fields!(diagnostics.final, diagnostics.backend_year_summary, YEAR_BUFFER_NAMES)
    end
    return diagnostics.final
end

function _finalize_integrator!(integrator::SimulationIntegrator)
    integrator.finalized && return integrator.result

    model_runtime = integrator.model_runtime
    diagnostics = integrator.diagnostics
    output = integrator.output
    history = diagnostics.history
    simulation_wall_sec = (time_ns() - integrator.clocks.simulation_wall_t0) * 1.0e-9
    status = _finished(integrator) ? :complete : :incomplete
    years_completed = completed_year_count(history, status, integrator.options.years)
    summary_path = ""
    history_csv_path = ""

    if integrator.options.write_outputs
        mkpath(integrator.options.output_dir)
        summary_path = joinpath(integrator.options.output_dir, "$(integrator.options.name)_summary.txt")
        history_csv_path = joinpath(integrator.options.output_dir, "$(integrator.options.name)_history.csv")
        time_block!(integrator.timings, :write_summary_text) do
            write_run_summary(summary_path, integrator.options, integrator.sim.forcing.time_values, model_runtime.ncol, history, status, integrator.timings)
        end
        time_block!(integrator.timings, :write_history_csv) do
            write_run_history_csv(history_csv_path, history)
        end
    end

    if output.writer !== nothing
        final_state = _current_year_summary!(integrator)
        final_grids = time_block!(integrator.timings, :scatter_final_outputs) do
            _scatter_model_final_grids(integrator.sim.model, integrator.sim.now, model_runtime.runtime, final_state, diagnostics.deltas, model_runtime.grid)
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
    end

    finalize_state!(integrator.sim.model, integrator.sim.now, model_runtime.runtime, integrator.options, integrator.timings)
    run_wall_sec = (time_ns() - integrator.clocks.run_wall_t0) * 1.0e-9
    print_run_report(
        integrator.io,
        integrator.options,
        integrator.sim.forcing.time_values,
        history,
        status,
        simulation_wall_sec,
        run_wall_sec,
        integrator.timings;
        nc_path=output.nc_path,
        summary_path=summary_path,
        history_csv_path=history_csv_path,
    )
    result = SimulationResult(
        history,
        status,
        years_completed,
        integrator.timings,
        simulation_wall_sec,
        run_wall_sec,
        output.nc_path,
        summary_path,
        history_csv_path,
    )
    integrator.finalized = true
    integrator.result = result
    return result
end
