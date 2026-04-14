@inline function _resolved_case_name(name::AbstractString, default_name::AbstractString)
    stripped = strip(name)
    return isempty(stripped) || stripped == "snowpack_case" ? String(default_name) : String(stripped)
end

function _resolved_case_run_config(
    definition::CaseDefinition,
    run::RunConfig,
    default_name::AbstractString,
)
    resolved_name = _resolved_case_name(run.name, default_name)
    resolved_input_label = isempty(run.input_label) ? definition.input_label : run.input_label
    resolved_output_dir = isempty(run.output_dir) ? _default_case_output_dir(resolved_name) : run.output_dir
    return RunConfig(
        name=resolved_name,
        input_label=resolved_input_label,
        output_dir=resolved_output_dir,
        netcdf_path=run.netcdf_path,
        write_outputs=run.write_outputs,
        write_netcdf=run.write_netcdf,
        netcdf_variables=run.netcdf_variables,
        cycles=run.cycles,
        backend=run.backend,
        history_stride=run.history_stride,
    )
end

function _resolved_build_run_config(
    definition::CaseDefinition,
    run::RunConfig;
    name::Union{Nothing, AbstractString}=nothing,
    input_label::Union{Nothing, AbstractString}=nothing,
    output_dir::Union{Nothing, AbstractString}=nothing,
    out_dir::Union{Nothing, AbstractString}=nothing,
    netcdf_path::Union{Nothing, AbstractString}=nothing,
    write_outputs::Union{Nothing, Bool}=nothing,
    write_netcdf::Union{Nothing, Bool}=nothing,
    netcdf_variables=nothing,
    cycles::Union{Nothing, Integer}=nothing,
    max_cycles::Union{Nothing, Integer}=nothing,
    backend=nothing,
    history_stride::Union{Nothing, Integer}=nothing,
    cycle_metrics_stride::Union{Nothing, Integer}=nothing,
    run_forcing_once::Union{Nothing, Bool}=nothing,
)
    isnothing(output_dir) || isnothing(out_dir) || error("Pass only one of `output_dir` or `out_dir`.")
    isnothing(cycles) || isnothing(max_cycles) || error("Pass only one of `cycles` or `max_cycles`.")
    isnothing(history_stride) || isnothing(cycle_metrics_stride) || error("Pass only one of `history_stride` or `cycle_metrics_stride`.")

    resolved_cycles = if !isnothing(cycles)
        Int(cycles)
    elseif !isnothing(max_cycles)
        Int(max_cycles)
    elseif run_forcing_once === true
        1
    else
        run.cycles
    end
    resolved_history_stride = isnothing(history_stride) ?
        (isnothing(cycle_metrics_stride) ? run.history_stride : Int(cycle_metrics_stride)) :
        Int(history_stride)

    return RunConfig(
        name=isnothing(name) ? run.name : String(name),
        input_label=isnothing(input_label) ? (isempty(run.input_label) ? definition.input_label : run.input_label) : String(input_label),
        output_dir=isnothing(output_dir) ? (isnothing(out_dir) ? run.output_dir : String(out_dir)) : String(output_dir),
        netcdf_path=isnothing(netcdf_path) ? run.netcdf_path : String(netcdf_path),
        write_outputs=isnothing(write_outputs) ? run.write_outputs : write_outputs,
        write_netcdf=isnothing(write_netcdf) ? run.write_netcdf : write_netcdf,
        netcdf_variables=isnothing(netcdf_variables) ? run.netcdf_variables : netcdf_variables,
        cycles=resolved_cycles,
        backend=isnothing(backend) ? run.backend : backend,
        history_stride=resolved_history_stride,
    )
end

function SnowpackCase(
    definition::CaseDefinition;
    run::RunConfig=RunConfig(),
    default_name::AbstractString="snowpack_case",
)
    resolved_run = _resolved_case_run_config(definition, run, default_name)
    resolved_run.write_netcdf && isnothing(definition.layout) &&
        error("`write_netcdf=true` requires a grid layout, but this case definition has `layout=nothing`.")
    return SnowpackCase(resolved_run.name, definition, resolved_run)
end

function SnowpackCase(
    domain::SnowpackDomain,
    forcing::ForcingData;
    layout::Union{Nothing, GridLayout}=nothing,
    input_label::AbstractString="",
    notes::AbstractVector{<:AbstractString}=String[],
    metadata::NamedTuple=(;),
    run::RunConfig=RunConfig(),
    default_name::AbstractString="snowpack_case",
)
    definition = CaseDefinition(
        domain,
        forcing;
        layout=layout,
        input_label=input_label,
        notes=notes,
        metadata=metadata,
    )
    return SnowpackCase(definition; run=run, default_name=default_name)
end

function build_case(
    definition::CaseDefinition;
    run::RunConfig=RunConfig(),
    kwargs...,
)
    resolved_run = isempty(kwargs) ? run : _resolved_build_run_config(definition, run; kwargs...)
    return SnowpackCase(definition; run=resolved_run, default_name=resolved_run.name)
end

function build_case(
    domain::SnowpackDomain,
    forcing::ForcingData;
    layout::Union{Nothing, GridLayout}=nothing,
    input_label::AbstractString="",
    notes::AbstractVector{<:AbstractString}=String[],
    metadata::NamedTuple=(;),
    kwargs...,
)
    definition = CaseDefinition(
        domain,
        forcing;
        layout=layout,
        input_label=input_label,
        notes=notes,
        metadata=metadata,
    )
    return build_case(definition; kwargs...)
end

@inline _reject_file_mode_kw(name::AbstractString, value) =
    isnothing(value) || error("`$name` is not supported when `forcing_file` is provided. Preprocess the file externally, then call `prescribed_case(; forcing_file=..., run=RunConfig(...))`.")

@inline function _resolved_prescribed_ntot(ntot, forcing_file)
    return isnothing(ntot) ? (isnothing(forcing_file) ? 5 : 20) : Int(ntot)
end

"""
    synthetic_case(; run=RunConfig(), ...)

Create a runnable synthetic `SnowpackCase`.
"""
function synthetic_case(;
    run::RunConfig=RunConfig(),
    physics::SnowpackPhysicalConstants{Float64}=physics(),
    ntot::Integer=5,
    variant::Symbol=:multi_column,
    ntime::Integer=12,
    nx::Union{Nothing, Integer}=nothing,
    ny::Union{Nothing, Integer}=nothing,
)
    definition = synthetic_definition(
        physics=physics,
        ntot=ntot,
        variant=variant,
        ntime=ntime,
        nx=nx,
        ny=ny,
    )
    return SnowpackCase(definition; run=run, default_name="synthetic_case")
end

"""
    prescribed_case(; run=RunConfig(), forcing_file=nothing, ...)

Create a runnable `SnowpackCase` from either direct prescribed forcing
arrays or a prepared external forcing file.

When `forcing_file` is provided, Chion assumes the file has already been masked
or subsetted externally. Cells are included only when the required prescribed
forcing fields remain finite across the loaded timeseries.
"""
function prescribed_case(;
    run::RunConfig=RunConfig(),
    forcing_file::Union{Nothing, AbstractString}=nothing,
    physics::SnowpackPhysicalConstants{Float64}=physics(),
    ntot::Union{Nothing, Integer}=nothing,
    nx::Union{Nothing, Integer}=nothing,
    ny::Union{Nothing, Integer}=nothing,
    dt_days=nothing,
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
    time_values=nothing,
    initial_surface_mass=nothing,
    initial_density=nothing,
    initial_temperature_c=nothing,
    initial_albedo=nothing,
    input_label::Union{Nothing, AbstractString}=nothing,
)
    resolved_ntot = _resolved_prescribed_ntot(ntot, forcing_file)
    definition = if isnothing(forcing_file)
        isnothing(dt_days) && error("`dt_days` is required unless `forcing_file` is provided.")
        isnothing(air_temperature_c) && error("`air_temperature_c` is required unless `forcing_file` is provided.")
        isnothing(snowfall_mm_day) && error("`snowfall_mm_day` is required unless `forcing_file` is provided.")
        isnothing(rainfall_mm_day) && error("`rainfall_mm_day` is required unless `forcing_file` is provided.")
        isnothing(shortwave_down) && error("`shortwave_down` is required unless `forcing_file` is provided.")
        prescribed_definition(
            physics=physics,
            ntot=resolved_ntot,
            nx=isnothing(nx) ? 1 : Int(nx),
            ny=isnothing(ny) ? 1 : Int(ny),
            dt_days=dt_days,
            air_temperature_c=air_temperature_c,
            snowfall_mm_day=snowfall_mm_day,
            rainfall_mm_day=rainfall_mm_day,
            shortwave_down=shortwave_down,
            wind_speed=isnothing(wind_speed) ? 5.0 : wind_speed,
            q_lw_down=q_lw_down,
            has_q_lw_down=has_q_lw_down,
            q_sh=q_sh,
            has_q_sh=has_q_sh,
            q_lh=q_lh,
            has_q_lh=has_q_lh,
            time_values=time_values,
            initial_surface_mass=isnothing(initial_surface_mass) ? 250.0 : initial_surface_mass,
            initial_density=isnothing(initial_density) ? 320.0 : initial_density,
            initial_temperature_c=isnothing(initial_temperature_c) ? -12.0 : initial_temperature_c,
            initial_albedo=isnothing(initial_albedo) ? physics.alpha_dry : initial_albedo,
            input_label=isnothing(input_label) ? "prescribed_forcing" : String(input_label),
        )
    else
        for (name, value) in (
            ("nx", nx),
            ("ny", ny),
            ("dt_days", dt_days),
            ("air_temperature_c", air_temperature_c),
            ("snowfall_mm_day", snowfall_mm_day),
            ("rainfall_mm_day", rainfall_mm_day),
            ("shortwave_down", shortwave_down),
            ("wind_speed", wind_speed),
            ("q_lw_down", q_lw_down),
            ("has_q_lw_down", has_q_lw_down),
            ("q_sh", q_sh),
            ("has_q_sh", has_q_sh),
            ("q_lh", q_lh),
            ("has_q_lh", has_q_lh),
            ("time_values", time_values),
            ("initial_surface_mass", initial_surface_mass),
            ("initial_density", initial_density),
            ("initial_temperature_c", initial_temperature_c),
            ("initial_albedo", initial_albedo),
            ("input_label", input_label),
        )
            _reject_file_mode_kw(name, value)
        end
        _prescribed_definition_from_forcing_file(
            forcing_file;
            physics=physics,
            ntot=resolved_ntot,
        )
    end
    return SnowpackCase(definition; run=run, default_name="prescribed_case")
end

"""
    run_case(case; io=stdout, copy_domain=true, timings=StepTimingStats())

Run a high-level `SnowpackCase`. By default the initial domain is
deep-copied so the same case can be rerun on CPU and GPU without reparsing the
forcing. Set `copy_domain=false` to reuse and mutate the stored domain.
"""
function run_case(
    case::SnowpackCase;
    io::IO=stdout,
    copy_domain::Bool=true,
    timings::StepTimingStats=StepTimingStats(),
    run_wall_t0::Integer=time_ns(),
)
    for note in case.definition.notes
        println(io, note)
    end
    domain = copy_domain ? deepcopy(case.definition.domain) : case.definition.domain
    return execute_case!(
        domain,
        case.definition.forcing;
        layout=case.definition.layout,
        options=case.run,
        io=io,
        timings=timings,
        run_wall_t0=run_wall_t0,
    )
end
