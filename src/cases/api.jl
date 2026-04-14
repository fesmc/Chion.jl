function SnowpackCase(definition::CaseDefinition; run::RunConfig=RunConfig(), kwargs...)
    resolved_run = isempty(kwargs) ? run : _resolved_run_config(definition, run; kwargs...)
    resolved_run.write_netcdf && isnothing(definition.layout) &&
        error("`write_netcdf=true` requires a grid layout, but this case definition has `layout=nothing`.")
    return SnowpackCase(resolved_run.name, definition, resolved_run)
end

function SnowpackCase(
    domain::SM.SnowpackDomain,
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

function _resolved_run_config(
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

    resolved_name = isnothing(name) ? run.name : String(name)
    resolved_input_label = isnothing(input_label) ? run.input_label : String(input_label)
    resolved_output_dir = isnothing(output_dir) ? out_dir : output_dir
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

    final_name = isempty(resolved_name) ? "snowpack_case" : resolved_name
    final_output_dir = let candidate = isnothing(resolved_output_dir) ? run.output_dir : String(resolved_output_dir)
        isempty(candidate) ? _default_case_output_dir(final_name) : candidate
    end
    final_input_label = isempty(resolved_input_label) ? definition.input_label : resolved_input_label

    return RunConfig(
        name=final_name,
        input_label=final_input_label,
        output_dir=final_output_dir,
        netcdf_path=isnothing(netcdf_path) ? run.netcdf_path : String(netcdf_path),
        write_outputs=isnothing(write_outputs) ? run.write_outputs : write_outputs,
        write_netcdf=isnothing(write_netcdf) ? run.write_netcdf : write_netcdf,
        netcdf_variables=isnothing(netcdf_variables) ? run.netcdf_variables : netcdf_variables,
        cycles=resolved_cycles,
        backend=isnothing(backend) ? run.backend : backend,
        history_stride=resolved_history_stride,
    )
end

"""
    build_case(definition; run=RunConfig(), kwargs...)

Create a runnable [`SnowpackCase`](@ref) from a reusable [`CaseDefinition`](@ref).
Pass `run=RunConfig(...)` for the baseline configuration and optional keyword
overrides when you want to tweak a few fields.
"""
function build_case(
    definition::CaseDefinition;
    run::RunConfig=RunConfig(),
    kwargs...,
)
    return SnowpackCase(definition; run=run, kwargs...)
end

function build_case(
    domain::SM.SnowpackDomain,
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

"""
    synthetic_case(; ...)

Create a reusable synthetic [`CaseDefinition`](@ref).
"""
function synthetic_case(;
    physics::SM.SnowpackPhysicalConstants{Float64}=physics(),
    ntot::Integer=5,
    variant::Symbol=:multi_column,
    ntime::Integer=12,
    nx::Union{Nothing, Integer}=nothing,
    ny::Union{Nothing, Integer}=nothing,
)
    return synthetic_definition(
        physics=physics,
        ntot=ntot,
        variant=variant,
        ntime=ntime,
        nx=nx,
        ny=ny,
    )
end

"""
    prescribed_case(; ...)

Create a reusable [`CaseDefinition`](@ref) directly from user-facing forcing
vectors and simple initial-condition keywords.
"""
function prescribed_case(;
    physics::SM.SnowpackPhysicalConstants{Float64}=physics(),
    ntot::Integer=5,
    nx::Integer=1,
    ny::Integer=1,
    dt_days,
    air_temperature_c,
    snowfall_mm_day,
    rainfall_mm_day,
    shortwave_down,
    wind_speed=5.0,
    q_lw_down=nothing,
    has_q_lw_down=nothing,
    q_sh=nothing,
    has_q_sh=nothing,
    q_lh=nothing,
    has_q_lh=nothing,
    time_values=nothing,
    initial_surface_mass=250.0,
    initial_density=320.0,
    initial_temperature_c=-12.0,
    initial_albedo=physics.alpha_dry,
    input_label::AbstractString="prescribed_forcing",
)
    return prescribed_definition(
        physics=physics,
        ntot=ntot,
        nx=nx,
        ny=ny,
        dt_days=dt_days,
        air_temperature_c=air_temperature_c,
        snowfall_mm_day=snowfall_mm_day,
        rainfall_mm_day=rainfall_mm_day,
        shortwave_down=shortwave_down,
        wind_speed=wind_speed,
        q_lw_down=q_lw_down,
        has_q_lw_down=has_q_lw_down,
        q_sh=q_sh,
        has_q_sh=has_q_sh,
        q_lh=q_lh,
        has_q_lh=has_q_lh,
        time_values=time_values,
        initial_surface_mass=initial_surface_mass,
        initial_density=initial_density,
        initial_temperature_c=initial_temperature_c,
        initial_albedo=initial_albedo,
        input_label=input_label,
    )
end

"""
    mar_case(path; ...)

Create a reusable [`CaseDefinition`](@ref) directly from a MAR NetCDF/HDF5 file.
"""
function mar_case(
    path::AbstractString;
    physics::SM.SnowpackPhysicalConstants{Float64}=physics(),
    ntot::Integer=20,
    mask_threshold::Real=50.0,
    turbulent_flux_sign::Real=1.0,
)
    return mar_definition(
        path;
        physics=physics,
        ntot=ntot,
        mask_threshold=mask_threshold,
        turbulent_flux_sign=turbulent_flux_sign,
    )
end

"""
    run_case(case; io=stdout, copy_domain=true, timings=TimingStats())

Run a high-level [`SnowpackCase`](@ref). By default the initial domain is
deep-copied so the same case can be rerun on CPU and GPU without reparsing the
forcing. Set `copy_domain=false` to reuse and mutate the stored domain.
"""
function run_case(
    case::SnowpackCase;
    io::IO=stdout,
    copy_domain::Bool=true,
    timings::TimingStats=TimingStats(),
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
