function SnowpackCase(
    definition::CaseDefinition;
    name::AbstractString="snowpack_case",
    input_label::Union{Nothing, AbstractString}=nothing,
    output_dir::AbstractString="",
    netcdf_path::AbstractString="",
    write_outputs::Bool=true,
    write_netcdf::Bool=true,
    netcdf_variables=copy(CASE_NETCDF_VARIABLES),
    cycles::Integer=10,
    backend=:threads,
    history_stride::Integer=1,
)
    resolved_name = String(name)
    resolved_output_dir = isempty(output_dir) ? _default_case_output_dir(resolved_name) : String(output_dir)
    resolved_input_label = isnothing(input_label) ? definition.input_label : String(input_label)
    run = RunConfig(
        name=resolved_name,
        input_label=resolved_input_label,
        output_dir=resolved_output_dir,
        netcdf_path=String(netcdf_path),
        write_outputs=write_outputs,
        write_netcdf=write_netcdf,
        netcdf_variables=netcdf_variables,
        cycles=Int(cycles),
        backend=backend,
        history_stride=Int(history_stride),
    )
    run.write_netcdf && isnothing(definition.layout) &&
        error("`write_netcdf=true` requires a grid layout, but this case definition has `layout=nothing`.")
    return SnowpackCase(run.name, definition, run)
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
    return SnowpackCase(definition; kwargs...)
end

"""
    synthetic_case(; ...)

Create a runnable synthetic [`SnowpackCase`](@ref) in one step.
"""
function synthetic_case(;
    physics::SM.SnowpackPhysicalConstants{Float64}=physics(),
    ntot::Integer=5,
    variant::Symbol=:multi_column,
    ntime::Integer=12,
    nx::Union{Nothing, Integer}=nothing,
    ny::Union{Nothing, Integer}=nothing,
    name::AbstractString="synthetic_case",
    input_label::Union{Nothing, AbstractString}=nothing,
    output_dir::AbstractString="",
    netcdf_path::AbstractString="",
    write_outputs::Bool=true,
    write_netcdf::Bool=true,
    netcdf_variables=copy(CASE_NETCDF_VARIABLES),
    cycles::Integer=10,
    backend=:threads,
    history_stride::Integer=1,
)
    definition = synthetic_definition(
        physics=physics,
        ntot=ntot,
        variant=variant,
        ntime=ntime,
        nx=nx,
        ny=ny,
    )
    return SnowpackCase(
        definition;
        name=name,
        input_label=input_label,
        output_dir=output_dir,
        netcdf_path=netcdf_path,
        write_outputs=write_outputs,
        write_netcdf=write_netcdf,
        netcdf_variables=netcdf_variables,
        cycles=cycles,
        backend=backend,
        history_stride=history_stride,
    )
end

"""
    prescribed_case(; ...)

Create a runnable [`SnowpackCase`](@ref) directly from user-facing forcing
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
    name::AbstractString="prescribed_case",
    output_dir::AbstractString="",
    netcdf_path::AbstractString="",
    write_outputs::Bool=true,
    write_netcdf::Bool=true,
    netcdf_variables=copy(CASE_NETCDF_VARIABLES),
    cycles::Integer=10,
    backend=:threads,
    history_stride::Integer=1,
)
    definition = prescribed_definition(
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
    return SnowpackCase(
        definition;
        name=name,
        input_label=input_label,
        output_dir=output_dir,
        netcdf_path=netcdf_path,
        write_outputs=write_outputs,
        write_netcdf=write_netcdf,
        netcdf_variables=netcdf_variables,
        cycles=cycles,
        backend=backend,
        history_stride=history_stride,
    )
end

"""
    mar_case(path; ...)

Create a runnable [`SnowpackCase`](@ref) directly from a MAR NetCDF/HDF5 file.
"""
function mar_case(
    path::AbstractString;
    physics::SM.SnowpackPhysicalConstants{Float64}=physics(),
    ntot::Integer=20,
    mask_threshold::Real=50.0,
    turbulent_flux_sign::Real=1.0,
    name::AbstractString="mar_case",
    input_label::Union{Nothing, AbstractString}=nothing,
    output_dir::AbstractString="",
    netcdf_path::AbstractString="",
    write_outputs::Bool=true,
    write_netcdf::Bool=true,
    netcdf_variables=copy(CASE_NETCDF_VARIABLES),
    cycles::Integer=10,
    backend=:threads,
    history_stride::Integer=1,
)
    definition = mar_definition(
        path;
        physics=physics,
        ntot=ntot,
        mask_threshold=mask_threshold,
        turbulent_flux_sign=turbulent_flux_sign,
    )
    return SnowpackCase(
        definition;
        name=name,
        input_label=input_label,
        output_dir=output_dir,
        netcdf_path=netcdf_path,
        write_outputs=write_outputs,
        write_netcdf=write_netcdf,
        netcdf_variables=netcdf_variables,
        cycles=cycles,
        backend=backend,
        history_stride=history_stride,
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
