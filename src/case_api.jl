using Dates
using HDF5
using Base.Threads: @threads

"""
    AbstractCaseSource

Extension point for case/forcing loaders used by [`load_case`](@ref) and
[`build_case`](@ref). Add a new forcing format by defining a subtype and a
matching `load_case(::YourSource; physics, ntot)` method that returns a
[`CaseDefinition`](@ref).
"""
abstract type AbstractCaseSource end

"""
    CaseDefinition

Reusable, execution-independent case inputs: an initial domain, case
forcing, an optional grid layout, and small metadata/notes for user-facing
workflows. Build this once with [`load_case`](@ref), then create CPU and GPU
cases from the same data with [`build_case`](@ref).
"""
struct CaseDefinition
    domain::SM.SnowpackDomain
    forcing::ForcingData
    layout::Union{Nothing, GridLayout}
    input_label::String
    notes::Vector{String}
    metadata::NamedTuple
end

function CaseDefinition(
    domain::SM.SnowpackDomain,
    forcing::ForcingData;
    layout::Union{Nothing, GridLayout}=nothing,
    input_label::AbstractString="",
    notes::AbstractVector{<:AbstractString}=String[],
    metadata::NamedTuple=(;),
)
    ncol = SM.column_count(domain)
    size(forcing.air_temperature, 1) == ncol ||
        error("Forcing column count ($(size(forcing.air_temperature, 1))) must match the domain column count ($ncol).")
    !isnothing(layout) && length(layout.js) != ncol &&
        error("Grid-layout point count ($(length(layout.js))) must match the domain column count ($ncol).")
    return CaseDefinition(
        domain,
        forcing,
        layout,
        String(input_label),
        String[String(note) for note in notes],
        metadata,
    )
end

"""
    SnowpackCase

High-level runnable case configuration. It keeps parsed inputs separate from run
configuration so the same case definition can be reused across CPU and GPU runs.
"""
struct SnowpackCase
    name::String
    definition::CaseDefinition
    run::RunConfig
end

function Base.show(io::IO, data::CaseDefinition)
    ncol = SM.column_count(data.domain)
    ntime = length(data.forcing.time_values)
    layout_state = isnothing(data.layout) ? "none" : "grid"
    print(io, "CaseDefinition(ncol=$(ncol), ntime=$(ntime), layout=$(layout_state))")
end

function Base.show(io::IO, case::SnowpackCase)
    ncol = SM.column_count(case.definition.domain)
    ntime = length(case.definition.forcing.time_values)
    layout_state = isnothing(case.definition.layout) ? "none" : "grid"
    print(
        io,
        "SnowpackCase(name=$(repr(case.name)), backend=$(case.run.backend), ncol=$(ncol), ntime=$(ntime), layout=$(layout_state))",
    )
end

"""
    SyntheticCaseSource

Descriptor for built-in synthetic forcing used by examples, smoke tests, and
documentation workflows.
"""
struct SyntheticCaseSource <: AbstractCaseSource
    variant::Symbol
    ntime::Int
    nx::Int
    ny::Int
end

"""
    SyntheticCaseSource(; variant=:multi_column, ntime=12, nx, ny)

Create a small synthetic forcing source for examples, smoke tests, and notebook
workflows. `variant` currently supports `:single_column` and `:multi_column`.
For multi-column cases, `nx` and `ny` control the grid size.
"""
function SyntheticCaseSource(;
    variant::Symbol=:multi_column,
    ntime::Integer=12,
    nx::Union{Nothing, Integer}=nothing,
    ny::Union{Nothing, Integer}=nothing,
)
    ntime > 0 || error("`ntime` must be positive.")
    variant in (:single_column, :multi_column) ||
        error("Unsupported synthetic forcing variant '$variant'. Use `:single_column` or `:multi_column`.")
    resolved_nx = isnothing(nx) ? (variant == :single_column ? 1 : 2) : Int(nx)
    resolved_ny = isnothing(ny) ? (variant == :single_column ? 1 : 2) : Int(ny)
    resolved_nx > 0 || error("`nx` must be positive.")
    resolved_ny > 0 || error("`ny` must be positive.")
    if variant == :single_column && (resolved_nx != 1 || resolved_ny != 1)
        error("`variant=:single_column` requires `nx=1` and `ny=1`.")
    end
    return SyntheticCaseSource(variant, Int(ntime), resolved_nx, resolved_ny)
end

struct MARCaseSource <: AbstractCaseSource
    path::String
    mask_threshold::Float64
    turbulent_flux_sign::Float64
end

"""
    MARCaseSource(path; mask_threshold=50.0, turbulent_flux_sign=1.0)

Create a forcing-source descriptor for MAR NetCDF/HDF5 forcing. Use it with
[`load_case`](@ref) or pass it directly to [`build_case`](@ref).
"""
function MARCaseSource(
    path::AbstractString;
    mask_threshold::Real=50.0,
    turbulent_flux_sign::Real=1.0,
)
    isempty(strip(path)) && error("`path` must point to a MAR NetCDF/HDF5 file.")
    return MARCaseSource(
        abspath(String(path)),
        Float64(mask_threshold),
        Float64(turbulent_flux_sign),
    )
end

@inline _case_symbol(value::Symbol) = value
@inline _case_symbol(value) = Symbol(lowercase(strip(String(value))))

"""
    physics(; albedo=:dynamic, densification=:bessi, fresh_snow_density=:constant, kwargs...)

Convenience constructor for [`SnowpackPhysicalConstants`](@ref) that accepts
the user-facing scheme keywords used by the case API.
"""
function physics(;
    albedo=:dynamic,
    densification=:bessi,
    fresh_snow_density=:constant,
    kwargs...,
)
    return SM.SnowpackPhysicalConstants(
        Float64;
        albedo_scheme=_case_symbol(albedo),
        low_density_densification=_case_symbol(densification),
        fresh_snow_density_scheme=_case_symbol(fresh_snow_density),
        kwargs...,
    )
end

function _case_slug(name::AbstractString)
    slug = replace(lowercase(strip(String(name))), r"[^a-z0-9]+" => "_")
    slug = strip(slug, '_')
    return isempty(slug) ? "snowpack_case" : slug
end

_default_case_output_dir(name::AbstractString) = joinpath(pwd(), "case_output", _case_slug(name))

function _regular_grid_layout(nx::Integer, ny::Integer)
    Int(nx) > 0 || error("`nx` must be positive.")
    Int(ny) > 0 || error("`ny` must be positive.")
    x = Float64[i - 1 for i in 1:Int(nx)]
    y = Float64[j - 1 for j in 1:Int(ny)]
    ncol = Int(nx) * Int(ny)
    js = Vector{Int}(undef, ncol)
    is = Vector{Int}(undef, ncol)
    idx = 1
    @inbounds for j in 1:Int(ny), i in 1:Int(nx)
        js[idx] = j
        is[idx] = i
        idx += 1
    end
    return GridLayout(x, y, js, is, ones(Float64, Int(ny), Int(nx)))
end

function _expand_column_vector(values, ncol::Int, name::AbstractString)
    if values isa Number
        return fill(Float64(values), ncol)
    end
    data = collect(values)
    if ndims(data) == 1
        length(data) == ncol || error("`$name` must be a scalar or a vector with length $ncol.")
        return Float64.(data)
    end
    error("`$name` must be a scalar or a vector with length $ncol.")
end

function _expand_numeric_timeseries(field, ncol::Int, ntime::Int, name::AbstractString)
    if field isa Number
        return fill(Float64(field), ncol, ntime)
    end
    data = collect(field)
    if ndims(data) == 1
        length(data) == ntime || error("`$name` must have length $ntime.")
        return repeat(reshape(Float64.(data), 1, ntime), ncol, 1)
    elseif ndims(data) == 2
        size(data) == (ncol, ntime) || error("`$name` must have size ($ncol, $ntime).")
        return Matrix{Float64}(data)
    end
    error("`$name` must be a scalar, a vector of length $ntime, or a matrix of size ($ncol, $ntime).")
end

function _expand_bool_timeseries(field, ncol::Int, ntime::Int, name::AbstractString)
    if field isa Bool
        return fill(field, ncol, ntime)
    end
    data = collect(field)
    if ndims(data) == 1
        length(data) == ntime || error("`$name` must have length $ntime.")
        return repeat(reshape(Bool.(data), 1, ntime), ncol, 1)
    elseif ndims(data) == 2
        size(data) == (ncol, ntime) || error("`$name` must have size ($ncol, $ntime).")
        return Bool.(data)
    end
    error("`$name` must be a Bool, a vector of length $ntime, or a matrix of size ($ncol, $ntime).")
end

function _synthetic_layout(variant::Symbol, nx::Int, ny::Int)
    if variant == :single_column
        return (
            x=[0.0],
            y=[0.0],
            js=[1],
            is=[1],
            mask=ones(Float64, 1, 1),
            surface_mass=[250.0],
            surface_temperature_offsets=[-12.0],
            surface_albedo_offsets=[0.0],
            label="Synthetic single-column seasonal forcing",
        )
    elseif variant == :multi_column
        ncol = nx * ny
        x = Float64[i - 1 for i in 1:nx]
        y = Float64[j - 1 for j in 1:ny]
        js = Vector{Int}(undef, ncol)
        is = Vector{Int}(undef, ncol)
        xfrac = Vector{Float64}(undef, ncol)
        yfrac = Vector{Float64}(undef, ncol)
        spatial_wave = Vector{Float64}(undef, ncol)
        surface_mass = Vector{Float64}(undef, ncol)
        surface_temperature_offsets = Vector{Float64}(undef, ncol)
        surface_albedo_offsets = Vector{Float64}(undef, ncol)
        idx = 1
        @inbounds for j in 1:ny, i in 1:nx
            xfrac_ij = nx == 1 ? 0.0 : (i - 1) / (nx - 1)
            yfrac_ij = ny == 1 ? 0.0 : (j - 1) / (ny - 1)
            wave = sinpi(xfrac_ij) * cospi(yfrac_ij)
            js[idx] = j
            is[idx] = i
            xfrac[idx] = xfrac_ij
            yfrac[idx] = yfrac_ij
            spatial_wave[idx] = wave
            surface_mass[idx] = 235.0 + 20.0 * xfrac_ij + 15.0 * yfrac_ij + 12.0 * wave
            surface_temperature_offsets[idx] = -13.0 + 1.8 * xfrac_ij - 1.2 * yfrac_ij + 0.6 * wave
            surface_albedo_offsets[idx] = clamp(0.015 * (0.5 - xfrac_ij) + 0.01 * wave, -0.03, 0.03)
            idx += 1
        end
        return (
            x=x,
            y=y,
            js=js,
            is=is,
            xfrac=xfrac,
            yfrac=yfrac,
            spatial_wave=spatial_wave,
            mask=ones(Float64, ny, nx),
            surface_mass=surface_mass,
            surface_temperature_offsets=surface_temperature_offsets,
            surface_albedo_offsets=surface_albedo_offsets,
            label="Synthetic $(nx)x$(ny) seasonal forcing",
        )
    end
    error("Unsupported synthetic forcing variant '$variant'. Use `:single_column` or `:multi_column`.")
end

"""
    load_case(source::AbstractCaseSource; physics=SnowpackPhysicalConstants(), ntot=20)

Parse or synthesize forcing inputs and return reusable [`CaseDefinition`](@ref).
This step is intentionally separate from [`run_case`](@ref) so users can inspect
or reuse the parsed inputs before execution.
"""
function load_case(
    source::SyntheticCaseSource;
    physics::SM.SnowpackPhysicalConstants{Float64}=SM.SnowpackPhysicalConstants(),
    ntot::Integer=5,
)
    ntot > 0 || error("`ntot` must be positive.")

    layout_data = _synthetic_layout(source.variant, source.nx, source.ny)
    ncol = length(layout_data.js)
    ntime = source.ntime
    N = ones(Int, ncol)
    mass = zeros(Float64, Int(ntot), ncol)
    mass[1, :] .= layout_data.surface_mass
    mass_w = zeros(Float64, Int(ntot), ncol)
    density = fill(320.0, Int(ntot), ncol)
    temperature = fill(physics.T0 - 12.0, Int(ntot), ncol)
    Tsrf = Vector{Float64}(undef, ncol)
    albedo_dynamic = Vector{Float64}(undef, ncol)

    @inbounds for col in 1:ncol
        temperature[1, col] = physics.T0 + layout_data.surface_temperature_offsets[col]
        Tsrf[col] = physics.T0 + layout_data.surface_temperature_offsets[col] + 2.0
        albedo_dynamic[col] = clamp(
            physics.alpha_dry + layout_data.surface_albedo_offsets[col],
            physics.alpha_ice,
            physics.alpha_dry,
        )
    end

    state = SnowpackStateFields(
        N,
        mass,
        mass_w,
        density,
        temperature;
        Tsrf=Tsrf,
        albedo_dynamic=albedo_dynamic,
        physics=physics,
    )

    time_values = [DateTime(2001, 1, 15, 12) + Day(30 * (t - 1)) for t in 1:ntime]
    dt_days = fill(30.0, ntime)
    air_temperature = zeros(Float64, ncol, ntime)
    snowfall_rate = zeros(Float64, ncol, ntime)
    rainfall_rate = zeros(Float64, ncol, ntime)
    shortwave_down = zeros(Float64, ncol, ntime)
    wind_speed = fill(5.0, ncol, ntime)
    q_lw_down = zeros(Float64, ncol, ntime)
    has_q_lw_down = fill(false, ncol, ntime)
    q_sh = zeros(Float64, ncol, ntime)
    has_q_sh = fill(false, ncol, ntime)
    q_lh = zeros(Float64, ncol, ntime)
    has_q_lh = fill(false, ncol, ntime)
    xfrac_values = hasproperty(layout_data, :xfrac) ? layout_data.xfrac : fill(0.0, ncol)
    yfrac_values = hasproperty(layout_data, :yfrac) ? layout_data.yfrac : fill(0.0, ncol)
    wave_values = hasproperty(layout_data, :spatial_wave) ? layout_data.spatial_wave : fill(0.0, ncol)

    @inbounds for col in 1:ncol, t in 1:ntime
        phase = 2π * (t - 1) / ntime
        xfrac = xfrac_values[col]
        yfrac = yfrac_values[col]
        wave = wave_values[col]
        temperature_offset = 1.8 * (xfrac - 0.5) - 2.4 * (yfrac - 0.5) + 0.8 * wave
        snowfall_scale = clamp(1.0 + 0.18 * (0.5 - yfrac) + 0.10 * wave, 0.75, 1.25)
        rainfall_scale = clamp(1.0 + 0.12 * xfrac - 0.08 * yfrac + 0.05 * wave, 0.80, 1.20)
        shortwave_scale = clamp(1.0 + 0.10 * (xfrac - 0.5) + 0.08 * wave, 0.85, 1.15)
        air_temperature[col, t] = physics.T0 - 14.0 + 11.0 * sin(phase) + temperature_offset
        snowfall_rate[col, t] = t in (1, 2, 3, max(ntime - 1, 1), ntime) ? 1.2e-5 * snowfall_scale : 1.0e-6 * snowfall_scale
        rainfall_rate[col, t] = t in max(1, ntime ÷ 2 - 1):min(ntime, ntime ÷ 2 + 1) ? 4.0e-6 * rainfall_scale : 0.0
        shortwave_down[col, t] = max(0.0, (180.0 + 130.0 * sin(phase - π / 2)) * shortwave_scale)
        wind_speed[col, t] = max(2.0, 4.5 + 1.0 * xfrac - 0.7 * yfrac + 0.6 * wave)
    end

    forcing = ForcingData(
        time_values=time_values,
        dt_days=dt_days,
        air_temperature=air_temperature,
        snowfall_rate=snowfall_rate,
        rainfall_rate=rainfall_rate,
        shortwave_down=shortwave_down,
        wind_speed=wind_speed,
        q_lw_down=q_lw_down,
        has_q_lw_down=has_q_lw_down,
        q_sh=q_sh,
        has_q_sh=has_q_sh,
        q_lh=q_lh,
        has_q_lh=has_q_lh,
    )
    layout = GridLayout(
        layout_data.x,
        layout_data.y,
        layout_data.js,
        layout_data.is,
        layout_data.mask,
    )
    metadata = (
        format=:synthetic,
        variant=source.variant,
        label=layout_data.label,
        nx=source.nx,
        ny=source.ny,
        grid_shape=(source.ny, source.nx),
        ncol=ncol,
        ntime=ntime,
        ntot=Int(ntot),
    )
    return CaseDefinition(
        SnowpackDomain(state),
        forcing;
        layout=layout,
        input_label=layout_data.label,
        metadata=metadata,
    )
end

function load_case(
    path::AbstractString;
    format=:mar,
    physics::SM.SnowpackPhysicalConstants{Float64}=SM.SnowpackPhysicalConstants(),
    ntot::Integer=20,
    kwargs...,
)
    format_symbol = _case_symbol(format)
    if format_symbol == :mar
        return load_case(MARCaseSource(path; kwargs...); physics=physics, ntot=ntot)
    end
    error(
        "Unsupported forcing format '$format'. " *
        "Available formats: `:mar`. " *
        "Add new formats by defining `load_case(::YourSource; physics, ntot)`.",
    )
end

load_case(data::CaseDefinition; kwargs...) = data

"""
    prescribed_case(; ...)

Create a case definition directly from user-supplied forcing arrays and simple
initial-condition keywords. High-level forcing inputs use Celsius and
mmWE/day.
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
    Int(ntot) > 0 || error("`ntot` must be positive.")
    dt_days_v = Float64.(collect(dt_days))
    isempty(dt_days_v) && error("`dt_days` must not be empty.")
    all(>(0.0), dt_days_v) || error("All `dt_days` entries must be positive.")
    ntime = length(dt_days_v)
    layout = _regular_grid_layout(nx, ny)
    ncol = length(layout.js)

    forcing = ForcingData(
        dt_days=dt_days_v,
        air_temperature_c=_expand_numeric_timeseries(air_temperature_c, ncol, ntime, "air_temperature_c"),
        snowfall_mm_day=_expand_numeric_timeseries(snowfall_mm_day, ncol, ntime, "snowfall_mm_day"),
        rainfall_mm_day=_expand_numeric_timeseries(rainfall_mm_day, ncol, ntime, "rainfall_mm_day"),
        shortwave_down=_expand_numeric_timeseries(shortwave_down, ncol, ntime, "shortwave_down"),
        wind_speed=_expand_numeric_timeseries(wind_speed, ncol, ntime, "wind_speed"),
        q_lw_down=isnothing(q_lw_down) ? nothing : _expand_numeric_timeseries(q_lw_down, ncol, ntime, "q_lw_down"),
        has_q_lw_down=isnothing(has_q_lw_down) ? nothing : _expand_bool_timeseries(has_q_lw_down, ncol, ntime, "has_q_lw_down"),
        q_sh=isnothing(q_sh) ? nothing : _expand_numeric_timeseries(q_sh, ncol, ntime, "q_sh"),
        has_q_sh=isnothing(has_q_sh) ? nothing : _expand_bool_timeseries(has_q_sh, ncol, ntime, "has_q_sh"),
        q_lh=isnothing(q_lh) ? nothing : _expand_numeric_timeseries(q_lh, ncol, ntime, "q_lh"),
        has_q_lh=isnothing(has_q_lh) ? nothing : _expand_bool_timeseries(has_q_lh, ncol, ntime, "has_q_lh"),
        time_values=time_values,
    )

    surface_mass = _expand_column_vector(initial_surface_mass, ncol, "initial_surface_mass")
    surface_density = _expand_column_vector(initial_density, ncol, "initial_density")
    surface_temperature_c = _expand_column_vector(initial_temperature_c, ncol, "initial_temperature_c")
    surface_albedo = _expand_column_vector(initial_albedo, ncol, "initial_albedo")

    N = Int[m > 0.0 ? 1 : 0 for m in surface_mass]
    mass = zeros(Float64, Int(ntot), ncol)
    mass_w = zeros(Float64, Int(ntot), ncol)
    density = zeros(Float64, Int(ntot), ncol)
    temperature = fill(physics.T0, Int(ntot), ncol)
    Tsrf = physics.T0 .+ surface_temperature_c
    snow_cover = Float64.(N .> 0)

    @inbounds for col in 1:ncol
        if N[col] > 0
            mass[1, col] = surface_mass[col]
            density[1, col] = surface_density[col]
            temperature[1, col] = physics.T0 + surface_temperature_c[col]
        end
    end

    state = SnowpackStateFields(
        N,
        mass,
        mass_w,
        density,
        temperature;
        Tsrf=Tsrf,
        snow_cover=snow_cover,
        albedo_dynamic=surface_albedo,
        physics=physics,
    )
    metadata = (
        format=:prescribed,
        nx=Int(nx),
        ny=Int(ny),
        grid_shape=(Int(ny), Int(nx)),
        ncol=ncol,
        ntime=ntime,
        ntot=Int(ntot),
    )
    return CaseDefinition(
        SnowpackDomain(state),
        forcing;
        layout=layout,
        input_label=String(input_label),
        metadata=metadata,
    )
end

"""
    synthetic_case(; physics=physics(), ntot=5, kwargs...)

Build a reusable synthetic [`CaseDefinition`](@ref) from
[`SyntheticCaseSource`](@ref) without constructing the source object manually.
"""
synthetic_case(;
    physics::SM.SnowpackPhysicalConstants{Float64}=physics(),
    ntot::Integer=5,
    kwargs...,
) = load_case(SyntheticCaseSource(; kwargs...); physics=physics, ntot=ntot)

"""
    mar_case(path; physics=physics(), ntot=20, kwargs...)

Load a reusable [`CaseDefinition`](@ref) from a MAR forcing file using
[`MARCaseSource`](@ref).
"""
mar_case(
    path::AbstractString;
    physics::SM.SnowpackPhysicalConstants{Float64}=physics(),
    ntot::Integer=20,
    kwargs...,
) = load_case(MARCaseSource(path; kwargs...); physics=physics, ntot=ntot)

function _resolve_build_case_value(provided, fallback)
    return isnothing(provided) ? fallback : provided
end

function _resolve_run_config(
    definition::CaseDefinition;
    name::Union{Nothing, AbstractString}=nothing,
    run::Union{Nothing, RunConfig}=nothing,
    backend=nothing,
    input_label::Union{Nothing, AbstractString}=nothing,
    output_dir::Union{Nothing, AbstractString}=nothing,
    netcdf_path::Union{Nothing, AbstractString}=nothing,
    write_outputs::Union{Nothing, Bool}=nothing,
    write_netcdf::Union{Nothing, Bool}=nothing,
    netcdf_variables=nothing,
    cycles::Union{Nothing, Integer}=nothing,
    history_stride::Union{Nothing, Integer}=nothing,
)
    base = isnothing(run) ? RunConfig() : run
    resolved_name = String(_resolve_build_case_value(name, base.name))
    resolved_output_dir = if !isnothing(output_dir)
        String(output_dir)
    elseif isempty(base.output_dir)
        _default_case_output_dir(resolved_name)
    else
        base.output_dir
    end
    return RunConfig(
        name=resolved_name,
        input_label=String(_resolve_build_case_value(input_label, isempty(base.input_label) ? definition.input_label : base.input_label)),
        output_dir=resolved_output_dir,
        netcdf_path=String(_resolve_build_case_value(netcdf_path, base.netcdf_path)),
        write_outputs=Bool(_resolve_build_case_value(write_outputs, base.write_outputs)),
        write_netcdf=Bool(_resolve_build_case_value(write_netcdf, base.write_netcdf)),
        netcdf_variables=_resolve_build_case_value(netcdf_variables, base.netcdf_variables),
        cycles=Int(_resolve_build_case_value(cycles, base.cycles)),
        history_stride=Int(_resolve_build_case_value(history_stride, base.history_stride)),
        backend=_resolve_build_case_value(backend, base.backend),
    )
end

"""
    build_case(definition::CaseDefinition; run=RunConfig(), ...)

Create a user-facing runnable case from a reusable [`CaseDefinition`](@ref).
Pass a [`RunConfig`](@ref) with the desired backend/output settings, and use
keyword overrides only when you need to tweak a small part of that config.
"""
function build_case(
    definition::CaseDefinition;
    name::Union{Nothing, AbstractString}=nothing,
    run::Union{Nothing, RunConfig}=nothing,
    backend=nothing,
    input_label::Union{Nothing, AbstractString}=nothing,
    output_dir::Union{Nothing, AbstractString}=nothing,
    netcdf_path::Union{Nothing, AbstractString}=nothing,
    write_outputs::Union{Nothing, Bool}=nothing,
    write_netcdf::Union{Nothing, Bool}=nothing,
    netcdf_variables=nothing,
    cycles::Union{Nothing, Integer}=nothing,
    history_stride::Union{Nothing, Integer}=nothing,
)
    run_config = _resolve_run_config(
        definition;
        name=name,
        run=run,
        backend=backend,
        input_label=input_label,
        output_dir=output_dir,
        netcdf_path=netcdf_path,
        write_outputs=write_outputs,
        write_netcdf=write_netcdf,
        netcdf_variables=netcdf_variables,
        cycles=cycles,
        history_stride=history_stride,
    )
    run_config.write_netcdf && isnothing(definition.layout) &&
        error("`write_netcdf=true` requires a grid layout, but this case definition has `layout=nothing`.")
    return SnowpackCase(run_config.name, definition, run_config)
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

const _MAR_FILL_THRESHOLD = -9.0e18

function _mar_read_dataset_shapes(path::AbstractString)
    shapes = Dict{String, Vector{Int}}()
    h5open(path, "r") do file
        for name in keys(file)
            obj = file[name]
            if obj isa HDF5.Dataset
                shapes[String(name)] = reverse(collect(size(obj)))
            end
        end
    end
    return shapes
end

function _mar_clean_fill!(A)
    @inbounds for i in eachindex(A)
        if A[i] <= _MAR_FILL_THRESHOLD
            A[i] = NaN
        end
    end
    return A
end

function _mar_read_hdf5_subset(
    path::AbstractString,
    varname::AbstractString,
    full_shape::Vector{Int};
    start::Union{Nothing, Vector{Int}}=nothing,
    count::Union{Nothing, Vector{Int}}=nothing,
)
    start_vec = isnothing(start) ? zeros(Int, length(full_shape)) : start
    count_vec = isnothing(count) ? copy(full_shape) : count
    length(start_vec) == length(full_shape) || error("start rank mismatch for variable '$varname'.")
    length(count_vec) == length(full_shape) || error("count rank mismatch for variable '$varname'.")
    data = h5open(path, "r") do file
        haskey(file, varname) || error("Variable '$varname' was not found in $(abspath(path)).")
        dataset = file[varname]
        nd = length(full_shape)
        file_ranges = ntuple(nd) do dim
            logical_dim = nd - dim + 1
            first_index = start_vec[logical_dim] + 1
            last_index = start_vec[logical_dim] + count_vec[logical_dim]
            first_index:last_index
        end
        raw = dataset[file_ranges...]
        logical = Float64.(raw)
        return nd > 1 ? permutedims(logical, nd:-1:1) : logical
    end
    _mar_clean_fill!(data)
    return data
end

function _mar_read_hdf5_full(path::AbstractString, varname::AbstractString, shapes::Dict{String, Vector{Int}})
    haskey(shapes, varname) || error("Variable '$varname' was not found in $(abspath(path)).")
    return _mar_read_hdf5_subset(path, varname, shapes[varname])
end

function _mar_read_timeslice_2d(
    path::AbstractString,
    varname::AbstractString,
    time_index::Int,
    shapes::Dict{String, Vector{Int}},
)
    haskey(shapes, varname) || error("Variable '$varname' was not found in $(abspath(path)).")
    shape = shapes[varname]
    if length(shape) == 3
        start = [time_index - 1, 0, 0]
        count = [1, shape[2], shape[3]]
    elseif length(shape) == 4
        start = [time_index - 1, 0, 0, 0]
        count = [1, 1, shape[3], shape[4]]
    else
        error("Variable '$varname' does not have a supported rank for 2D slicing.")
    end
    data = _mar_read_hdf5_subset(path, varname, shape; start=start, count=count)
    return dropdims(data; dims=Tuple(findall(==(1), size(data))))
end

function _mar_read_timeslice_3d(
    path::AbstractString,
    varname::AbstractString,
    time_index::Int,
    shapes::Dict{String, Vector{Int}},
)
    haskey(shapes, varname) || error("Variable '$varname' was not found in $(abspath(path)).")
    shape = shapes[varname]
    length(shape) == 4 || error("Variable '$varname' does not have the expected 4D layout.")
    start = [time_index - 1, 0, 0, 0]
    count = [1, shape[2], shape[3], shape[4]]
    data = _mar_read_hdf5_subset(path, varname, shape; start=start, count=count)
    return dropdims(data; dims=(1,))
end

@inline _mar_valid_or(default::Float64, x::Float64) = isfinite(x) ? x : default
@inline _mar_mmwe_day_to_kgm2s(x::Float64) = isfinite(x) ? max(x, 0.0) / 86_400.0 : 0.0

function _mar_read_times(path::AbstractString, shapes::Dict{String, Vector{Int}})
    yyyy = round.(Int, vec(_mar_read_hdf5_full(path, "YYYY", shapes)))
    mm = round.(Int, vec(_mar_read_hdf5_full(path, "MM", shapes)))
    dd = round.(Int, vec(_mar_read_hdf5_full(path, "DD", shapes)))
    hh = round.(Int, vec(_mar_read_hdf5_full(path, "HH", shapes)))
    ntime = length(yyyy)
    times = Vector{DateTime}(undef, ntime)
    for i in 1:ntime
        times[i] = DateTime(yyyy[i], mm[i], dd[i], hh[i])
    end
    return times
end

function _mar_infer_dt_days(time_values::Vector{DateTime}, time_index::Int)
    if length(time_values) == 1
        return 1.0
    elseif time_index < length(time_values)
        return Dates.value(time_values[time_index + 1] - time_values[time_index]) / (1000 * 60 * 60 * 24)
    else
        return Dates.value(time_values[time_index] - time_values[time_index - 1]) / (1000 * 60 * 60 * 24)
    end
end

function _mar_extract_layers(
    total_height::Float64,
    density_profile::AbstractVector{<:Real},
    temperature_profile_c::AbstractVector{<:Real},
    liquid_water_profile::AbstractVector{<:Real},
    outlay_bounds::AbstractMatrix{<:Real};
    ntot::Int=80,
    c::SM.SnowpackPhysicalConstants=SM.SnowpackPhysicalConstants(),
    mass_split::Float64=SM.DEFAULT_MASS_SPLIT,
)
    if !isfinite(total_height) || total_height <= 0.0
        return (
            N=0,
            mass=Float64[],
            mass_w=Float64[],
            density=Float64[],
            temperature=Float64[],
        )
    end

    layer_mass = Float64[]
    layer_mass_w = Float64[]
    layer_density = Float64[]
    layer_temperature = Float64[]

    function append_restart_layer!(
        snow_mass::Float64,
        liquid_mass::Float64,
        density::Float64,
        temperature::Float64,
    )
        snow_mass <= 0.0 && return
        push!(layer_mass, snow_mass)
        push!(layer_mass_w, liquid_mass)
        push!(layer_density, density)
        push!(layer_temperature, temperature)
        return
    end

    function merge_restart_bottom_pair!()
        n = length(layer_mass)
        n >= 2 || return

        lower_mass = layer_mass[n - 1]
        bottom_mass = layer_mass[n]
        combined_mass = lower_mass + bottom_mass
        if combined_mass <= 0.0
            layer_mass[n - 1] = 0.0
            layer_mass_w[n - 1] += layer_mass_w[n]
            layer_density[n - 1] = max(layer_density[n - 1], layer_density[n])
            layer_temperature[n - 1] = min(layer_temperature[n - 1], layer_temperature[n])
        else
            layer_mass[n - 1] = combined_mass
            layer_mass_w[n - 1] += layer_mass_w[n]
            layer_density[n - 1] =
                (lower_mass * layer_density[n - 1] + bottom_mass * layer_density[n]) / combined_mass
            layer_temperature[n - 1] =
                (lower_mass * layer_temperature[n - 1] + bottom_mass * layer_temperature[n]) / combined_mass
        end

        pop!(layer_mass)
        pop!(layer_mass_w)
        pop!(layer_density)
        pop!(layer_temperature)
        return
    end

    nlayers = size(outlay_bounds, 1)
    for k in 1:nlayers
        lower = max(outlay_bounds[k, 1], 0.0)
        upper = outlay_bounds[k, 2]
        depth_cap = if k == nlayers && total_height > upper
            total_height
        else
            min(total_height, upper)
        end
        layer_thickness = max(depth_cap - lower, 0.0)
        layer_thickness <= 0.0 && continue

        rho = clamp(_mar_valid_or(300.0, Float64(density_profile[k])), 50.0, c.rho_i)
        temp_k = clamp(_mar_valid_or(-10.0, Float64(temperature_profile_c[k])) + c.T0, 200.0, c.T0)
        water_fraction = max(_mar_valid_or(0.0, Float64(liquid_water_profile[k])), 0.0)
        snow_mass = rho * layer_thickness
        liquid_mass = water_fraction * snow_mass
        append_restart_layer!(snow_mass, liquid_mass, rho, temp_k)
    end

    while length(layer_mass) > ntot
        merge_restart_bottom_pair!()
    end

    return (
        N=length(layer_mass),
        mass=layer_mass,
        mass_w=layer_mass_w,
        density=layer_density,
        temperature=layer_temperature,
    )
end

function _mar_populate_domain_column_from_restart!(
    domain::SM.SnowpackDomain,
    idx::Int,
    total_height::Float64,
    density_profile::AbstractVector{<:Real},
    temperature_profile_c::AbstractVector{<:Real},
    liquid_water_profile::AbstractVector{<:Real},
    outlay_bounds::AbstractMatrix{<:Real},
)
    extracted = _mar_extract_layers(
        total_height,
        density_profile,
        temperature_profile_c,
        liquid_water_profile,
        outlay_bounds;
        ntot=domain.Ntot,
        c=domain.c,
        mass_split=domain.mass_split,
    )
    domain.N[idx] = extracted.N
    if extracted.N > 0
        domain.mass[1:extracted.N, idx] .= extracted.mass
        domain.mass_w[1:extracted.N, idx] .= extracted.mass_w
        domain.density[1:extracted.N, idx] .= extracted.density
        domain.temperature[1:extracted.N, idx] .= extracted.temperature
        if extracted.N < domain.Ntot
            domain.mass[(extracted.N + 1):end, idx] .= 0.0
            domain.mass_w[(extracted.N + 1):end, idx] .= 0.0
            domain.density[(extracted.N + 1):end, idx] .= SM.DEFAULT_DENSITY_INIT
            domain.temperature[(extracted.N + 1):end, idx] .= domain.c.T0 - 10.0
        end
        domain.Tsrf[idx] = extracted.temperature[1]
    else
        domain.mass[:, idx] .= 0.0
        domain.mass_w[:, idx] .= 0.0
        domain.density[:, idx] .= SM.DEFAULT_DENSITY_INIT
        domain.temperature[:, idx] .= domain.c.T0 - 10.0
        domain.Tsrf[idx] = domain.c.T0 - 10.0
    end
    domain.mass_base[idx] = 0.0
    domain.smb_ice[idx] = 0.0
    domain.runoff[idx] = 0.0
    domain.snow_cover[idx] = 0.0
    domain.albedo_dynamic[idx] = domain.c.alpha_dry
    return domain
end

function _mar_read_full_timeseries_3d(
    path::AbstractString,
    varname::AbstractString,
    shapes::Dict{String, Vector{Int}},
)
    haskey(shapes, varname) || error("Variable '$varname' was not found in $(abspath(path)).")
    shape = shapes[varname]
    if length(shape) == 4
        data = _mar_read_hdf5_full(path, varname, shapes)
        size(data, 2) == 1 || error("Variable '$varname' has an unexpected non-singleton vertical dimension.")
        return dropdims(data; dims=(2,))
    elseif length(shape) == 3
        return _mar_read_hdf5_full(path, varname, shapes)
    end
    error("Variable '$varname' does not have a supported timeseries layout.")
end

function _mar_read_first_available_timeseries_3d(
    path::AbstractString,
    candidate_names::Vector{String},
    shapes::Dict{String, Vector{Int}},
)
    for name in candidate_names
        if haskey(shapes, name)
            return (name=name, data=_mar_read_full_timeseries_3d(path, name, shapes))
        end
    end
    return nothing
end

function load_case(
    source::MARCaseSource;
    physics::SM.SnowpackPhysicalConstants{Float64}=SM.SnowpackPhysicalConstants(),
    ntot::Integer=20,
)
    ntot > 0 || error("`ntot` must be positive.")
    isfile(source.path) || error("Forcing file was not found: $(source.path)")

    shapes = _mar_read_dataset_shapes(source.path)
    required_variables = ("x", "y", "MSK", "OUTLAY_bnds", "TT", "SF", "RF", "SWD", "LWD", "SHF", "LHF", "ZN3", "RO1", "TI1", "WA1", "YYYY", "MM", "DD", "HH")
    missing = String[var for var in required_variables if !haskey(shapes, var)]
    isempty(missing) || error(
        "MAR forcing file $(abspath(source.path)) is missing required variables: $(join(missing, ", "))."
    )

    time_values = _mar_read_times(source.path, shapes)
    dt_days = [_mar_infer_dt_days(time_values, t) for t in eachindex(time_values)]

    x = vec(_mar_read_hdf5_full(source.path, "x", shapes))
    y = vec(_mar_read_hdf5_full(source.path, "y", shapes))
    mask = _mar_read_hdf5_full(source.path, "MSK", shapes)
    outlay_bounds = _mar_read_hdf5_full(source.path, "OUTLAY_bnds", shapes)

    tt_full = _mar_read_full_timeseries_3d(source.path, "TT", shapes)
    sf_full = _mar_read_full_timeseries_3d(source.path, "SF", shapes)
    rf_full = _mar_read_full_timeseries_3d(source.path, "RF", shapes)
    swd_full = _mar_read_full_timeseries_3d(source.path, "SWD", shapes)
    lwd_full = _mar_read_full_timeseries_3d(source.path, "LWD", shapes)
    shf_full = _mar_read_full_timeseries_3d(source.path, "SHF", shapes) .* source.turbulent_flux_sign
    lhf_full = _mar_read_full_timeseries_3d(source.path, "LHF", shapes) .* source.turbulent_flux_sign

    u_wind_info = _mar_read_first_available_timeseries_3d(source.path, ["UU", "U10"], shapes)
    v_wind_info = _mar_read_first_available_timeseries_3d(source.path, ["VV", "V10"], shapes)
    wind_full = if !isnothing(u_wind_info) && !isnothing(v_wind_info)
        hypot.(u_wind_info.data, v_wind_info.data)
    else
        nothing
    end
    wind_note = if isnothing(wind_full)
        "Wind forcing: MAR wind components not found; using default 5.0 m s^-1."
    else
        @sprintf(
            "Wind forcing: |V| from MAR components %s and %s.",
            u_wind_info.name,
            v_wind_info.name,
        )
    end

    zn3_init = _mar_read_timeslice_2d(source.path, "ZN3", 1, shapes)
    ro1_init = _mar_read_timeslice_3d(source.path, "RO1", 1, shapes)
    ti1_init = _mar_read_timeslice_3d(source.path, "TI1", 1, shapes)
    wa1_init = _mar_read_timeslice_3d(source.path, "WA1", 1, shapes)

    ny, nx = size(mask)
    valid_mask = falses(ny, nx)
    @inbounds for j in 1:ny, i in 1:nx
        valid_mask[j, i] =
            isfinite(mask[j, i]) &&
            mask[j, i] >= source.mask_threshold &&
            isfinite(tt_full[1, j, i])
    end
    valid_indices = findall(valid_mask)
    nvalid = length(valid_indices)
    nvalid > 0 || error(
        "No valid MAR grid cells remain after applying `mask_threshold=$(source.mask_threshold)`."
    )
    ntime = length(time_values)

    js = Vector{Int}(undef, nvalid)
    is = Vector{Int}(undef, nvalid)
    domain = SM.SnowpackDomain(ncol=nvalid, Ntot=Int(ntot), c=physics)
    tair_k = Matrix{Float64}(undef, nvalid, ntime)
    snow_rate = Matrix{Float64}(undef, nvalid, ntime)
    rain_rate = Matrix{Float64}(undef, nvalid, ntime)
    s_boa = Matrix{Float64}(undef, nvalid, ntime)
    q_lw = Matrix{Float64}(undef, nvalid, ntime)
    has_q_lw = fill(false, nvalid, ntime)
    q_sh = Matrix{Float64}(undef, nvalid, ntime)
    has_q_sh = fill(false, nvalid, ntime)
    q_lh = Matrix{Float64}(undef, nvalid, ntime)
    has_q_lh = fill(false, nvalid, ntime)
    wind_speed = Matrix{Float64}(undef, nvalid, ntime)

    @threads :static for idx in eachindex(valid_indices)
        j, i = Tuple(valid_indices[idx])
        js[idx] = j
        is[idx] = i
        _mar_populate_domain_column_from_restart!(
            domain,
            idx,
            Float64(zn3_init[j, i]),
            @view(ro1_init[:, j, i]),
            @view(ti1_init[:, j, i]),
            @view(wa1_init[:, j, i]),
            outlay_bounds,
        )
        for t in 1:ntime
            tair_k[idx, t] = _mar_valid_or(-15.0, Float64(tt_full[t, j, i])) + domain.c.T0
            snow_rate[idx, t] = _mar_mmwe_day_to_kgm2s(Float64(sf_full[t, j, i]))
            rain_rate[idx, t] = _mar_mmwe_day_to_kgm2s(Float64(rf_full[t, j, i]))
            s_boa[idx, t] = _mar_valid_or(0.0, Float64(swd_full[t, j, i]))
            q_lw_ij = Float64(lwd_full[t, j, i])
            q_sh_ij = Float64(shf_full[t, j, i])
            q_lh_ij = Float64(lhf_full[t, j, i])
            has_q_lw[idx, t] = isfinite(q_lw_ij)
            has_q_sh[idx, t] = isfinite(q_sh_ij)
            has_q_lh[idx, t] = isfinite(q_lh_ij)
            q_lw[idx, t] = has_q_lw[idx, t] ? q_lw_ij : 0.0
            q_sh[idx, t] = has_q_sh[idx, t] ? q_sh_ij : 0.0
            q_lh[idx, t] = has_q_lh[idx, t] ? q_lh_ij : 0.0
            wind_speed[idx, t] = isnothing(wind_full) ? 5.0 : _mar_valid_or(5.0, Float64(wind_full[t, j, i]))
        end
    end

    forcing = ForcingData(
        time_values=time_values,
        dt_days=dt_days,
        air_temperature=tair_k,
        snowfall_rate=snow_rate,
        rainfall_rate=rain_rate,
        shortwave_down=s_boa,
        wind_speed=wind_speed,
        q_lw_down=q_lw,
        has_q_lw_down=has_q_lw,
        q_sh=q_sh,
        has_q_sh=has_q_sh,
        q_lh=q_lh,
        has_q_lh=has_q_lh,
    )
    layout = GridLayout(x, y, js, is, mask)
    metadata = (
        format=:mar,
        path=source.path,
        mask_threshold=source.mask_threshold,
        turbulent_flux_sign=source.turbulent_flux_sign,
        ncol=nvalid,
        ntime=ntime,
        ntot=Int(ntot),
    )
    return CaseDefinition(
        domain,
        forcing;
        layout=layout,
        input_label=source.path,
        notes=[wind_note],
        metadata=metadata,
    )
end
