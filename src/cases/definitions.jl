"""
    CaseDefinition

Reusable, execution-independent case inputs: an initial domain, case forcing,
an optional grid layout, and small metadata/notes for user-facing workflows.
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

High-level runnable case configuration. Parsed inputs live in
[`CaseDefinition`](@ref); execution settings live in [`RunConfig`](@ref).
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

@inline _case_symbol(value::Symbol) = value
@inline _case_symbol(value) = Symbol(lowercase(strip(String(value))))

"""
    physics(; albedo=:dynamic, densification=:bessi, fresh_snow_density=:constant, kwargs...)

Convenience constructor for [`SnowpackPhysicalConstants`](@ref) that accepts
the user-facing scheme keywords used by the case helpers.
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
            xfrac=[0.0],
            yfrac=[0.0],
            spatial_wave=[0.0],
            mask=ones(Float64, 1, 1),
            surface_mass=[250.0],
            surface_temperature_offsets=[-12.0],
            surface_albedo_offsets=[0.0],
            label="Synthetic single-column forcing",
        )
    elseif variant == :multi_column
        x = Float64[i - 1 for i in 1:nx]
        y = Float64[j - 1 for j in 1:ny]
        ncol = nx * ny
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
    synthetic_definition(; physics=physics(), ntot=5, variant=:multi_column, ntime=12, nx, ny)

Build a reusable synthetic [`CaseDefinition`](@ref).
"""
function synthetic_definition(;
    physics::SM.SnowpackPhysicalConstants{Float64}=physics(),
    ntot::Integer=5,
    variant::Symbol=:multi_column,
    ntime::Integer=12,
    nx::Union{Nothing, Integer}=nothing,
    ny::Union{Nothing, Integer}=nothing,
)
    Int(ntot) > 0 || error("`ntot` must be positive.")
    Int(ntime) > 0 || error("`ntime` must be positive.")
    variant in (:single_column, :multi_column) ||
        error("Unsupported synthetic forcing variant '$variant'. Use `:single_column` or `:multi_column`.")
    resolved_nx = isnothing(nx) ? (variant == :single_column ? 1 : 2) : Int(nx)
    resolved_ny = isnothing(ny) ? (variant == :single_column ? 1 : 2) : Int(ny)
    resolved_nx > 0 || error("`nx` must be positive.")
    resolved_ny > 0 || error("`ny` must be positive.")
    if variant == :single_column && (resolved_nx != 1 || resolved_ny != 1)
        error("`variant=:single_column` requires `nx=1` and `ny=1`.")
    end

    layout_data = _synthetic_layout(variant, resolved_nx, resolved_ny)
    ncol = length(layout_data.js)
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

    time_values = [DateTime(2001, 1, 15, 12) + Day(30 * (t - 1)) for t in 1:Int(ntime)]
    dt_days = fill(30.0, Int(ntime))
    air_temperature = zeros(Float64, ncol, Int(ntime))
    snowfall_rate = zeros(Float64, ncol, Int(ntime))
    rainfall_rate = zeros(Float64, ncol, Int(ntime))
    shortwave_down = zeros(Float64, ncol, Int(ntime))
    wind_speed = fill(5.0, ncol, Int(ntime))
    q_lw_down = zeros(Float64, ncol, Int(ntime))
    has_q_lw_down = fill(false, ncol, Int(ntime))
    q_sh = zeros(Float64, ncol, Int(ntime))
    has_q_sh = fill(false, ncol, Int(ntime))
    q_lh = zeros(Float64, ncol, Int(ntime))
    has_q_lh = fill(false, ncol, Int(ntime))
    xfrac_values = hasproperty(layout_data, :xfrac) ? layout_data.xfrac : fill(0.0, ncol)
    yfrac_values = hasproperty(layout_data, :yfrac) ? layout_data.yfrac : fill(0.0, ncol)
    wave_values = hasproperty(layout_data, :spatial_wave) ? layout_data.spatial_wave : fill(0.0, ncol)

    @inbounds for col in 1:ncol, t in 1:Int(ntime)
        phase = 2π * (t - 1) / Int(ntime)
        xfrac = xfrac_values[col]
        yfrac = yfrac_values[col]
        wave = wave_values[col]
        temperature_offset = 1.8 * (xfrac - 0.5) - 2.4 * (yfrac - 0.5) + 0.8 * wave
        snowfall_scale = clamp(1.0 + 0.18 * (0.5 - yfrac) + 0.10 * wave, 0.75, 1.25)
        rainfall_scale = clamp(1.0 + 0.12 * xfrac - 0.08 * yfrac + 0.05 * wave, 0.80, 1.20)
        shortwave_scale = clamp(1.0 + 0.10 * (xfrac - 0.5) + 0.08 * wave, 0.85, 1.15)
        air_temperature[col, t] = physics.T0 - 14.0 + 11.0 * sin(phase) + temperature_offset
        snowfall_rate[col, t] = t in (1, 2, 3, max(Int(ntime) - 1, 1), Int(ntime)) ? 1.2e-5 * snowfall_scale : 1.0e-6 * snowfall_scale
        rainfall_rate[col, t] = t in max(1, Int(ntime) ÷ 2 - 1):min(Int(ntime), Int(ntime) ÷ 2 + 1) ? 4.0e-6 * rainfall_scale : 0.0
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
        variant=variant,
        label=layout_data.label,
        nx=resolved_nx,
        ny=resolved_ny,
        grid_shape=(resolved_ny, resolved_nx),
        ncol=ncol,
        ntime=Int(ntime),
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

"""
    prescribed_definition(; ...)

Create a case definition directly from user-supplied forcing arrays and simple
initial-condition keywords. High-level forcing inputs use Celsius and mmWE/day.
"""
function prescribed_definition(;
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
