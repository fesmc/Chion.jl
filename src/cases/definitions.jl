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
    isnothing(layout) || length(layout.js) == ncol ||
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

struct SnowpackCase
    name::String
    definition::CaseDefinition
    run::RunConfig
end

function Base.show(io::IO, data::CaseDefinition)
    layout_state = isnothing(data.layout) ? "none" : "grid"
    print(io, "CaseDefinition(ncol=$(SM.column_count(data.domain)), ntime=$(length(data.forcing.time_values)), layout=$(layout_state))")
end

function Base.show(io::IO, case::SnowpackCase)
    layout_state = isnothing(case.definition.layout) ? "none" : "grid"
    print(
        io,
        "SnowpackCase(name=$(repr(case.name)), backend=$(case.run.backend), ncol=$(SM.column_count(case.definition.domain)), ntime=$(length(case.definition.forcing.time_values)), layout=$(layout_state))",
    )
end

@inline _case_symbol(value::Symbol) = value
@inline _case_symbol(value) = Symbol(lowercase(strip(String(value))))

function physics(; albedo=:dynamic, densification=:bessi, fresh_snow_density=:constant, kwargs...)
    return SM.SnowpackPhysicalConstants(
        Float64;
        albedo_scheme=_case_symbol(albedo),
        low_density_densification=_case_symbol(densification),
        fresh_snow_density_scheme=_case_symbol(fresh_snow_density),
        kwargs...,
    )
end

function _grid_axes(nx::Int, ny::Int)
    nx > 0 || error("`nx` must be positive.")
    ny > 0 || error("`ny` must be positive.")
    x = Float64.(0:nx-1)
    y = Float64.(0:ny-1)
    js = repeat(collect(1:ny), inner=nx)
    is = repeat(collect(1:nx), outer=ny)
    return x, y, js, is
end

function _regular_grid_layout(nx::Integer, ny::Integer)
    x, y, js, is = _grid_axes(Int(nx), Int(ny))
    return GridLayout(x, y, js, is, ones(Float64, length(y), length(x)))
end

function _column_vector(values, ncol::Int, name::AbstractString)
    values isa Number && return fill(Float64(values), ncol)
    data = collect(values)
    ndims(data) == 1 && length(data) == ncol || error("`$name` must be a scalar or a vector with length $ncol.")
    return Float64.(data)
end

function _surface_state(;
    ntot::Integer,
    surface_mass,
    surface_density,
    layer_temperature_c,
    surface_Tsrf_c=layer_temperature_c,
    surface_albedo,
    physics::SM.SnowpackPhysicalConstants{Float64},
)
    ncol = length(surface_mass)
    N = Int.(surface_mass .> 0.0)
    mass = zeros(Float64, Int(ntot), ncol)
    mass_w = zeros(Float64, Int(ntot), ncol)
    density = zeros(Float64, Int(ntot), ncol)
    temperature = fill(physics.T0, Int(ntot), ncol)
    @inbounds for col in eachindex(surface_mass)
        N[col] == 0 && continue
        mass[1, col] = surface_mass[col]
        density[1, col] = surface_density[col]
        temperature[1, col] = physics.T0 + layer_temperature_c[col]
    end
    return SnowpackStateFields(
        N,
        mass,
        mass_w,
        density,
        temperature;
        Tsrf=physics.T0 .+ surface_Tsrf_c,
        snow_cover=Float64.(N .> 0),
        albedo_dynamic=surface_albedo,
        physics=physics,
    )
end

function _synthetic_layout(variant::Symbol, nx::Int, ny::Int)
    variant in (:single_column, :multi_column) ||
        error("Unsupported synthetic forcing variant '$variant'. Use `:single_column` or `:multi_column`.")
    variant == :single_column && (nx == 1 && ny == 1) ||
        variant != :single_column ||
        error("`variant=:single_column` requires `nx=1` and `ny=1`.")

    x, y, js, is = _grid_axes(nx, ny)
    ncol = length(js)
    xfrac = nx == 1 ? zeros(Float64, ncol) : repeat(range(0.0, 1.0; length=nx), outer=ny)
    yfrac = ny == 1 ? zeros(Float64, ncol) : repeat(range(0.0, 1.0; length=ny), inner=nx)
    spatial_wave = sinpi.(xfrac) .* cospi.(yfrac)
    label = variant == :single_column ? "Synthetic single-column forcing" : "Synthetic $(nx)x$(ny) seasonal forcing"
    return (
        x=x,
        y=y,
        js=js,
        is=is,
        xfrac=xfrac,
        yfrac=yfrac,
        spatial_wave=spatial_wave,
        mask=ones(Float64, ny, nx),
        surface_mass=variant == :single_column ? fill(250.0, ncol) : 235.0 .+ 20.0 .* xfrac .+ 15.0 .* yfrac .+ 12.0 .* spatial_wave,
        surface_temperature_offsets=variant == :single_column ? fill(-12.0, ncol) : -13.0 .+ 1.8 .* xfrac .- 1.2 .* yfrac .+ 0.6 .* spatial_wave,
        surface_albedo_offsets=variant == :single_column ? zeros(Float64, ncol) : clamp.(0.015 .* (0.5 .- xfrac) .+ 0.01 .* spatial_wave, -0.03, 0.03),
        label=label,
    )
end

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
    resolved_nx = isnothing(nx) ? (variant == :single_column ? 1 : 2) : Int(nx)
    resolved_ny = isnothing(ny) ? (variant == :single_column ? 1 : 2) : Int(ny)
    layout_data = _synthetic_layout(variant, resolved_nx, resolved_ny)
    ncol = length(layout_data.js)

    surface_mass = layout_data.surface_mass
    surface_temperature_c = layout_data.surface_temperature_offsets
    state = _surface_state(
        ntot=ntot,
        surface_mass=surface_mass,
        surface_density=fill(320.0, ncol),
        layer_temperature_c=surface_temperature_c,
        surface_Tsrf_c=surface_temperature_c .+ 2.0,
        surface_albedo=clamp.(physics.alpha_dry .+ layout_data.surface_albedo_offsets, physics.alpha_ice, physics.alpha_dry),
        physics=physics,
    )

    phases = 2π .* (0:Int(ntime)-1) ./ Int(ntime)
    xfrac = reshape(layout_data.xfrac, ncol, 1)
    yfrac = reshape(layout_data.yfrac, ncol, 1)
    wave = reshape(layout_data.spatial_wave, ncol, 1)
    temperature_offset = 1.8 .* (xfrac .- 0.5) .- 2.4 .* (yfrac .- 0.5) .+ 0.8 .* wave
    snowfall_scale = clamp.(1.0 .+ 0.18 .* (0.5 .- yfrac) .+ 0.10 .* wave, 0.75, 1.25)
    rainfall_scale = clamp.(1.0 .+ 0.12 .* xfrac .- 0.08 .* yfrac .+ 0.05 .* wave, 0.80, 1.20)
    shortwave_scale = clamp.(1.0 .+ 0.10 .* (xfrac .- 0.5) .+ 0.08 .* wave, 0.85, 1.15)
    seasonal_temperature = reshape(11.0 .* sin.(phases), 1, Int(ntime))
    snow_profile = reshape([t in (1, 2, 3, max(Int(ntime) - 1, 1), Int(ntime)) ? 1.2e-5 : 1.0e-6 for t in 1:Int(ntime)], 1, Int(ntime))
    rain_profile = reshape([t in max(1, Int(ntime) ÷ 2 - 1):min(Int(ntime), Int(ntime) ÷ 2 + 1) ? 4.0e-6 : 0.0 for t in 1:Int(ntime)], 1, Int(ntime))
    shortwave_profile = reshape(max.(0.0, 180.0 .+ 130.0 .* sin.(phases .- π / 2)), 1, Int(ntime))
    wind_profile = max.(2.0, 4.5 .+ xfrac .- 0.7 .* yfrac .+ 0.6 .* wave)

    forcing = ForcingData(
        time_values=[DateTime(2001, 1, 15, 12) + Day(30 * (t - 1)) for t in 1:Int(ntime)],
        dt_days=fill(30.0, Int(ntime)),
        air_temperature=physics.T0 .- 14.0 .+ temperature_offset .+ seasonal_temperature,
        snowfall_rate=snowfall_scale .* snow_profile,
        rainfall_rate=rainfall_scale .* rain_profile,
        shortwave_down=shortwave_scale .* shortwave_profile,
        wind_speed=repeat(wind_profile, 1, Int(ntime)),
        q_lw_down=zeros(Float64, ncol, Int(ntime)),
        has_q_lw_down=fill(false, ncol, Int(ntime)),
        q_sh=zeros(Float64, ncol, Int(ntime)),
        has_q_sh=fill(false, ncol, Int(ntime)),
        q_lh=zeros(Float64, ncol, Int(ntime)),
        has_q_lh=fill(false, ncol, Int(ntime)),
    )

    return CaseDefinition(
        SnowpackDomain(state),
        forcing;
        layout=GridLayout(layout_data.x, layout_data.y, layout_data.js, layout_data.is, layout_data.mask),
        input_label=layout_data.label,
        metadata=(
            format=:synthetic,
            variant=variant,
            label=layout_data.label,
            nx=resolved_nx,
            ny=resolved_ny,
            grid_shape=(resolved_ny, resolved_nx),
            ncol=ncol,
            ntime=Int(ntime),
            ntot=Int(ntot),
        ),
    )
end

function prescribed_definition(;
    physics::SM.SnowpackPhysicalConstants{Float64}=physics(),
    ntot::Integer=15,
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
    layout = _regular_grid_layout(nx, ny)
    ncol = length(layout.js)
    forcing = ForcingData(
        dt_days=dt_days,
        ncol=ncol,
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
    )
    state = _surface_state(
        ntot=ntot,
        surface_mass=_column_vector(initial_surface_mass, ncol, "initial_surface_mass"),
        surface_density=_column_vector(initial_density, ncol, "initial_density"),
        layer_temperature_c=_column_vector(initial_temperature_c, ncol, "initial_temperature_c"),
        surface_albedo=_column_vector(initial_albedo, ncol, "initial_albedo"),
        physics=physics,
    )
    return CaseDefinition(
        SnowpackDomain(state),
        forcing;
        layout=layout,
        input_label=String(input_label),
        metadata=(
            format=:prescribed,
            source=:direct,
            nx=Int(nx),
            ny=Int(ny),
            grid_shape=(Int(ny), Int(nx)),
            ncol=ncol,
            ntime=length(forcing.time_values),
            ntot=Int(ntot),
        ),
    )
end
