"""
Core snowpack types, Terrarium-style state containers, and reusable step workspaces.
"""

struct SnowpackPhysicalConstants{NF <: AbstractFloat}
    # Densities (kg/m^3)
    rho_s::NF
    rho_i::NF
    rho_w::NF
    rho_s_a::NF
    rho_s_b::NF
    rho_s_c::NF
    fresh_snow_density_scheme::Symbol

    # Thermal properties
    Ki::NF
    ci::NF
    cw::NF
    Lm::NF

    # Heat flux and albedo parameters
    D_sh::NF
    alpha_dry::NF
    alpha_wet::NF
    alpha_ice::NF
    max_lwc_albedo::NF
    albedo_scheme::Symbol

    # Emissivity
    ϵ_air::NF
    ϵ_snow::NF

    # Universal constants
    σ::NF
    R::NF
    T0::NF
    seconds_per_day::NF
    seconds_per_month::NF
    seconds_per_year::NF
    low_density_densification::Symbol
end

Base.eltype(::SnowpackPhysicalConstants{NF}) where {NF} = NF
@inline number_type(::SnowpackPhysicalConstants{NF}) where {NF} = NF

@inline function _normalize_low_density_densification(scheme::Symbol)
    scheme in (:bessi, :htessel) ||
        error("Unsupported low-density densification scheme '$scheme'. Use :bessi or :htessel.")
    return scheme
end

@inline function _normalize_fresh_snow_density_scheme(scheme::Symbol)
    normalized_scheme = if scheme == :bessi
        :constant
    elseif scheme == :htessel
        :parameterized
    else
        scheme
    end
    normalized_scheme in (:constant, :parameterized) ||
        error(
            "Unsupported fresh-snow density scheme '$scheme'. " *
            "Use :constant, :parameterized, or the aliases :bessi / :htessel.",
        )
    return normalized_scheme
end

@inline function _normalize_albedo_scheme(scheme::Symbol)
    normalized_scheme = if scheme == :bessi
        :constant
    elseif scheme == :legacy
        :constant
    else
        scheme
    end
    normalized_scheme in (:constant, :dynamic) ||
        error(
            "Unsupported albedo scheme '$scheme'. " *
            "Use :constant, :dynamic, or the aliases :legacy / :bessi.",
        )
    return normalized_scheme
end

function SnowpackPhysicalConstants(::Type{NF};
    rho_s::Real=315.0,
    rho_i::Real=917.0,
    rho_w::Real=1000.0,
    rho_s_a::Real=109.0,
    rho_s_b::Real=6.0,
    rho_s_c::Real=26.0,
    fresh_snow_density_scheme::Symbol=:constant,
    Ki::Real=2.1,
    ci::Real=2110.0,
    cw::Real=4181.0,
    Lm::Real=334000.0,
    D_sh::Real=10.0,
    alpha_dry::Real=0.85,
    alpha_wet::Real=0.72,
    alpha_ice::Real=0.3,
    max_lwc_albedo::Real=0.1,
    albedo_scheme::Symbol=:dynamic,
    ϵ_air::Real=0.75,
    ϵ_snow::Real=0.98,
    σ::Real=5.670373e-8,
    R::Real=8.314,
    T0::Real=273.15,
    seconds_per_day::Real=DEFAULT_SECONDS_PER_DAY,
    seconds_per_month::Real=DEFAULT_SECONDS_PER_MONTH,
    seconds_per_year::Real=DEFAULT_SECONDS_PER_YEAR,
    low_density_densification::Symbol=:bessi,
) where {NF <: AbstractFloat}
    return SnowpackPhysicalConstants(
        convert(NF, rho_s),
        convert(NF, rho_i),
        convert(NF, rho_w),
        convert(NF, rho_s_a),
        convert(NF, rho_s_b),
        convert(NF, rho_s_c),
        _normalize_fresh_snow_density_scheme(fresh_snow_density_scheme),
        convert(NF, Ki),
        convert(NF, ci),
        convert(NF, cw),
        convert(NF, Lm),
        convert(NF, D_sh),
        convert(NF, alpha_dry),
        convert(NF, alpha_wet),
        convert(NF, alpha_ice),
        convert(NF, max_lwc_albedo),
        _normalize_albedo_scheme(albedo_scheme),
        convert(NF, ϵ_air),
        convert(NF, ϵ_snow),
        convert(NF, σ),
        convert(NF, R),
        convert(NF, T0),
        convert(NF, seconds_per_day),
        convert(NF, seconds_per_month),
        convert(NF, seconds_per_year),
        _normalize_low_density_densification(low_density_densification),
    )
end

SnowpackPhysicalConstants(; kwargs...) = SnowpackPhysicalConstants(Float64; kwargs...)

mutable struct SnowpackDomain{
        NF <: AbstractFloat,
        NI <: AbstractVector{<:Integer},
        MT <: AbstractMatrix{NF},
        VT <: AbstractVector{NF},
    } <: AbstractSnowpackDomain{NF}
    c::SnowpackPhysicalConstants{NF}
    Ntot::Int
    ncol::Int
    mass_max::NF
    mass_split::NF
    mass_min::NF
    rho_max::NF
    f_base_max::NF
    N::NI
    mass::MT
    mass_w::MT
    density::MT
    temperature::MT
    mass_base::VT
    smb_ice::VT
    runoff::VT
    Tsrf::VT
    snow_cover::VT
    albedo_dynamic::VT
end

struct SnowpackStepForcing{NF <: AbstractFloat}
    air_temperature::NF
    precipitation_rate::NF
    dt_days::NF
    snowfall_rate::NF
    rainfall_rate::NF
    shortwave_down::NF
    wind_speed::NF
    q_sw_net::NF
    q_lw_down::NF
    q_sh::NF
    q_lh::NF
    has_q_sw_net::Bool
    has_q_lw_down::Bool
    has_q_sh::Bool
    has_q_lh::Bool
    diurnal_shortwave::Bool
    latitude::NF
    day_of_year::NF
end

Base.eltype(::AbstractSnowpackState{NF}) where {NF} = NF

function SnowpackDomain(;
    c::SnowpackPhysicalConstants=SnowpackPhysicalConstants(),
    Ntot::Int=DEFAULT_NTOT,
    ncol::Int=1,
    mass_max::Real=DEFAULT_MASS_MAX,
    mass_split::Real=DEFAULT_MASS_SPLIT,
    mass_min::Real=DEFAULT_MASS_MIN,
    rho_max::Real=DEFAULT_RHO_MAX,
    f_base_max::Real=DEFAULT_F_BASE_MAX,
    density_init::Real=DEFAULT_DENSITY_INIT,
    temperature_init::Real=DEFAULT_TEMPERATURE_INIT,
)
    ncol > 0 || error("`ncol` must be positive.")
    mass_split < mass_max || error("`mass_split` must be smaller than `mass_max`.")
    mass_min < mass_split || error("`mass_min` must be smaller than `mass_split`.")
    mass_split / mass_max >= 0.5 || error("`mass_split / mass_max` must be at least 0.5.")

    NF = number_type(c)
    N = zeros(Int, ncol)
    mass = zeros(NF, Ntot, ncol)
    mass_w = zeros(NF, Ntot, ncol)
    density = fill(convert(NF, density_init), Ntot, ncol)
    temperature = fill(convert(NF, temperature_init), Ntot, ncol)
    mass_base = zeros(NF, ncol)
    smb_ice = zeros(NF, ncol)
    runoff = zeros(NF, ncol)
    Tsrf = fill(c.T0, ncol)
    snow_cover = zeros(NF, ncol)
    albedo_dynamic = fill(c.alpha_dry, ncol)
    return SnowpackDomain(
        c,
        Ntot,
        ncol,
        convert(NF, mass_max),
        convert(NF, mass_split),
        convert(NF, mass_min),
        convert(NF, rho_max),
        convert(NF, f_base_max),
        N,
        mass,
        mass_w,
        density,
        temperature,
        mass_base,
        smb_ice,
        runoff,
        Tsrf,
        snow_cover,
        albedo_dynamic,
    )
end

function SnowpackDomain(
    N::AbstractVector{<:Integer},
    mass::AbstractMatrix{NF},
    mass_w::AbstractMatrix{NF},
    density::AbstractMatrix{NF},
    temperature::AbstractMatrix{NF},
    mass_base::AbstractVector{NF},
    smb_ice::AbstractVector{NF},
    runoff::AbstractVector{NF},
    Tsrf::AbstractVector{NF},
    snow_cover::AbstractVector{NF},
    albedo_dynamic::AbstractVector{NF};
    c::SnowpackPhysicalConstants{NF}=SnowpackPhysicalConstants(NF),
    mass_max::Real=DEFAULT_MASS_MAX,
    mass_split::Real=DEFAULT_MASS_SPLIT,
    mass_min::Real=DEFAULT_MASS_MIN,
    rho_max::Real=DEFAULT_RHO_MAX,
    f_base_max::Real=DEFAULT_F_BASE_MAX,
) where {NF <: AbstractFloat}
    ncol = length(N)
    size(mass, 2) == ncol || error("`mass` must have one column per entry of `N`.")
    size(mass_w) == size(mass) || error("`mass_w` must match `mass`.")
    size(density) == size(mass) || error("`density` must match `mass`.")
    size(temperature) == size(mass) || error("`temperature` must match `mass`.")
    length(mass_base) == ncol || error("`mass_base` must match `N`.")
    length(smb_ice) == ncol || error("`smb_ice` must match `N`.")
    length(runoff) == ncol || error("`runoff` must match `N`.")
    length(Tsrf) == ncol || error("`Tsrf` must match `N`.")
    length(snow_cover) == ncol || error("`snow_cover` must match `N`.")
    length(albedo_dynamic) == ncol || error("`albedo_dynamic` must match `N`.")
    mass_split < mass_max || error("`mass_split` must be smaller than `mass_max`.")
    mass_min < mass_split || error("`mass_min` must be smaller than `mass_split`.")
    mass_split / mass_max >= 0.5 || error("`mass_split / mass_max` must be at least 0.5.")

    return SnowpackDomain(
        c,
        size(mass, 1),
        ncol,
        convert(NF, mass_max),
        convert(NF, mass_split),
        convert(NF, mass_min),
        convert(NF, rho_max),
        convert(NF, f_base_max),
        N,
        mass,
        mass_w,
        density,
        temperature,
        mass_base,
        smb_ice,
        runoff,
        Tsrf,
        snow_cover,
        albedo_dynamic,
    )
end

@inline column_count(domain::AbstractSnowpackDomain) = domain.ncol

struct EnergyWorkspace{NF <: AbstractFloat, VT <: AbstractVector{NF}}
    lower::VT
    diag::VT
    upper::VT
    rhs::VT
    interface_conductance::VT
    previous_temperature::VT
    layer_thickness::VT
    thermal_conductivity::VT
end

function EnergyWorkspace(::Type{NF}, Ntot::Int) where {NF <: AbstractFloat}
    allocate() = zeros(NF, Ntot)
    return EnergyWorkspace(
        allocate(),
        allocate(),
        allocate(),
        allocate(),
        allocate(),
        allocate(),
        allocate(),
        allocate(),
    )
end

EnergyWorkspace(domain::AbstractSnowpackDomain{NF}) where {NF <: AbstractFloat} =
    EnergyWorkspace(NF, domain.Ntot)

struct StepWorkspace{NF <: AbstractFloat, VT <: AbstractVector{NF}, ET}
    liquid_water_before_energy::VT
    energy::ET
end

function StepWorkspace(::Type{NF}, Ntot::Int) where {NF <: AbstractFloat}
    return StepWorkspace(zeros(NF, Ntot), EnergyWorkspace(NF, Ntot))
end

StepWorkspace(domain::AbstractSnowpackDomain{NF}) where {NF <: AbstractFloat} =
    StepWorkspace(NF, domain.Ntot)

function threaded_workspaces(domain::AbstractSnowpackDomain{NF}) where {NF <: AbstractFloat}
    return [StepWorkspace(NF, domain.Ntot) for _ in 1:Threads.maxthreadid()]
end

function SnowpackStepForcing(
    c::SnowpackPhysicalConstants,
    air_temperature,
    precipitation_rate,
    dt_days;
    snowfall_rate=zero(air_temperature),
    rainfall_rate=zero(air_temperature),
    shortwave_down=oftype(air_temperature, 400.0),
    wind_speed=oftype(air_temperature, 10.0),
    q_sw_net=nothing,
    q_lw_down=nothing,
    q_sh=nothing,
    q_lh=nothing,
    diurnal_shortwave::Bool=false,
    latitude=zero(air_temperature),
    day_of_year=zero(air_temperature),
)
    NF = number_type(c)
    return SnowpackStepForcing(
        convert(NF, air_temperature),
        convert(NF, precipitation_rate),
        convert(NF, dt_days),
        convert(NF, snowfall_rate),
        convert(NF, rainfall_rate),
        convert(NF, shortwave_down),
        convert(NF, wind_speed),
        convert(NF, isnothing(q_sw_net) ? zero(NF) : q_sw_net),
        convert(NF, isnothing(q_lw_down) ? zero(NF) : q_lw_down),
        convert(NF, isnothing(q_sh) ? zero(NF) : q_sh),
        convert(NF, isnothing(q_lh) ? zero(NF) : q_lh),
        !isnothing(q_sw_net),
        !isnothing(q_lw_down),
        !isnothing(q_sh),
        !isnothing(q_lh),
        diurnal_shortwave,
        convert(NF, latitude),
        convert(NF, day_of_year),
    )
end

function cpu_domain(domain::SnowpackDomain)
    return adapt(Array, domain)
end

function gpu_domain(domain::SnowpackDomain)
    cuda_available() || error("CUDA is not functional in the current environment.")
    return adapt(CUDA.CuArray, domain)
end

variables(::AbstractSnowpackDomain) = (
    PrognosticVariable{:mass}("Solid snow and firn mass per layer."),
    PrognosticVariable{:mass_w}("Liquid water mass per layer."),
    PrognosticVariable{:density}("Bulk snow density per layer."),
    PrognosticVariable{:temperature}("Layer temperature."),
    AuxiliaryVariable{:snow_cover}("Diagnosed snow-cover fraction."),
    AuxiliaryVariable{:albedo_dynamic}("Surface albedo used for radiative forcing."),
)

@adapt_structure SnowpackPhysicalConstants
@adapt_structure SnowpackDomain
@adapt_structure SnowpackStepForcing
@adapt_structure EnergyWorkspace
@adapt_structure StepWorkspace
