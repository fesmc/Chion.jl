"""
Core snowpack types and constructors.
"""

"""
Physical constants for snow/ice model.
"""
struct SnowpackPhysicalConstants
    # Densities (kg/m^3)
    rho_s::Float64
    rho_i::Float64
    rho_w::Float64
    rho_s_a::Float64
    rho_s_b::Float64
    rho_s_c::Float64
    fresh_snow_density_scheme::Symbol

    # Thermal properties
    Ki::Float64
    ci::Float64
    cw::Float64
    Lm::Float64

    # Heat flux and albedo parameters
    D_sh::Float64
    alpha_dry::Float64
    alpha_wet::Float64
    alpha_ice::Float64
    max_lwc_albedo::Float64
    albedo_scheme::Symbol

    # Emissivity
    ϵ_air::Float64
    ϵ_snow::Float64

    # Universal constants
    σ::Float64
    R::Float64
    T0::Float64
    seconds_per_day::Float64
    seconds_per_month::Float64
    seconds_per_year::Float64
    low_density_densification::Symbol
end

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
        :dynamic
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

"""
    SnowpackPhysicalConstants(; kwargs...)

Initialize physical constants with default or custom values.

# Keyword Arguments
- `D_sh`: Coefficient for sensible heat flux, default=10 W/(m^2 K), range=[5, 20]
- `alpha_dry`: Albedo of fresh snow, default=0.85, range=[0.75, 0.9]
- `alpha_wet`: Albedo of wet snow, default=0.72, range=[0.5, 0.8]
- `max_lwc_albedo`: Liquid-water-content scale used by the BESSI Aoki albedo
  reduction, default=`0.1`
- `albedo_scheme`: `:dynamic` for the BESSI-style persistent Aoki albedo, or
  `:constant` for the legacy dry-snow / wet-snow / ice switch
- `ϵ_air`: Emissivity of air, default=0.75, range=[0.6, 0.9]
- `rho_s_a`: Fresh-snow density parameter `a`, default=109
- `rho_s_b`: Fresh-snow density parameter `b`, default=6
- `rho_s_c`: Fresh-snow density parameter `c`, default=26
- `fresh_snow_density_scheme`: `:constant` for the original BESSI-style constant
  `rho_s`, or `:parameterized` for `a + b*(T_air - T0) + c*sqrt(wind)`
- `low_density_densification`: Scheme for `rho < 550 kg m^-3`, one of
  `:bessi` or `:htessel`
"""
function SnowpackPhysicalConstants(;
    # Densities (kg/m^3)
    rho_s::Float64=315.0,
    rho_i::Float64=917.0,
    rho_w::Float64=1000.0,
    rho_s_a::Float64=109.0,
    rho_s_b::Float64=6.0,
    rho_s_c::Float64=26.0,
    fresh_snow_density_scheme::Symbol=:constant,

    # Thermal properties
    Ki::Float64=2.1,
    ci::Float64=2110.0,
    cw::Float64=4181.0,
    Lm::Float64=334000.0,

    # Heat flux and albedo
    D_sh::Float64=10.0,
    alpha_dry::Float64=0.85,
    alpha_wet::Float64=0.72,
    alpha_ice::Float64=0.3,
    max_lwc_albedo::Float64=0.1,
    albedo_scheme::Symbol=:constant,

    # Emissivity
    ϵ_air::Float64=0.75,
    ϵ_snow::Float64=0.98,

    # Universal constants
    σ::Float64=5.670373e-8,
    R::Float64=8.314,
    T0::Float64=273.15,
    seconds_per_day::Float64=DEFAULT_SECONDS_PER_DAY,
    seconds_per_month::Float64=DEFAULT_SECONDS_PER_MONTH,
    seconds_per_year::Float64=DEFAULT_SECONDS_PER_YEAR,
    low_density_densification::Symbol=:htessel,
)
    return SnowpackPhysicalConstants(
        rho_s,
        rho_i,
        rho_w,
        rho_s_a,
        rho_s_b,
        rho_s_c,
        _normalize_fresh_snow_density_scheme(fresh_snow_density_scheme),
        Ki,
        ci,
        cw,
        Lm,
        D_sh,
        alpha_dry,
        alpha_wet,
        alpha_ice,
        max_lwc_albedo,
        _normalize_albedo_scheme(albedo_scheme),
        ϵ_air,
        ϵ_snow,
        σ,
        R,
        T0,
        seconds_per_day,
        seconds_per_month,
        seconds_per_year,
        _normalize_low_density_densification(low_density_densification),
    )
end

"""
    SnowpackColumn

A column-based snowpack model with mass-following dynamic grid.

# Grid parameters
- `Ntot::Int`: Maximum number of vertical layers (default: 7)
- `N::Int`: Number of currently active layers

# Parameters (from Born et al. 2019)
- `mass_max::Float64`: Maximum mass before layer split [kg/m^2] (default: 500)
- `mass_split::Float64`: Target mass for split layers [kg/m^2] (default: 300)
- `mass_min::Float64`: Minimum mass before layer merge [kg/m^2] (default: 100)
- `rho_i::Float64`: Ice density [kg/m^3] (default: 917)

# State variables
- `mass::Vector{Float64}`: Solid snow/ice mass in each layer [kg/m^2]
- `mass_w::Vector{Float64}`: Liquid water mass in each layer [kg/m^2]
- `density::Vector{Float64}`: Bulk snow density in each layer [kg/m^3]
"""
mutable struct SnowpackColumn
    # Constants
    c::SnowpackPhysicalConstants

    # Grid parameters
    Ntot::Int
    N::Int

    # Model parameters
    mass_max::Float64
    mass_split::Float64
    mass_min::Float64
    rho_max::Float64
    f_base_max::Float64

    # State variables
    mass::Vector{Float64}
    mass_w::Vector{Float64}
    density::Vector{Float64}
    temperature::Vector{Float64}
    mass_base::Float64
    smb_ice::Float64
    runoff::Float64
    Tsrf::Float64
    snow_cover::Float64
    albedo_dynamic::Float64

    function SnowpackColumn(;
        c::SnowpackPhysicalConstants=SnowpackPhysicalConstants(),
        Ntot::Int=DEFAULT_NTOT,
        N::Int=DEFAULT_N_ACTIVE,
        mass_max::Float64=DEFAULT_MASS_MAX,
        mass_split::Float64=DEFAULT_MASS_SPLIT,
        mass_min::Float64=DEFAULT_MASS_MIN,
        rho_max::Float64=DEFAULT_RHO_MAX,
        f_base_max::Float64=DEFAULT_F_BASE_MAX,
        density_init::Float64=DEFAULT_DENSITY_INIT,
        temperature_init::Float64=DEFAULT_TEMPERATURE_INIT,
    )
        @assert 0 <= N <= Ntot

        mass = zeros(Float64, Ntot)
        mass_w = zeros(Float64, Ntot)
        density = fill(density_init, Ntot)
        temperature = fill(temperature_init, Ntot)
        mass_base = 0.0
        smb_ice = 0.0
        runoff = 0.0
        Tsrf = c.T0
        snow_cover = 0.0
        albedo_dynamic = c.alpha_dry

        @assert mass_split < mass_max
        @assert mass_min < mass_split

        # Keep the split target in the upper half of the admissible range so the
        # split creates a lighter surface layer than the subsurface remainder.
        @assert mass_split / mass_max >= 0.5

        new(
            c,
            Ntot,
            N,
            mass_max,
            mass_split,
            mass_min,
            rho_max,
            f_base_max,
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
end
