"""
Column-based snowpack model with dynamic layering.
Based on Born et al. (2019) algorithm.

This initial version focuses on mass conservation and layer dynamics.
"""

module SnowpackModel

using Printf
include("model_constants.jl")
export SnowpackPhysicalConstants
export SnowpackColumn
export step!
export go_percolation!
export go_refreezing!
export continuous_bottom_deplete!
export get_state
export print_state
export go_densification!

"""
Physical constants for snow/ice model
"""
struct SnowpackPhysicalConstants
    # Densities (kg/m³)
    rho_s::Float64      # Legacy constant fresh-snow density [kg/m³]
    rho_i::Float64      # Density of ice
    rho_w::Float64      # Density of water
    rho_s_a::Float64    # Fresh-snow density parameter a [kg/m³]
    rho_s_b::Float64    # Fresh-snow density parameter b [kg/m³/K]
    rho_s_c::Float64    # Fresh-snow density parameter c [kg/m³/(m/s)^0.5]
    
    # Thermal properties
    Ki::Float64         # Thermal conductivity of ice (W/(m·K))
    ci::Float64         # Heat capacity of ice (J/(kg·K))
    cw::Float64         # Heat capacity of water (J/(kg·K))
    Lm::Float64         # Latent heat of melting (J/kg)
    
    # Heat flux and albedo parameters
    D_sh::Float64       # Coefficient for sensible heat flux (W/(m²·K))
    alpha_dry::Float64  # Albedo of fresh snow
    alpha_wet::Float64  # Albedo of wet snow
    alpha_ice::Float64  # Albedo of ice
    
    # Emissivity
    ϵ_air::Float64      # Emissivity of air
    ϵ_snow::Float64     # Emissivity of snow
    
    # Universal constants
    σ::Float64                  # Stefan-Boltzmann constant (W/(m²·K⁴))
    R::Float64                  # Universal gas constant (J/(K·mol))
    T0::Float64                 # Freezing point of water (K)
    seconds_per_day::Float64    # Seconds per day
    seconds_per_month::Float64  # Seconds per month
    seconds_per_year::Float64   # Seconds per year
    low_density_densification::Symbol # :bessi or :htessel for rho < 550 kg m^-3
end

@inline function _normalize_low_density_densification(scheme::Symbol)
    scheme in (:bessi, :htessel) ||
        error("Unsupported low-density densification scheme '$scheme'. Use :bessi or :htessel.")
    return scheme
end

"""
    SnowpackPhysicalConstants(; kwargs...)

Initialize physical constants with default or custom values.

# Keyword Arguments
- `D_sh`: Coefficient for sensible heat flux, default=10 W/(m²·K), range=[5, 20]
- `alpha_dry`: Albedo of fresh snow, default=0.8, range=[0.75, 0.9]
- `alpha_wet`: Albedo of wet snow, default=0.6, range=[0.5, 0.7]
- `ϵ_air`: Emissivity of air, default=0.75, range=[0.6, 0.9]
- `rho_s_a`: Fresh-snow density parameter `a`, default=109
- `rho_s_b`: Fresh-snow density parameter `b`, default=6
- `rho_s_c`: Fresh-snow density parameter `c`, default=26
- `low_density_densification`: Scheme for `rho < 550 kg m^-3`, one of `:bessi` or `:htessel`

# Example
```julia
# Use default values
c = SnowpackPhysicalConstants()

# Customize specific parameters
c = SnowpackPhysicalConstants(D_sh=20.0, alpha_dry=0.85, ϵ_air=0.8)
```
"""
function SnowpackPhysicalConstants(;
    # Densities (kg/m³)
    rho_s::Float64=250.0,
    rho_i::Float64=917.0,
    rho_w::Float64=1000.0,
    rho_s_a::Float64=109.0,
    rho_s_b::Float64=6.0,
    rho_s_c::Float64=26.0,
    
    # Thermal properties
    Ki::Float64=2.1,
    ci::Float64=2110.0,
    cw::Float64=4181.0,
    Lm::Float64=334000.0,
    
    # Heat flux and albedo
    D_sh::Float64=10.0,
    alpha_dry::Float64=0.8,
    alpha_wet::Float64=0.6,
    alpha_ice::Float64=0.35,
    
    # Emissivity
    ϵ_air::Float64=0.75,
    ϵ_snow::Float64=0.98,
    
    # Universal constants
    σ::Float64=5.670373e-8,
    R::Float64=8.314,
    T0::Float64=273.15,
    seconds_per_day::Float64 = DEFAULT_SECONDS_PER_DAY,
    seconds_per_month::Float64 = DEFAULT_SECONDS_PER_MONTH,
    seconds_per_year::Float64 = DEFAULT_SECONDS_PER_YEAR,
    low_density_densification::Symbol=:bessi,
)
    return SnowpackPhysicalConstants(
        # Densities
        rho_s,
        rho_i,
        rho_w,
        rho_s_a,
        rho_s_b,
        rho_s_c,
        
        # Thermal properties
        Ki,
        ci,
        cw,
        Lm,
        
        # Heat flux and albedo
        D_sh,
        alpha_dry,
        alpha_wet,
        alpha_ice,
        
        # Emissivity
        ϵ_air,
        ϵ_snow,
        
        # Universal constants
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
- `mass_max::Float64`: Maximum mass before layer split [kg/m²] (default: 500)
- `mass_split::Float64`: Target mass for split layers [kg/m²] (default: 300)
- `mass_min::Float64`: Minimum mass before layer merge [kg/m²] (default: 100)
- `rho_i::Float64`: Ice density [kg/m³] (default: 917)

# State variables
- `mass::Vector{Float64}`: Solid snow/ice mass in each layer [kg/m²]
- `mass_w::Vector{Float64}`: Liquid water mass in each layer [kg/m²]
- `density::Vector{Float64}`: Bulk snow density in each layer [kg/m³]

"""
mutable struct SnowpackColumn
    # Constants
    c::SnowpackPhysicalConstants

    # Grid parameters
    Ntot::Int
    N::Int

    # Model parameters
    mass_max::Float64           # kg/m²
    mass_split::Float64         # kg/m²
    mass_min::Float64           # kg/m²
    rho_max::Float64            # kg/m³
    f_base_max::Float64         # 1

    #ζmax::Float64   # Maximum liquid water content

    # State variables
    mass::Vector{Float64}           # Solid snow/ice mass [kg/m²]
    mass_w::Vector{Float64}         # Liquid water mass [kg/m²]
    density::Vector{Float64}        # kg/m³
    temperature::Vector{Float64}    # K
    mass_base::Float64              # kg/m²
    runoff::Float64                 # kg/m²
    Tsrf::Float64                   # Legacy surface temperature state [K]
    snow_cover::Float64             # 1

    function SnowpackColumn(;
        c::SnowpackPhysicalConstants = SnowpackPhysicalConstants(),
        Ntot::Int = DEFAULT_NTOT,
        N::Int = DEFAULT_N_ACTIVE,
        mass_max::Float64 = DEFAULT_MASS_MAX,
        mass_split::Float64 = DEFAULT_MASS_SPLIT,
        mass_min::Float64 = DEFAULT_MASS_MIN,
        rho_max::Float64 = DEFAULT_RHO_MAX,
        f_base_max::Float64 = DEFAULT_F_BASE_MAX,
        density_init::Float64 = DEFAULT_DENSITY_INIT,
        temperature_init::Float64 = DEFAULT_TEMPERATURE_INIT,
    )   

        @assert 0 <= N <= Ntot

        # Initialize with no initial mass
        mass = zeros(Float64, Ntot)
        mass_w = zeros(Float64, Ntot)
        density = fill(density_init, Ntot)
        temperature = fill(temperature_init, Ntot)
        mass_base = 0.0
        runoff = 0.0
        Tsrf = c.T0
        snow_cover = 0.0

        # Consistency check
        @assert mass_split < mass_max
        @assert mass_min < mass_split

        # Make sure mass_split is more than 50% of mass_max, so that when
        # surface layer splits, the surface contains less mass than the subsurface layer
        @assert mass_split / mass_max >= 0.5
        
        new(c, Ntot, N, mass_max, mass_split, mass_min, rho_max, f_base_max,
            mass, mass_w, density, temperature, mass_base, runoff, Tsrf, snow_cover)
    end
end

include("energy_flux.jl")
include("densification.jl")
include("mass_balance.jl")
include("percolation.jl")
include("refreezing.jl")

@inline function _resolve_keyword_alias(
    preferred_value,
    legacy_value,
    preferred_name::AbstractString,
    legacy_name::AbstractString,
)
    if !isnothing(preferred_value) && !isnothing(legacy_value) &&
       !isequal(preferred_value, legacy_value)
        error(
            "Received both `$preferred_name` and `$legacy_name` with different values. " *
            "Use one or provide matching values.",
        )
    end
    return isnothing(preferred_value) ? legacy_value : preferred_value
end

@inline function _bulk_snow_density(column::SnowpackColumn)
    if column.N <= 0
        return 0.0
    end

    total_mass = 0.0
    total_thickness = 0.0
    @inbounds for i in 1:column.N
        m = column.mass[i]
        ρ = column.density[i]
        if m > 0.0 && ρ > EPS_TINY
            total_mass += m
            total_thickness += m / ρ
        end
    end

    if total_mass <= 0.0 || total_thickness <= EPS_TINY
        return 0.0
    end
    return total_mass / total_thickness
end

@inline function _total_snow_water_mass(column::SnowpackColumn)
    if column.N <= 0
        return 0.0
    end

    total_wet_mass = 0.0
    @inbounds for i in 1:column.N
        total_wet_mass += max(column.mass[i], 0.0) + max(column.mass_w[i], 0.0)
    end
    return total_wet_mass
end

@inline function _snow_cover_fraction(column::SnowpackColumn)
    if column.N <= 0
        return 0.0
    end

    W_s = _total_snow_water_mass(column)
    W_s <= 0.0 && return 0.0

    rho_sn = _bulk_snow_density(column)
    rho_sn <= EPS_TINY && return 0.0

    return min(1.0, (W_s / rho_sn) / 0.1)
end

@inline function update_snow_cover!(column::SnowpackColumn)
    column.snow_cover = _snow_cover_fraction(column)
    return column.snow_cover
end

@inline function _column_has_liquid_water(column::SnowpackColumn)
    @inbounds for layer_index in 1:column.N
        if column.mass_w[layer_index] > EPS_TINY
            return true
        end
    end
    return false
end


"""
    step!(column::SnowpackColumn, air_temperature::Float64, precipitation_rate::Float64, dt_days::Float64)

Advance the snowpack column by one time step.

# Arguments
- `column`: The snowpack column to update
- `air_temperature`: Near-surface air temperature [K]
- `precipitation_rate`: Total precipitation rate at surface [kg/m²/s]
- `dt_days`: Time step [d]
- `snow_fraction`: Optional fraction of precipitation that falls as snow [1]
- `snowfall_rate`: Optional direct snowfall rate [kg/m²/s]
- `rainfall_rate`: Optional direct rainfall rate [kg/m²/s]
- `shortwave_down`: Optional downward shortwave forcing [W m^-2], used when `q_sw_net` is not provided
- `wind_speed`: Optional near-surface wind speed [m/s], default = `5.0`

Legacy keyword aliases `f_s`, `p_snow`, `p_rain`, and `s_boa` are still accepted.

# Process
1. Apply surface mass flux
2. Handle layer splitting/merging
3. Solve temperature evolution
4. Apply melt, percolation, and refreezing when needed
"""
function step!(
    column::SnowpackColumn,
    air_temperature::Float64,
    precipitation_rate::Float64,
    dt_days::Float64;
    snow_fraction=nothing,
    f_s=nothing,
    P_ave=precipitation_rate,
    snowfall_rate::Union{Nothing, Float64}=nothing,
    rainfall_rate::Union{Nothing, Float64}=nothing,
    shortwave_down::Union{Nothing, Float64}=nothing,
    p_snow::Union{Nothing, Float64}=nothing,
    p_rain::Union{Nothing, Float64}=nothing,
    s_boa::Union{Nothing, Float64}=nothing,
    wind_speed::Float64=5.0,
    q_sw_net::Union{Nothing, Float64}=nothing,
    q_lw_down::Union{Nothing, Float64}=nothing,
    q_sh::Union{Nothing, Float64}=nothing,
    q_lh::Union{Nothing, Float64}=nothing,
)
    resolved_snow_fraction = _resolve_keyword_alias(snow_fraction, f_s, "snow_fraction", "f_s")
    resolved_snowfall_rate = _resolve_keyword_alias(snowfall_rate, p_snow, "snowfall_rate", "p_snow")
    resolved_rainfall_rate = _resolve_keyword_alias(rainfall_rate, p_rain, "rainfall_rate", "p_rain")
    resolved_shortwave_down = _resolve_keyword_alias(shortwave_down, s_boa, "shortwave_down", "s_boa")

    if !isnothing(resolved_snowfall_rate) || !isnothing(resolved_rainfall_rate)
        # Direct forcing path: caller provides separated rain/snow rates.
        snowfall_rate = isnothing(resolved_snowfall_rate) ? 0.0 : resolved_snowfall_rate
        rainfall_rate = isnothing(resolved_rainfall_rate) ? 0.0 : resolved_rainfall_rate
    else
        if isnothing(resolved_snow_fraction)
            # Determine the snowfall fraction following Born et al. (2019).
            if air_temperature > column.c.T0
                resolved_snow_fraction = 0.0
            else
                resolved_snow_fraction = 1.0
            end
        end

        rainfall_rate = precipitation_rate * (1.0 - resolved_snow_fraction)
        snowfall_rate = precipitation_rate - rainfall_rate
    end

    # Convert timestep to seconds internally
    dt_seconds = dt_days * column.c.seconds_per_day
    started_without_surface_snow = column.N == 0 || column.mass[1] <= EPS_EMPTY_LAYER

    apply_accumulation!(
        column,
        snowfall_rate,
        rainfall_rate,
        dt_seconds;
        air_temperature=air_temperature,
        wind_speed=wind_speed,
    )
    if snowfall_rate > 0.0 && started_without_surface_snow && column.N > 0
        column.temperature[1] = air_temperature
    end
    update_snow_cover!(column)
    liquid_water_before_energy = column.c.low_density_densification == :htessel ?
        copy(@view column.mass_w[1:column.N]) : Float64[]
    accumulation_rate = (snowfall_rate > 0.0 ? snowfall_rate : 0.0) +
                        ((column.N > 0 && column.mass[1] > 0.0) ? rainfall_rate : 0.0)
    if column.N >= 1 && column.mass[1] > 0.0
        go_densification!(column, accumulation_rate, dt_seconds)
    end

    # Calculate the energy balance.
    diagnosed_shortwave_down = isnothing(resolved_shortwave_down) ? 400.0 : max(resolved_shortwave_down, 0.0)
    energy = go_energy_flux!(
        column,
        air_temperature,
        diagnosed_shortwave_down,
        nothing,
        nothing,
        dt_seconds;
        snowfall_rate=snowfall_rate,
        rainfall_rate=rainfall_rate,
        diffusion_model=1,
        q_sw_net=q_sw_net,
        q_lw_down=q_lw_down,
        q_sh=q_sh,
        q_lh=q_lh,
    )

    if energy.needs_melt
        surface_temperature = column.temperature[1]
        if isnothing(q_sw_net) && isnothing(q_lw_down) && isnothing(q_sh) && isnothing(q_lh)
            # Backward-compatible melt energy diagnosis for default parameterized forcing.
            longwave_flux = column.c.σ * (
                column.c.ϵ_air * air_temperature^4 - column.c.ϵ_snow * surface_temperature^4
            )
            sensible_heat_flux = column.c.D_sh * (air_temperature - surface_temperature)
            latent_heat_flux = energy.latent_heat_constant_term -
                               energy.latent_heat_linear_coefficient * surface_temperature
            melt_energy = max(
                (diagnosed_shortwave_down + longwave_flux + sensible_heat_flux + latent_heat_flux) *
                dt_seconds - energy.energy_to_melting,
                0.0,
            )
        else
            # For externally prescribed fluxes, use the same linearized net-flux form
            # as in the temperature solve.
            melt_energy = max(
                (energy.surface_flux_constant - energy.surface_flux_linear * surface_temperature) *
                dt_seconds - energy.energy_to_melting,
                0.0,
            )
        end
        melt_mass = melt_energy / column.c.Lm
        apply_melt!(column, melt_mass)
    end

    has_liquid_water = _column_has_liquid_water(column)
    if has_liquid_water
        go_percolation!(column)
        has_liquid_water = _column_has_liquid_water(column)
    end

    if column.c.low_density_densification == :htessel && !isempty(liquid_water_before_energy) &&
       has_liquid_water
        _apply_htessel_liquid_water_compaction!(column, liquid_water_before_energy, dt_seconds)
        has_liquid_water = _column_has_liquid_water(column)
    end

    if has_liquid_water
        go_refreezing!(column)
    end

    update_snow_cover!(column)
    return nothing
end



"""
    get_state(column::SnowpackColumn) -> Dict

Get the current state of the snowpack column.

Returns a dictionary with:
- `N`: Number of active layers
- `mass` / `solid_mass`: Solid snow/ice mass in each active layer [kg/m²]
- `mass_w` / `liquid_water_mass`: Liquid water mass in each active layer [kg/m²]
- `density`: Density in each active layer [kg/m³]
- `total_mass`: Total mass in column [kg/m²]
- `total_liquid_water`: Total liquid water mass in column [kg/m²]
- `total_wet_mass`: Total snow plus liquid water mass in column [kg/m²]
- `thickness`: Thickness of each active layer [m]
- `total_thickness`: Total column thickness [m]
- `surface_temperature`: Surface temperature [K]
- `snow_cover`: Diagnosed snow cover fraction [1]
"""
function get_state(column::SnowpackColumn)
    snow_cover = _snow_cover_fraction(column)
    if column.N == 0
        return Dict(
            "N" => 0,
            "n_active" => 0,
            "mass" => Float64[],
            "solid_mass" => Float64[],
            "mass_w" => Float64[],
            "liquid_water_mass" => Float64[],
            "density" => Float64[],
            "total_mass" => 0.0,
            "total_liquid_water" => 0.0,
            "total_wet_mass" => 0.0,
            "thickness" => Float64[],
            "total_thickness" => 0.0,
            "surface_temperature" => column.c.T0,
            "snow_cover" => snow_cover,
        )
    end

    @views active_solid_mass = column.mass[1:column.N]
    @views active_liquid_water_mass = column.mass_w[1:column.N]
    @views active_density = column.density[1:column.N]
    thickness = active_solid_mass ./ active_density

    return Dict(
        "N" => column.N,
        "n_active" => column.N,
        "mass" => active_solid_mass,
        "solid_mass" => active_solid_mass,
        "mass_w" => active_liquid_water_mass,
        "liquid_water_mass" => active_liquid_water_mass,
        "density" => active_density,
        "total_mass" => sum(active_solid_mass),
        "total_liquid_water" => sum(active_liquid_water_mass),
        "total_wet_mass" => sum(active_solid_mass) + sum(active_liquid_water_mass),
        "thickness" => thickness,
        "total_thickness" => sum(thickness),
        "surface_temperature" => column.temperature[1],
        "snow_cover" => snow_cover,
    )
end


"""
    print_state(column::SnowpackColumn)

Print a formatted summary of the current snowpack state.
"""
function print_state(column::SnowpackColumn)
    state = get_state(column)
    
    println("=" ^ 60)
    println("Snowpack Column State")
    println("=" ^ 60)
    println("Active layers: ", state["N"])
    println("Total mass: ", round(state["total_mass"], digits=2), " kg/m²")
    println("Total thickness: ", round(state["total_thickness"], digits=3), " m")
    println("Snow cover: ", round(state["snow_cover"], digits=3))
    println()
    
    if state["N"] > 0
        println("Layer details (surface = 1):")
        println("-" ^ 60)
        println("Layer | Mass (kg/m²) | Density (kg/m³) | Thickness (m)")
        println("-" ^ 60)
        for i in 1:state["N"]
            @printf("%5d | %12.2f | %15.1f | %13.4f\n", 
                    i, state["mass"][i], state["density"][i], state["thickness"][i])
        end
        println("=" ^ 60)
    end
end

end # module
