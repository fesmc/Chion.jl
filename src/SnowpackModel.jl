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
export get_state
export print_state

export calc_density_gradient_HL80
export calc_density_gradient_powerlaw_ref
export go_densification!

"""
Physical constants for snow/ice model
"""
struct SnowpackPhysicalConstants
    # Densities (kg/m³)
    rho_s::Float64      # Density of fresh snow
    rho_i::Float64      # Density of ice
    rho_w::Float64      # Density of water
    
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
end

"""
    SnowpackPhysicalConstants(; kwargs...)

Initialize physical constants with default or custom values.

# Keyword Arguments
- `D_sh`: Coefficient for sensible heat flux, default=10 W/(m²·K), range=[5, 20]
- `alpha_dry`: Albedo of fresh snow, default=0.8, range=[0.75, 0.9]
- `alpha_wet`: Albedo of wet snow, default=0.6, range=[0.5, 0.7]
- `ϵ_air`: Emissivity of air, default=0.75, range=[0.6, 0.9]

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
    rho_s::Float64=150.0,
    rho_i::Float64=917.0,
    rho_w::Float64=1000.0,
    
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
    seconds_per_year::Float64 = DEFAULT_SECONDS_PER_YEAR
)
    return SnowpackPhysicalConstants(
        # Densities
        rho_s,
        rho_i,
        rho_w,
        
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
        seconds_per_year
    )
end

"""
    SnowpackColumn

A column-based snowpack model with mass-following dynamic grid.

# Grid parameters
- `Ntot::Int`: Maximum number of vertical layers (default: 15)
- `N::Int`: Number of currently active (filled) layers
- `kbase::Int`: Index of base active layer (Ntot-N+1)

# Parameters (from Born et al. 2019)
- `mass_max::Float64`: Maximum mass before layer split [kg/m²] (default: 500)
- `mass_split::Float64`: Target mass for split layers [kg/m²] (default: 300)
- `mass_min::Float64`: Minimum mass before layer merge [kg/m²] (default: 100)
- `rho_i::Float64`: Ice density [kg/m³] (default: 917)

# State variables
- `mass::Vector{Float64}`: Mass of snow+water in each layer [kg/m²]
- `mass_snow::Vector{Float64}`: Mass of snow in each layer [kg/m²]
- `mass_w::Vector{Float64}`: Mass of water in each layer [kg/m²]
- `density::Vector{Float64}`: Density of snow in each layer [kg/m³]

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
    mass::Vector{Float64}           # kg/m²
    mass_w::Vector{Float64}         # kg/m²
    density::Vector{Float64}        # kg/m³
    temperature::Vector{Float64}    # K
    mass_base::Float64              # kg/m²
    runoff::Float64                 # kg/m²
    Tsrf::Float64                   # K

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

        # Initialize with no initial mass
        mass = zeros(Float64, Ntot)
        mass_w = zeros(Float64, Ntot)
        density = fill(density_init, Ntot)
        temperature = fill(temperature_init, Ntot)
        mass_base = 0.0
        runoff = 0.0
        Tsrf = c.T0

        # Consistency check
        @assert mass_split < mass_max
        @assert mass_min < mass_split

        # Make sure mass_split is more than 50% of mass_max, so that when
        # surface layer splits, the surface contains less mass than the subsurface layer
        @assert mass_split / mass_max >= 0.5
        
        new(c, Ntot, N, mass_max, mass_split, mass_min, rho_max, f_base_max,
            mass, mass_w, density, temperature, mass_base, runoff, Tsrf)
    end
end

include("energy_flux.jl")
include("densification.jl")
include("mass_balance.jl")
include("percolation.jl")
include("refreezing.jl")


"""
    step!(column::SnowpackColumn, mdot::Float64, dt::Float64) -> Float64

Advance the snowpack column by one time step.

# Arguments
- `column`: The snowpack column to update
- `T2m` : Near-surface air temperature [K]
- `P`: Precipitation rate at surface [kg/m²/s]
- `dt`: Time step [d]
- `f_s`: Fraction of precipitation that is snow [1], default nothing, calculate internally
- `p_snow`: Optional direct snowfall rate [kg/m²/s], overrides `P/f_s` partition when provided
- `p_rain`: Optional direct rainfall rate [kg/m²/s], overrides `P/f_s` partition when provided
- `s_boa`: Optional downward shortwave forcing [W m^-2], used when `q_sw_net` is not provided

# Process
1. Apply surface mass flux
2. Handle layer splitting/merging
3. Propagate melt through layers if negative
4. Check for ice formation at base
"""
function step!(
    column::SnowpackColumn,
    T2m::Float64,
    P::Float64,
    dt::Float64;
    f_s=nothing,
    P_ave=P,
    p_snow::Union{Nothing, Float64}=nothing,
    p_rain::Union{Nothing, Float64}=nothing,
    s_boa::Union{Nothing, Float64}=nothing,
    q_sw_net::Union{Nothing, Float64}=nothing,
    q_lw_down::Union{Nothing, Float64}=nothing,
    q_sh::Union{Nothing, Float64}=nothing,
    q_lh::Union{Nothing, Float64}=nothing,
)
    if !isnothing(p_snow) || !isnothing(p_rain)
        # Direct forcing path: caller provides separated rain/snow rates.
        P_snow = isnothing(p_snow) ? 0.0 : p_snow
        P_rain = isnothing(p_rain) ? 0.0 : p_rain
    else
        if isnothing(f_s)
            # Determine fraction of snow and rain as a function of T2m
            # following Born et al. (2019)
            if T2m > column.c.T0
                f_s = 0.0
            else
                f_s = 1.0
            end
        end

        # Backward-compatible partitioning from total precipitation and f_s.
        P_rain = P * (1.0-f_s)
        P_snow = P - P_rain
    end

    # Convert timestep to seconds internally
    dt_sec = dt * column.c.seconds_per_day

    # For first snowfall, seed surface temperature with air temperature
    # (Fortran behavior when first box is empty and snow starts).
    if P_snow > 0.0 && column.N > 0 && column.mass[1] == 0.0
        column.temperature[1] = T2m
    end

    # Handle accumulation first
    apply_accumulation!(column, P_snow, P_rain, dt_sec)

    # Fortran passes At = accum + rainman [kg m^-2 s^-1] to densification,
    # and only runs densification when at least 3 boxes are snow-filled.
    At = (P_snow > 0.0 ? P_snow : 0.0) + ((column.N > 0 && column.mass[1] > 0.0) ? P_rain : 0.0)
    if column.N >= 3 && column.mass[3] > 0.0
        go_densification!(column, At, dt_sec)
    end

    # Caculate energy balance
    S_boa = isnothing(s_boa) ? 400.0 : max(s_boa, 0.0)
    energy = go_energy_flux!(
        column, T2m, S_boa, nothing, nothing, dt_sec;
        P_snow=P_snow,
        P_rain=P_rain,
        diff_model=1,
        q_sw_net=q_sw_net,
        q_lw_down=q_lw_down,
        q_sh=q_sh,
        q_lh=q_lh,
    )

    # For now set a linear temperature profile in the firn to depth
    #column.Tsrf = min(T2m,column.c.T0)
    #column.temperature[1] = column.Tsrf
    #column.temperature[column.N] = column.Tsrf - 10.0
    # Handle melt
    if energy.china_syndrome
        Ts = column.temperature[1]
        if isnothing(q_sw_net) && isnothing(q_lw_down) && isnothing(q_sh) && isnothing(q_lh)
            # Backward-compatible melt energy diagnosis for default parameterized forcing.
            Qp_lw = column.c.σ * (column.c.ϵ_air * T2m^4 - column.c.ϵ_snow * Ts^4)
            Qp_sh = column.c.D_sh * (T2m - Ts)
            Qp_lh = energy.K_lh - energy.H_lh * Ts
            QQ = max((S_boa + Qp_lw + Qp_sh + Qp_lh) * dt_sec - energy.Q_heat, 0.0)
        else
            # For externally prescribed fluxes, use the same linearized net-flux form
            # as in the temperature solve.
            QQ = max((energy.F_const - energy.F_lin * Ts) * dt_sec - energy.Q_heat, 0.0)
        end
        melt_mass = QQ / column.c.Lm
        apply_melt!(column, melt_mass)
    end

    # Fortran flow: melting -> percolation -> refreezing.
    # Apply percolation first if liquid water exists.
    if column.N > 0 && maximum(@view column.mass_w[1:column.N]) > 0.0
        go_percolation!(column)
    end

    # Then refreeze liquid water into the cold content of each active layer.
    if column.N > 0 && maximum(@view column.mass_w[1:column.N]) > 0.0
        go_refreezing!(column)
    end

    
    return
end



"""
    get_state(column::SnowpackColumn) -> Dict

Get the current state of the snowpack column.

Returns a dictionary with:
- `N`: Number of active layers
- `mass`: Mass in each active layer [kg/m²]
- `density`: Density in each active layer [kg/m³]
- `total_mass`: Total mass in column [kg/m²]
- `thickness`: Thickness of each active layer [m]
- `total_thickness`: Total column thickness [m]
"""
function get_state(column::SnowpackColumn)
    if column.N == 0
        return Dict(
            "N" => 0,
            "mass" => Float64[],
            "density" => Float64[],
            "total_mass" => 0.0,
            "thickness" => Float64[],
            "total_thickness" => 0.0
        )
    end
    
    active_mass = column.mass[1:column.N]
    active_density = column.density[1:column.N]
    thickness = active_mass ./ active_density
    
    return Dict(
        "N" => column.N,
        "mass" => active_mass,
        "density" => active_density,
        "total_mass" => sum(active_mass),
        "thickness" => thickness,
        "total_thickness" => sum(thickness)
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
