"""
Column-based snowpack model with dynamic layering.
Based on Born et al. (2019) algorithm.

This initial version focuses on mass conservation and layer dynamics.
"""

module SnowpackModel

using Printf

export SnowpackPhysicalConstants
export SnowpackColumn
export step!
export get_state
export print_state

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
    rho_s::Float64=350.0,
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
    seconds_per_day::Float64 = 86400.0
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
        seconds_per_day
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
    
    function SnowpackColumn(;
        c::SnowpackPhysicalConstants = SnowpackPhysicalConstants(),
        Ntot::Int = 7,
        N::Int = 1,
        mass_max::Float64 = 500.0,
        mass_split::Float64 = 300.0,
        mass_min::Float64 = 100.0,
        rho_max::Float64 = 900.0,
        f_base_max::Float64 = 0.1,
        density_init::Float64 = 300.0,
        temperature_init::Float64 = 273.0,
    )   

        # Initialize with no initial mass
        mass = zeros(Float64, Ntot)
        mass_w = zeros(Float64, Ntot)
        density = fill(density_init, Ntot)
        temperature = fill(temperature_init, Ntot)
        mass_base = 0.0
        runoff = 0.0
        
        # Consistency check
        @assert mass_split < mass_max
        @assert mass_min < mass_split

        # Make sure mass_split is more than 50% of mass_max, so that when
        # surface layer splits, the surface contains less mass than the subsurface layer
        @assert mass_split / mass_max >= 0.5
        
        new(c, Ntot, N, mass_max, mass_split, mass_min, rho_max, f_base_max,
            mass, mass_w, density, temperature, mass_base, runoff)
    end
end


"""
    step!(column::SnowpackColumn, mdot::Float64, dt::Float64) -> Float64

Advance the snowpack column by one time step.

# Arguments
- `column`: The snowpack column to update
- `T2m` : Near-surface air temperature [K]
- `P`: Precipitation rate at surface [kg/m²/s]
- `dt`: Time step [d]
- `f_s`: Fraction of precipitation that is snow [1], default nothing, calculate internally

# Process
1. Apply surface mass flux
2. Handle layer splitting/merging
3. Propagate melt through layers if negative
4. Check for ice formation at base
"""
function step!(column::SnowpackColumn, T2m::Float64, P::Float64, dt::Float64; f_s=nothing)

    if isnothing(f_s)
        # Determine fraction of snow and rain as a function of T2m
        # following Born et al. (2019)
        if T2m > column.c.T0
            f_s = 0.0
        else
            f_s = 1.0
        end
    end

    # Get separate contributions of rain and snow depending on arguments
    P_rain = P * (1.0-f_s)
    P_snow = P - P_rain

    # Convert timestep to seconds internally
    dt_sec = dt * column.c.seconds_per_day

    # Handle accumulation first
    apply_accumulation!(column, P_snow, P_rain, dt_sec)

    # Caculate energy balance
    # to do

    # Handle melt
    #apply_melt!(column, -dmass)
    
    return
end


"""
    apply_accumulation!(column::SnowpackColumn, P_snow::Float64, P_rain::Float64, dt::Float64) -> Float64

Add mass to the surface layer and handle layer dynamics.

"""
function apply_accumulation!(column::SnowpackColumn, P_snow::Float64, P_rain::Float64, dt::Float64)
    
    # If no active layers, create the first one
    if column.N == 0
        column.N = 1
    end
    
    # Add mass to surface layer (first layer)
    column.mass[1] += P_snow * dt
    column.mass_w[1] += P_rain * dt

    # If all layers are full and surface exceeds mass_max, first merge bottom layers
    if column.mass[1] > column.mass_max && column.N == column.Ntot
        merge_bottom_layer!(column)
    end

    # Check if surface layer needs splitting or merging
    if column.mass[1] > column.mass_max
        split_surface_layer!(column)
    elseif column.mass[1] < column.mass_min
        merge_surface_layer!(column)
    end
    
    # Check if mass should be removed from basal layer due to saturation
    if column.N == column.Ntot && column.mass[column.N] > (2*column.mass_max)
        f = column.f_base_max
        mass_to_base = f * column.mass[column.N]
        column.mass[column.N] -= mass_to_base
        column.mass_base += mass_to_base
    end

    return
end

"""
    split_surface_layer!(column::SnowpackColumn)

Split the surface layer when it exceeds mass_max.
Lower part contains mass_split, upper part contains remainder.
Shifts all layers down by one position.

Returns mass flux at base if bottom layers need to be merged.
"""
function split_surface_layer!(column::SnowpackColumn)
    
    # Only split the suface layer if the mass is too high
    @assert column.mass[1] > column.mass_max
    
    # This routine assumes that the bottom layer is currently empty (bottom-layer merging already applied)
    @assert column.N < column.Ntot

    surface_mass = column.mass[1]
    surface_mass_w = column.mass_w[1]
    surface_density = column.density[1]
    surface_temperature = column.temperature[1]

    # Activate another layer
    column.N += 1

    # Shift all layers down
    for i in column.N:-1:3
        column.mass[i] = column.mass[i-1]
        column.mass_w[i] = column.mass_w[i-1]
        column.density[i] = column.density[i-1]
        column.temperature[i] = column.temperature[i-1]
    end
    
    # Subsurface layer (now at index 2) gets mass_split
    column.mass[2] = column.mass_split
    column.mass[1] = surface_mass - column.mass_split
    
    # Then split water proportionally
    water_fraction = column.mass_split / surface_mass
    column.mass_w[2] = surface_mass_w * water_fraction
    column.mass_w[1] = surface_mass_w * (1 - water_fraction)

    # Density is equal in both layers
    column.density[1] = surface_density
    column.density[2] = surface_density
    
    # Temperature is equal in both layers
    column.temperature[1] = surface_temperature
    column.temperature[2] = surface_temperature
    
    return
end


"""
    merge_surface_layer!(column::SnowpackColumn)

Merge surface layer with the second layer when surface mass < mass_min.
If combined mass > 2*mass_split, only partial transfer to keep surface ≤ mass_split.
"""
function merge_surface_layer!(column::SnowpackColumn)

    if column.N == 1 && column.mass[1] < 1e-10
        # Only one empty layer, disable it
        column.N = 0
        reset_column_at_index!(column,1)
        return
    elseif column.N == 1
        # Only one layer, but nothing to merge with, skip merging
        return
    end
    
    surface_mass = column.mass[1]
    subsurface_mass = column.mass[2]
    combined_mass = surface_mass + subsurface_mass
    
    if combined_mass > 2 * column.mass_split
        total = column.mass_split
        transferred_to_surface = column.mass_split - surface_mass
        transferred_water = transferred_to_surface / column.mass[2] * column.mass_w[2]
        
        # Partial transfer: keep surface at mass_split
        column.mass[1] = column.mass_split
        column.mass[2] = combined_mass - column.mass_split
        
        # Proportional transfer of water to the surface
        column.mass_w[1] += transferred_water
        column.mass_w[2] -= transferred_water

        # Weighted average density for surface layer
        column.density[1] = (surface_mass * column.density[1] + 
                            transferred_to_surface * column.density[2]) / total
        
        # Weighted average temperature for surface layer
        column.temperature[1] = (surface_mass * column.temperature[1] + 
                            transferred_to_surface * column.temperature[2]) / total
        
    else
        # Full merge into surface layer
        column.mass[1] = combined_mass
        column.mass_w[1] += column.mass_w[2]

        # Mass-weighted average density
        column.density[1] = (surface_mass * column.density[1] + 
                            subsurface_mass * column.density[2]) / combined_mass
        
        # Mass-weighted average temperature
        column.temperature[1] = (surface_mass * column.temperature[1] + 
                            subsurface_mass * column.temperature[2]) / combined_mass
        
        # Reduce number of active layers by 1
        column.N -= 1

        # Shift all layers up
        for i in 2:(column.N)
            column.mass[i] = column.mass[i+1]
            column.mass_w[i] = column.mass_w[i+1]
            column.density[i] = column.density[i+1]
            column.temperature[i] = column.temperature[i+1]
        end
        
        # Clear the lower, now inactive layer
        reset_column_at_index!(column,column.N+1)
    end

    return
end


"""
    merge_bottom_layer!(column::SnowpackColumn) -> Float64

Merge the two lowest layers to make room for a new surface layer.
This happens when all Ntot are full and surface needs to split.

Returns mass flux at base if merged layer exceeds ice density.
"""
function merge_bottom_layer!(column::SnowpackColumn)

    # Only possible to merge bottom layers if all layers are active
    @assert column.N == column.Ntot
    
    # Deactivate lowest active layer
    column.N -= 1

    # Merge layers N and N+1
    N = column.N
    Np1 = column.N + 1
    
    # Get total masses of combined layers
    combined_mass = column.mass[N] + column.mass[Np1]
    combined_mass_w = column.mass_w[N] + column.mass_w[Np1]

    # Get mass-weighted average density
    combined_density = (column.mass[N] * column.density[N] + 
                       column.mass[Np1] * column.density[Np1]) / 
                       combined_mass
    
    # Get mass-weighted average temperature
    combined_temperature = (column.mass[N] * column.temperature[N] + 
                       column.mass[Np1] * column.temperature[Np1]) / 
                       combined_mass
    
    if combined_density > column.c.rho_i
        # Ice forms - excess mass exits at base

        # Calculate mass of bottom layer at the density of ice (maximum allowed)
        # If layer is too dense, then the excess mass will be sent to the mass_base buffer
        mass_limited_to_ice_density = combined_mass * (column.c.rho_i / combined_density)

        column.mass_base += combined_mass - mass_limited_to_ice_density
        column.mass[N] = mass_limited_to_ice_density
        column.mass_w[N] = combined_mass_w
        column.density[N] = column.c.rho_i
        column.temperature[N] = combined_temperature
    else
        column.mass[N] = combined_mass
        column.mass_w[N] = combined_mass_w
        column.density[N] = combined_density
        column.temperature[N] = combined_temperature
    end
    
    # Reset lower, now unactive, layer to zero values
    reset_column_at_index!(column,Np1)

    return
end

function reset_column_at_index!(column,i::Int)

    column.mass[i] = 0.0
    column.mass_w[i] = 0.0
    column.density[i] = 0.0
    column.temperature[i] = column.c.T0

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
