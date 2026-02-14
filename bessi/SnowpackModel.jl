"""
Column-based snowpack model with dynamic layering.
Based on Born et al. (2019) algorithm.

This initial version focuses on mass conservation and layer dynamics.
"""

module SnowpackModel

using Printf

export SnowpackColumn, step!, get_state
export SnowpackPhysicalConstants

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
- `rho_ice::Float64`: Ice density [kg/m³] (default: 917)

# State variables
- `mass::Vector{Float64}`: Mass of snow+water in each layer [kg/m²]
- `mass_snow::Vector{Float64}`: Mass of snow in each layer [kg/m²]
- `mass_water::Vector{Float64}`: Mass of water in each layer [kg/m²]
- `density::Vector{Float64}`: Density of snow in each layer [kg/m³]

"""

mutable struct SnowpackColumn
    # Constants
    c::SnowpackPhysicalConstants

    # Grid parameters
    Ntot::Int
    N::Int
    kbase::Int

    # Model parameters
    mass_max::Float64           # kg/m²
    mass_split::Float64         # kg/m²
    mass_min::Float64           # kg/m²
    
    #ζmax::Float64   # Maximum liquid water content

    # State variables
    mass::Vector{Float64}           # kg/m²
    mass_snow::Vector{Float64}      # kg/m²
    mass_ice::Vector{Float64}       # kg/m²
    density::Vector{Float64}        # kg/m³
    temperature::Vector{Float64}    # K

    function SnowpackColumn(;
        c::SnowpackPhysicalConstants = SnowpackPhysicalConstants(),
        Ntot::Int = 15,
        N::Int = 1,
        mass_max::Float64 = 500.0,
        mass_split::Float64 = 300.0,
        mass_min::Float64 = 100.0,
        density_init::Float64 = 300.0,
        temperature_init::Float64 = 273.0,
    )   

        # Get index of base layer
        kbase = Ntot - N + 1

        # Initialize with no initial mass
        mass = zeros(Float64, Ntot)
        mass_snow = zeros(Float64, Ntot)
        mass_water = zeros(Float64, Ntot)
        density = fill(density_init, Ntot)
        temperature = fill(temperature_init, Ntot)

        new(c, Ntot, N, kbase, mass_max, mass_split, mass_min, 
            mass, mass_snow, mass_water, density, temperature)
    end
end


"""
    step!(column::SnowpackColumn, mdot::Float64, dt::Float64) -> Float64

Advance the snowpack column by one time step.

# Arguments
- `column`: The snowpack column to update
- `mdot`: Mass rate at surface [kg/m²/s] (positive = accumulation, negative = melt)
- `dt`: Time step [d]

# Returns
- `mdot_base`: Mass flux at the base [kg/m²/s] (positive when ice accumulates)

# Process
1. Apply surface mass flux
2. Handle layer splitting/merging
3. Propagate melt through layers if negative
4. Check for ice formation at base
"""
function step!(column::SnowpackColumn, mdot::Float64, dt::Float64)

    # Convert timestep to seconds internally
    dt_sec = dt * column.seconds_per_day

    dmass = mdot * dt_sec   # Total mass change for this timestep [kg/m²]
    mdot_base = 0.0         # Mass flux out at base
    
    # Handle accumulation (positive mdot)
    if dmass > 0
        mdot_base = apply_accumulation!(column, dmass)
    # Handle melt (negative mdot)
    elseif dmass < 0
        mdot_base = apply_melt!(column, -dmass)
    end
    
    return mdot_base / dt_sec   # Convert back to rate
end


"""
    apply_accumulation!(column::SnowpackColumn, dmass::Float64) -> Float64

Add mass to the surface layer and handle layer dynamics.

Returns mass flux at base if layers need to be merged.
"""
function apply_accumulation!(column::SnowpackColumn, dmass::Float64)
    mdot_base = 0.0
    
    # If no active layers, create the first one
    if column.N == 0
        column.N = 1
    end
    
    # Add mass to surface layer (layer end)
    column.mass[end] += dmass
    
    # Check if surface layer needs splitting
    while column.mass[end] > column.mass_max && column.N < column.Ntot
        mdot_base += split_surface_layer!(column)
    end
    
    # If all layers are full and surface still exceeds mass_max, merge bottom layers
    if column.mass[end] > column.mass_max && column.N == column.Ntot
        mdot_base += merge_bottom_layers!(column)
        # Try splitting again after merging
        if column.mass[end] > column.mass_max
            mdot_base += split_surface_layer!(column)
        end
    end
    
    return mdot_base
end


"""
    apply_melt!(column::SnowpackColumn, dmass::Float64) -> Float64

Remove mass from the column, propagating melt through layers.
Mass is removed from surface but redistributed through column until base.

Returns mass flux at base if ice forms.
"""
function apply_melt!(column::SnowpackColumn, dmass::Float64)
    mdot_base = 0.0
    
    if column.N == 0
        return mdot_base
    end
    
    remaining_melt = dmass
    
    # Remove from surface
    if column.mass[1] >= remaining_melt
        column.mass[1] -= remaining_melt
        remaining_melt = 0.0
    else
        remaining_melt -= column.mass[1]
        column.mass[1] = 0.0
    end
    
    # Propagate remaining melt through layers (mass redistribution)
    # In reality, melt percolates down and redistributes
    # For now, we'll move mass from surface to deeper layers
    if remaining_melt > 0
        # Redistribute the melt mass downward through the column
        layer = 2
        while remaining_melt > 0 && layer <= column.N
            # Add melt mass to this layer (it's moving down)
            # This maintains total column mass
            column.mass[layer] += remaining_melt
            
            # Check if this layer can hold it without exceeding ice density
            thickness = column.mass[layer] / column.density[layer]
            max_mass_at_ice = thickness * column.rho_ice
            
            if column.mass[layer] > max_mass_at_ice
                # Layer has densified to ice - mass must exit at base
                excess = column.mass[layer] - max_mass_at_ice
                column.mass[layer] = max_mass_at_ice
                column.density[layer] = column.rho_ice
                mdot_base += excess
                remaining_melt = 0.0
            else
                remaining_melt = 0.0
            end
            
            layer += 1
        end
        
        # If we've gone through all layers, add to base layer
        if remaining_melt > 0 && column.N > 0
            base_layer = column.N
            column.mass[base_layer] += remaining_melt
            
            # Check for ice formation at base
            thickness = column.mass[base_layer] / column.density[base_layer]
            max_mass_at_ice = thickness * column.rho_ice
            
            if column.mass[base_layer] > max_mass_at_ice
                excess = column.mass[base_layer] - max_mass_at_ice
                column.mass[base_layer] = max_mass_at_ice
                column.density[base_layer] = column.rho_ice
                mdot_base += excess
            end
        end
    end
    
    # Handle layer merging if surface layer is too small
    while column.N > 0 && column.mass[1] < column.mass_min
        merge_surface_layer!(column)
    end
    
    return mdot_base
end


"""
    split_surface_layer!(column::SnowpackColumn) -> Float64

Split the surface layer when it exceeds mass_max.
Lower part contains mass_split, upper part contains remainder.
Shifts all layers down by one position.

Returns mass flux at base if bottom layers need to be merged.
"""
function split_surface_layer!(column::SnowpackColumn)
    mdot_base = 0.0
    
    if column.mass[end] <= column.mass_max
        return mdot_base
    end
    
    surface_mass = column.mass[end]
    surface_density = column.density[end]
    
    # If all layers are full, merge bottom two first
    if column.N == column.Ntot
        mdot_base = merge_bottom_layers!(column)
    end
    
    # Shift all layers down
    for i in column.N:-1:1
        column.mass[i+1] = column.mass[i]
        column.density[i+1] = column.density[i]
    end
    
    # Split: lower layer (now at index 2) gets mass_split
    column.mass[2] = column.mass_split
    column.density[2] = surface_density
    
    # Upper layer (index 1) gets remainder
    column.mass[1] = surface_mass - column.mass_split
    column.density[1] = surface_density
    
    column.N += 1
    
    return mdot_base
end


"""
    merge_surface_layer!(column::SnowpackColumn)

Merge surface layer with the second layer when surface mass < mass_min.
If combined mass > 2*mass_split, only partial transfer to keep surface ≤ mass_split.
"""
function merge_surface_layer!(column::SnowpackColumn)
    if column.N <= 1
        # Only one or no layers, nothing to merge with
        if column.N == 1 && column.mass[1] < 1e-10
            column.N = 0
        end
        return
    end
    
    surface_mass = column.mass[1]
    second_mass = column.mass[2]
    combined_mass = surface_mass + second_mass
    
    if combined_mass > 2 * column.mass_split
        # Partial transfer: keep surface at mass_split
        transfer_mass = column.mass_split - surface_mass
        column.mass[1] = column.mass_split
        column.mass[2] -= transfer_mass
        
        # Weighted average density for surface layer
        total = column.mass_split
        column.density[1] = (surface_mass * column.density[1] + 
                            transfer_mass * column.density[2]) / total
    else
        # Full merge
        # Mass-weighted average density
        column.density[1] = (surface_mass * column.density[1] + 
                            second_mass * column.density[2]) / combined_mass
        column.mass[1] = combined_mass
        
        # Shift all layers up
        for i in 2:(column.N-1)
            column.mass[i] = column.mass[i+1]
            column.density[i] = column.density[i+1]
        end
        
        # Clear the last active layer
        column.mass[column.N] = 0.0
        column.N -= 1
    end
end


"""
    merge_bottom_layers!(column::SnowpackColumn) -> Float64

Merge the two lowest layers to make room for a new surface layer.
This happens when all Ntot are full and surface needs to split.

Returns mass flux at base if merged layer exceeds ice density.
"""
function merge_bottom_layers!(column::SnowpackColumn)
    mdot_base = 0.0
    
    if column.N < 2
        return mdot_base
    end
    
    # Merge layers N-1 and N
    bottom = column.N
    second_bottom = column.N - 1
    
    combined_mass = column.mass[bottom] + column.mass[second_bottom]
    
    # Mass-weighted average density
    combined_density = (column.mass[bottom] * column.density[bottom] + 
                       column.mass[second_bottom] * column.density[second_bottom]) / 
                       combined_mass
    
    # Check if combined layer exceeds ice density
    thickness = combined_mass / combined_density
    max_mass_at_ice = thickness * column.rho_ice
    
    if combined_mass > max_mass_at_ice
        # Ice forms - excess mass exits at base
        mdot_base = combined_mass - max_mass_at_ice
        column.mass[second_bottom] = max_mass_at_ice
        column.density[second_bottom] = column.rho_ice
    else
        column.mass[second_bottom] = combined_mass
        column.density[second_bottom] = combined_density
    end
    
    # Clear the bottom layer
    column.mass[bottom] = 0.0
    column.N -= 1
    
    return mdot_base
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
