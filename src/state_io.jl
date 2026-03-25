"""
State accessors and formatted state output.
"""

"""
    get_state(column::SnowpackColumn) -> Dict

Get the current state of the snowpack column.

Returns a dictionary with:
- `N`: Number of active layers
- `mass` / `solid_mass`: Solid snow/ice mass in each active layer [kg/m^2]
- `mass_w` / `liquid_water_mass`: Liquid water mass in each active layer [kg/m^2]
- `density`: Density in each active layer [kg/m^3]
- `total_mass`: Total mass in column [kg/m^2]
- `total_liquid_water`: Total liquid water mass in column [kg/m^2]
- `total_wet_mass`: Total snow plus liquid water mass in column [kg/m^2]
- `thickness`: Thickness of each active layer [m]
- `total_thickness`: Total column thickness [m]
- `surface_temperature`: Surface temperature [K]
- `snow_cover`: Diagnosed snow cover fraction [1]
- `surface_albedo`: Current surface albedo used for shortwave absorption [1]
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
            "surface_albedo" => column.albedo_dynamic,
            "albedo_dynamic" => column.albedo_dynamic,
            "smb_ice" => column.smb_ice,
            "ice_sheet_smb" => column.smb_ice,
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
        "surface_albedo" => column.albedo_dynamic,
        "albedo_dynamic" => column.albedo_dynamic,
        "smb_ice" => column.smb_ice,
        "ice_sheet_smb" => column.smb_ice,
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
    println("Total mass: ", round(state["total_mass"], digits=2), " kg/m^2")
    println("Total thickness: ", round(state["total_thickness"], digits=3), " m")
    println("Snow cover: ", round(state["snow_cover"], digits=3))
    println("Surface albedo: ", round(state["surface_albedo"], digits=3))
    println()

    if state["N"] > 0
        println("Layer details (surface = 1):")
        println("-" ^ 60)
        println("Layer | Mass (kg/m^2) | Density (kg/m^3) | Thickness (m)")
        println("-" ^ 60)
        for layer_index in 1:state["N"]
            @printf(
                "%5d | %12.2f | %15.1f | %13.4f\n",
                layer_index,
                state["mass"][layer_index],
                state["density"][layer_index],
                state["thickness"][layer_index],
            )
        end
        println("=" ^ 60)
    end
end
