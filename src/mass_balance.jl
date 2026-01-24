module MassBalance

export compute_mass_balance, apply_melt

using ..SnowLayers: Snowpack, SnowLayer

# Compute snow accumulation and melt per time step
function compute_mass_balance(snowfall::Float64, melt_energy::Float64)
    latent_heat_fusion = 334000  # J/kg

    melt_mass = melt_energy > 0 ? melt_energy / latent_heat_fusion : 0.0
    net_snow = snowfall - melt_mass

    return net_snow, melt_mass
end

# Apply melt to the snowpack layers
function apply_melt(pack::Snowpack, melt_mass::Float64)
    remaining_melt = melt_mass
    
    # Melt active layers from top to bottom
    for layer in pack.layers
        if layer.active
            layer_mass = layer.thickness * layer.density
            
            if remaining_melt >= layer_mass
                # Entire layer melts
                remaining_melt -= layer_mass
                deactivate_layer!(layer)
            else
                # Partial melting
                melted_thickness = remaining_melt / layer.density
                layer.thickness -= melted_thickness
                remaining_melt = 0.0
            end
            
            if remaining_melt <= 0.0
                break
            end
        end
    end
end

# Deactivate a fully melted layer
function deactivate_layer!(layer::SnowLayer)
    layer.thickness = 0.0
    layer.density = 0.0
    layer.temperature = 0.0
    layer.liquid_water = 0.0
    layer.active = false
end

end
