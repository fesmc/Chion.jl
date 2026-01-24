module SnowLayers

export SnowLayer, Snowpack, add_snowfall, update_temperature

mutable struct SnowLayer
    thickness::Float64          # Layer thickness (m)
    density::Float64            # Layer density (kg/m³)
    temperature::Float64        # Layer temperature (°C)
    liquid_water::Float64       # Liquid water content (m³)
    active::Bool                # Flag to indicate if the layer is active
end

mutable struct Snowpack
    layers::Vector{SnowLayer}
    max_layers::Int
end

# Initialize a snowpack with inactive layers
function Snowpack(max_layers::Int)
    layers = [SnowLayer(0.0, 0.0, 0.0, 0.0, false) for _ in 1:max_layers]
    return Snowpack(layers, max_layers)
end

# Add snowfall to the snowpack, activating an inactive layer if available
function add_snowfall(pack::Snowpack, snowfall::Float64, snow_temp::Float64)
    if snowfall <= 0.0
        return
    end
    
    # Find the first inactive layer (or overwrite the oldest active layer if all are active)
    idx = findfirst(l -> !l.active, pack.layers)
    
    if isnothing(idx)
        # All layers are active, shift layers down to make space
        shift_layers_down(pack)
        idx = 1
    end
    
    # Add new snow to the selected layer
    pack.layers[idx] = SnowLayer(snowfall, 100.0, snow_temp, 0.0, true)
end

# Shift layers down to discard the oldest and make room for new snowfall
function shift_layers_down(pack::Snowpack)
    for i in reverse(2:pack.max_layers)
        pack.layers[i] = pack.layers[i - 1]
    end
    # Deactivate the topmost layer
    pack.layers[1] = SnowLayer(0.0, 0.0, 0.0, 0.0, false)
end

# Update temperature of each active layer based on energy input
function update_temperature(pack::Snowpack, energy_inputs::Vector{Float64})
    specific_heat_snow = 2.1e3  # J/(kg·K)
    
    for (i, layer) in enumerate(pack.layers)
        if layer.active
            energy_input = energy_inputs[i]
            mass = layer.density * layer.thickness
            delta_temp = energy_input / (mass * specific_heat_snow)
            layer.temperature += delta_temp
        end
    end
end

end
