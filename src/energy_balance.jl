module EnergyBalance

export compute_energy_balance, distribute_energy

using ..SnowLayers: Snowpack, SnowLayer

# Compute net surface energy flux
function compute_energy_balance(shortwave_in::Float64, longwave_in::Float64, 
                                longwave_out::Float64, sensible_heat::Float64,
                                latent_heat::Float64, ground_flux::Float64)

    net_radiation = shortwave_in + longwave_in - longwave_out
    net_flux = net_radiation + sensible_heat + latent_heat + ground_flux

    return net_flux  # Positive = warming, Negative = cooling
end

# Distribute surface energy to active snowpack layers (simplified for conduction)
function distribute_energy(pack::Snowpack, surface_energy::Float64)
    energy_inputs = zeros(Float64, length(pack.layers))
    
    for (i, layer) in enumerate(pack.layers)
        if layer.active
            energy_inputs[i] = surface_energy * exp(-0.5 * (i - 1))  # Exponential decay into snow
        end
    end
    
    return energy_inputs
end

end
