cd(@__DIR__)
import Pkg; Pkg.activate("$(homedir())/.JuliaEnvironments/myanalysis")
using Revise

include("src/snow_layers.jl")
include("src/mass_balance.jl")
include("src/energy_balance.jl")

using .SnowLayers
using .MassBalance
using .EnergyBalance
using CairoMakie

# Initialize snowpack with static max layers
N = 5
snowpack = SnowLayers.Snowpack(N)

# Example meteorological inputs over 24 hours
hours = 0:1:24
snowfall_series = [0.01 * rand() for _ in hours]  # Random snowfall in m
shortwave_in_series = [200.0 + 50.0 * sin(π * h / 12) for h in hours]  # Diurnal shortwave
net_energy_series = Float64[]

for (t, snowfall) in enumerate(snowfall_series)
    # Energy inputs
    shortwave_in = shortwave_in_series[t]
    longwave_in = 250.0
    longwave_out = 300.0
    sensible_heat = 30.0
    latent_heat = -10.0
    ground_flux = 5.0

    # Compute energy and mass balance
    surface_energy = EnergyBalance.compute_energy_balance(shortwave_in, longwave_in, longwave_out,
                                                          sensible_heat, latent_heat, ground_flux)
    push!(net_energy_series, surface_energy)

    SnowLayers.add_snowfall(snowpack, snowfall, -5.0)

    # Apply melt
    _, melt_mass = MassBalance.compute_mass_balance(snowfall, max(surface_energy, 0) * 3600)
    MassBalance.apply_melt(snowpack, melt_mass)
    
    # Distribute energy and update temperatures
    energy_inputs = EnergyBalance.distribute_energy(snowpack, surface_energy * 3600)
    SnowLayers.update_temperature(snowpack, energy_inputs)
end

# Plotting snowpack depth over time
total_thickness = [sum(l.thickness for l in snowpack.layers if l.active) for _ in hours]

begin
    f = Figure(size = (800, 400))
    ax1 = Axis(f[1, 1], xlabel = "Time (hours)", ylabel = "Shortwave radiation (W/m^2)")
    lines!(ax1, hours, shortwave_in_series, color = :blue)
    ax2 = Axis(f[2, 1], xlabel = "Time (hours)", ylabel = "Snowfall (m)")
    lines!(ax2, hours, snowfall_series, color = :blue)
    #ax3 = Axis(f[3, 1], xlabel = "Time (hours)", ylabel = "Snow Depth (m)")
    #lines!(ax3, hours, total_thickness, color = :blue)

    save("plots/snowpack_evolution_static_layers.png", f)
    f
end
