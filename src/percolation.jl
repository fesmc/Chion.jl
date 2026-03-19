"""
Liquid-water percolation following the original BESSI-style routine.
"""

"""
    go_percolation!(
        solid_mass::Vector{Float64},
        liquid_water_mass::Vector{Float64},
        snow_density::Vector{Float64},
        ice_density::Float64,
        water_density::Float64;
        max_lwc::Float64 = 0.05,
        rho_i_tol::Float64 = 10.0,
    ) -> Float64

Translate the Fortran `go_percolation` logic to Julia.

- `solid_mass`: layer snow/ice mass [kg m^-2]
- `liquid_water_mass`: layer liquid-water mass [kg m^-2]
- `snow_density`: layer snow density [kg m^-3]
- returns `runoff` produced by percolation in this call [kg m^-2]
"""
function go_percolation!(
    solid_mass::AbstractVector{Float64},
    liquid_water_mass::AbstractVector{Float64},
    snow_density::AbstractVector{Float64},
    ice_density::Float64,
    water_density::Float64;
    max_lwc::Float64 = 0.05,
    rho_i_tol::Float64 = 10.0,
)
    n_layers = length(solid_mass)
    @assert length(liquid_water_mass) == n_layers
    @assert length(snow_density) == n_layers

    runoff = 0.0
    layer_index = 1
    while layer_index <= n_layers
        if solid_mass[layer_index] > 0.0
            if snow_density[layer_index] > ice_density - rho_i_tol
                # Very dense snow: push all liquid water downward.
                percolating_liquid_water = liquid_water_mass[layer_index]
                liquid_water_mass[layer_index] = 0.0
                if layer_index < n_layers
                    if solid_mass[layer_index + 1] > 0.0
                        liquid_water_mass[layer_index + 1] += percolating_liquid_water
                    else
                        runoff += percolating_liquid_water
                    end
                else
                    runoff += percolating_liquid_water
                end
            else
                liquid_water_content = liquid_water_mass[layer_index] /
                                       solid_mass[layer_index] /
                                       water_density /
                                       (1.0 / snow_density[layer_index] - 1.0 / ice_density)
                if liquid_water_content > max_lwc
                    percolating_liquid_water = (
                        liquid_water_content - max_lwc
                    ) * water_density * solid_mass[layer_index] *
                        (1.0 / snow_density[layer_index] - 1.0 / ice_density)
                    liquid_water_mass[layer_index] -= percolating_liquid_water
                    if layer_index < n_layers
                        if solid_mass[layer_index + 1] > 0.0
                            liquid_water_mass[layer_index + 1] += percolating_liquid_water
                        else
                            runoff += percolating_liquid_water
                        end
                    else
                        runoff += percolating_liquid_water
                    end
                end
            end
        else
            # Keep BESSI control flow: stop at the first empty layer.
            break
        end
        layer_index += 1
    end

    return runoff
end

"""
    go_percolation!(column::SnowpackColumn; max_lwc=0.05, rho_i_tol=10.0) -> Float64

Apply percolation to active layers of a `SnowpackColumn`.
Returned runoff is also added to `column.runoff`.
"""
function go_percolation!(
    column::SnowpackColumn;
    max_lwc::Float64 = 0.05,
    rho_i_tol::Float64 = 10.0,
)
    if column.N <= 0 || column.mass[1] <= 0.0
        return 0.0
    end

    @views runoff = go_percolation!(
        column.mass[1:column.N],
        column.mass_w[1:column.N],
        column.density[1:column.N],
        column.c.rho_i,
        column.c.rho_w;
        max_lwc=max_lwc,
        rho_i_tol=rho_i_tol,
    )
    column.runoff += runoff
    return runoff
end
