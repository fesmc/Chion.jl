"""
Liquid-water percolation following the current BESSI overflow routine.
"""

"""
    go_percolation!(
        solid_mass::Vector{Float64},
        liquid_water_mass::Vector{Float64},
        snow_density::Vector{Float64},
        ice_density::Float64,
        water_density::Float64;
        max_lwc::Float64 = 0.1,
        rho_i_tol::Float64 = 10.0,
    ) -> Float64

Translate the active Fortran `percolate` logic to Julia.

- `solid_mass`: layer snow/ice mass [kg m^-2]
- `liquid_water_mass`: layer liquid-water mass [kg m^-2]
- `snow_density`: layer snow density [kg m^-3]
- `rho_i_tol`: accepted for backward compatibility but unused by the current
  BESSI overflow routine
- returns `runoff` produced by percolation in this call [kg m^-2]
"""
@inline function _is_lowest_active_snow_layer(
    solid_mass::AbstractVector{Float64},
    layer_index::Int,
)
    return layer_index == length(solid_mass) || solid_mass[layer_index + 1] <= 0.0
end

function go_percolation!(
    solid_mass::AbstractVector{Float64},
    liquid_water_mass::AbstractVector{Float64},
    snow_density::AbstractVector{Float64},
    ice_density::Float64,
    water_density::Float64;
    max_lwc::Float64 = 0.1,
    rho_i_tol::Float64 = 10.0,
)
    n_layers = length(solid_mass)
    @assert length(liquid_water_mass) == n_layers
    @assert length(snow_density) == n_layers

    runoff = 0.0
    for layer_index in 1:n_layers
        if solid_mass[layer_index] <= 0.0
            runoff += liquid_water_mass[layer_index]
            liquid_water_mass[layer_index] = 0.0
            continue
        end

        pore_volume = solid_mass[layer_index] / snow_density[layer_index] -
                      solid_mass[layer_index] / ice_density
        if pore_volume <= EPS_TINY
            excess_water = liquid_water_mass[layer_index]
            liquid_water_mass[layer_index] = 0.0
            if _is_lowest_active_snow_layer(solid_mass, layer_index)
                runoff += excess_water
            else
                liquid_water_mass[layer_index + 1] += excess_water
            end
            continue
        end

        liquid_water_content = liquid_water_mass[layer_index] / water_density / pore_volume
        if liquid_water_content > max_lwc
            excess_water = (liquid_water_content - max_lwc) * pore_volume * water_density
            liquid_water_mass[layer_index] -= excess_water

            if _is_lowest_active_snow_layer(solid_mass, layer_index)
                runoff += excess_water
            else
                liquid_water_mass[layer_index + 1] += excess_water
            end
        end
    end

    return runoff
end

"""
    go_percolation!(column::SnowpackColumn; max_lwc=0.1, rho_i_tol=10.0) -> Float64

Apply percolation to active layers of a `SnowpackColumn`.
Returned runoff is also added to `column.runoff`.
"""
function go_percolation!(
    column::SnowpackColumn;
    max_lwc::Float64 = 0.1,
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
