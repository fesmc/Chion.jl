#=
Column diagnostics of a BESSI state and a per-column snapshot for inspection.
=#

Base.@propagate_inbounds function _column_summary(N, mass, mass_w, density, idx, sample)
    n = @inbounds N[idx]
    thickness_local = zero(eltype(sample))
    wet_mass_local = zero(eltype(sample))
    solid_mass_local = zero(eltype(sample))
    liquid_water_local = zero(eltype(sample))
    @inbounds begin
        for layer_index in 1:n
            solid = mass[layer_index, idx]
            liquid = mass_w[layer_index, idx]
            rho = density[layer_index, idx]
            solid_mass_local += solid
            wet_mass_local += solid + liquid
            liquid_water_local += liquid
            valid_layer = (solid > zero(solid)) & (rho > EPS_TINY)
            safe_density = ifelse(valid_layer, rho, one(rho))
            thickness_local += ifelse(valid_layer, solid / safe_density, zero(solid))
        end
    end
    valid_thickness = thickness_local > EPS_TINY
    safe_thickness = ifelse(valid_thickness, thickness_local, one(thickness_local))
    bulk_density_local = ifelse(
        valid_thickness,
        solid_mass_local / safe_thickness,
        zero(eltype(sample)),
    )
    return thickness_local, wet_mass_local, bulk_density_local, liquid_water_local
end

"""Thickness, wet mass, bulk density and liquid water of every column (see `update_diagnostics!`)."""
@kernel function _summarize_column_kernel!(output, fields)
    idx = @index(Global)
    @inbounds begin
        thickness, wet_mass, bulk_density, liquid_water =
            _column_summary(fields.N, fields.mass, fields.mass_w, fields.density, idx, output.thickness)
        output.thickness[idx] = thickness
        output.wet_mass[idx] = wet_mass
        output.bulk_density[idx] = bulk_density
        output.liquid_water[idx] = liquid_water
    end
end

"""
    column_profile(model::BESSIModel, idx = 1) -> NamedTuple

The active layers of column `idx` (host copies): `N`, `mass`, `mass_w`, `density`, `temperature`,
`thickness` (m), and the column's `albedo`, `Tsrf` and `smb_ice`.
"""
function column_profile(model::BESSIModel, idx::Int = 1)
    state = cpu_state(model)
    n = state.N[idx]
    mass = state.mass[1:n, idx]
    density = state.density[1:n, idx]
    return (N = n, mass = mass, mass_w = state.mass_w[1:n, idx], density = density,
        temperature = state.temperature[1:n, idx], thickness = mass ./ density,
        albedo = state.albedo[idx], Tsrf = state.Tsrf[idx], smb_ice = state.smb_ice[idx])
end
