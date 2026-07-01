"""
State access helpers shared across column, domain, and array-backed kernels.
"""

"""
    _resolve_keyword_alias(preferred_value, legacy_value, preferred_name, legacy_name)

Resolve a preferred keyword value against a legacy alias. Returns the preferred
value when present, the legacy value otherwise, and throws if both are
provided with different values.
"""
@inline function _resolve_keyword_alias(
    preferred_value,
    legacy_value,
    preferred_name::AbstractString,
    legacy_name::AbstractString,
)
    if !isnothing(preferred_value) && !isnothing(legacy_value) &&
       !isequal(preferred_value, legacy_value)
        error(
            "Received both `$preferred_name` and `$legacy_name` with different values. " *
            "Use one or provide matching values.",
        )
    end
    return isnothing(preferred_value) ? legacy_value : preferred_value
end

@inline _n_active(N::Base.RefValue{<:Integer}, ::Int) = Int(N[])
@inline _n_active(N::AbstractVector{<:Integer}, idx::Int) = @inbounds Int(N[idx])

@inline _set_n_active!(N::Base.RefValue{<:Integer}, ::Int, value::Int) = (N[] = value)
@inline _set_n_active!(N::AbstractVector{<:Integer}, idx::Int, value::Int) = (@inbounds N[idx] = value)

@inline _get_scalar(x::Base.RefValue, ::Int) = x[]
@inline _get_scalar(x::AbstractVector, idx::Int) = @inbounds x[idx]

@inline _set_scalar!(x::Base.RefValue, ::Int, value) = (x[] = value)
@inline _set_scalar!(x::AbstractVector, idx::Int, value) = (@inbounds x[idx] = value)

@inline _get_layer(x::AbstractVector, layer_index::Int, ::Int) = @inbounds x[layer_index]
@inline _get_layer(x::AbstractMatrix, layer_index::Int, idx::Int) = @inbounds x[layer_index, idx]

@inline _set_layer!(x::AbstractVector, layer_index::Int, ::Int, value) = (@inbounds x[layer_index] = value)
@inline _set_layer!(x::AbstractMatrix, layer_index::Int, idx::Int, value) = (@inbounds x[layer_index, idx] = value)

"""
    _bulk_snow_density(N_storage, mass, density, idx)

Compute the bulk density of active snow in column `idx` from solid mass and
layer thickness. Returns zero when the column has no valid snow mass.
"""
@inline function _bulk_snow_density(
    N_storage,
    mass,
    density,
    idx::Int,
)
    n = _n_active(N_storage, idx)
    if n <= 0
        return zero(eltype(density))
    end

    total_mass = zero(eltype(density))
    total_thickness = zero(eltype(density))
    @inbounds for layer_index in 1:n
        layer_mass = _get_layer(mass, layer_index, idx)
        layer_density = _get_layer(density, layer_index, idx)
        if layer_mass > zero(layer_mass) && layer_density > EPS_TINY
            total_mass += layer_mass
            total_thickness += layer_mass / layer_density
        end
    end

    if total_mass <= zero(total_mass) || total_thickness <= EPS_TINY
        return zero(total_mass)
    end
    return total_mass / total_thickness
end

"""
    _total_snow_water_mass(N_storage, mass, mass_w, idx)

Return the total wet mass in column `idx`, including solid snow and liquid
water. Negative layer masses are clipped to zero in the accumulation.
"""
@inline function _total_snow_water_mass(
    N_storage,
    mass,
    mass_w,
    idx::Int,
)
    n = _n_active(N_storage, idx)
    if n <= 0
        return zero(eltype(mass))
    end

    total_wet_mass = zero(eltype(mass))
    @inbounds for layer_index in 1:n
        total_wet_mass += max(_get_layer(mass, layer_index, idx), zero(eltype(mass))) +
                          max(_get_layer(mass_w, layer_index, idx), zero(eltype(mass)))
    end
    return total_wet_mass
end

"""
    _column_has_liquid_water(N_storage, mass_w, idx)

Return `true` when any active layer in column `idx` contains liquid water
above the empty-layer tolerance.
"""
@inline function _column_has_liquid_water(
    N_storage,
    mass_w,
    idx::Int,
)
    n = _n_active(N_storage, idx)
    @inbounds for layer_index in 1:n
        if _get_layer(mass_w, layer_index, idx) > EPS_TINY
            return true
        end
    end
    return false
end

@inline _column_has_liquid_water(domain, idx::Int) =
    _column_has_liquid_water(domain.N, domain.mass_w, idx)
