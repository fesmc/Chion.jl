"""
State accessors and formatted state output.
"""

@inline function _active_column_profile(values, n::Int, idx::Int)
    n == 0 && return Float64[]
    return [@inbounds _get_layer(values, layer_index, idx) for layer_index in 1:n]
end

"""
    _state_dict(N_storage, mass, mass_w, density, temperature, smb_ice, albedo_dynamic, idx, c)

Build a dictionary snapshot for column `idx`. The result includes active-layer
profiles, bulk totals, and surface diagnostics, and allocates new Julia arrays
for the returned profile data.
"""
function _state_dict(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    smb_ice,
    albedo,
    idx::Int,
    c::SnowpackPhysicalConstants,
)
    n = _n_active(N_storage, idx)
    active_solid_mass = _active_column_profile(mass, n, idx)
    active_liquid_water_mass = _active_column_profile(mass_w, n, idx)
    active_density = _active_column_profile(density, n, idx)
    thickness = active_solid_mass ./ active_density
    total_mass = sum(active_solid_mass)
    total_liquid_water = sum(active_liquid_water_mass)

    return Dict(
        "N" => n,
        "mass" => active_solid_mass,
        "mass_w" => active_liquid_water_mass,
        "density" => active_density,
        "total_mass" => total_mass,
        "total_liquid_water" => total_liquid_water,
        "total_wet_mass" => total_mass + total_liquid_water,
        "thickness" => thickness,
        "total_thickness" => sum(thickness),
        "surface_temperature" => n == 0 ? c.T0 : _get_layer(temperature, 1, idx),
        "albedo" => _get_scalar(albedo, idx),
        "smb_ice" => _get_scalar(smb_ice, idx),
    )
end

"""
    get_state(state, idx=1)

Return a dictionary snapshot for column `idx` of `state`. Values are copied
into plain Julia containers so callers can inspect state without mutating the
state.
"""
function get_state(state, idx::Int=1)
    return _state_dict(
        state.N,
        state.mass,
        state.mass_w,
        state.density,
        state.temperature,
        state.smb_ice,
        state.albedo,
        idx,
        state.c,
    )
end

"""
    print_state(state, idx=1)

Print a short formatted summary of column `idx` to standard output. This is a
diagnostic convenience wrapper around `get_state`.
"""
function print_state(state, idx::Int=1)
    snapshot = get_state(state, idx)
    println("=" ^ 60)
    println("Snowpack Column State")
    println("=" ^ 60)
    println("Column index: ", idx)
    println("Active layers: ", snapshot["N"])
    println("Total mass: ", round(snapshot["total_mass"], digits=2), " kg/m^2")
    println("Total thickness: ", round(snapshot["total_thickness"], digits=3), " m")
    println("Surface albedo: ", round(snapshot["albedo"], digits=3))
    println()
end

"""
    compute_auxiliary!(state, idx)

Recompute derived diagnostics for column `idx` in-place. Currently this
updates surface albedo.
"""
function compute_auxiliary!(state, idx::Int)
    update_surface_albedo!(state, idx)
    return nothing
end

"""
    compute_auxiliary!(state)

Recompute derived diagnostics for every column in `state`. Mutates the
state’s auxiliary fields in-place and returns `nothing`.
"""
function compute_auxiliary!(state)
    for idx in 1:ncols(state)
        compute_auxiliary!(state, idx)
    end
    return nothing
end

"""
Domain-wide summary helpers.
"""

const _DOMAIN_SUMMARY_FIELDS =
    (:thickness, :wet_mass, :bulk_density, :base_mass, :smb_ice, :liquid_water, :runoff, :melt, :refreezing, :vapor_mass, :sublimation, :latent_heat_flux_sum, :albedo)

@inline function _column_summary(N, mass, mass_w, density, idx, sample)
    n = N[idx]
    thickness_local = zero(eltype(sample))
    wet_mass_local = zero(eltype(sample))
    solid_mass_local = zero(eltype(sample))
    liquid_water_local = zero(eltype(sample))
    for layer_index in 1:n
        solid = mass[layer_index, idx]
        liquid = mass_w[layer_index, idx]
        rho = density[layer_index, idx]
        solid_mass_local += solid
        wet_mass_local += solid + liquid
        liquid_water_local += liquid
        if solid > zero(solid) && rho > EPS_TINY
            thickness_local += solid / rho
        end
    end
    bulk_density_local =
        thickness_local > EPS_TINY ? solid_mass_local / thickness_local : zero(eltype(sample))
    return thickness_local, wet_mass_local, bulk_density_local, liquid_water_local
end

@inline function _summary_buffers(state, names::NTuple{N,Symbol}) where {N}
    NF = number_type(state.c)
    ncol = ncols(state)
    return NamedTuple{names}(ntuple(_ -> similar(state.mass, NF, ncol), N))
end

@inline function _launch_summary_kernel!(kernel, state, args...)
    kernel! = kernel(_ka_backend(state.mass))
    event = kernel!(args...; ndrange=ncols(state))
    _wait_kernel(event)
    return nothing
end

"""
    _summarize_domain_state_kernel!(...)

KernelAbstractions kernel that summarizes each column into thickness, wet
mass, bulk density, basal mass, ice SMB, liquid water, and runoff arrays. All
output arrays are mutated in-place.
"""
@kernel function _summarize_domain_state_kernel!(
    output,
    state::BESSIState,
)
    idx = @index(Global)
    if idx <= length(state.N)
        thickness_local, wet_mass_local, bulk_density_local, liquid_water_local =
            _column_summary(state.N, state.mass, state.mass_w, state.density, idx, output.thickness)
        output.thickness[idx] = thickness_local
        output.wet_mass[idx] = wet_mass_local
        output.bulk_density[idx] = bulk_density_local
        output.base_mass[idx] = state.mass_base[idx]
        output.smb_ice[idx] = state.smb_ice[idx]
        output.liquid_water[idx] = liquid_water_local
        output.runoff[idx] = state.runoff[idx]
        output.melt[idx] = state.melt[idx]
        output.refreezing[idx] = state.refreezing[idx]
        output.vapor_mass[idx] = state.vapor_mass[idx]
        output.sublimation[idx] = state.sublimation[idx]
        output.latent_heat_flux_sum[idx] = state.latent_heat_flux_sum[idx]
        output.albedo[idx] = state.albedo[idx]
    end
end

"""
    _summarize_year_state_kernel!(...)

KernelAbstractions kernel that summarizes the subset of column diagnostics
needed for equilibrium-year tracking. Mutates the supplied output arrays
in-place.
"""
@kernel function _summarize_year_state_kernel!(
    output,
    state::BESSIState,
)
    idx = @index(Global)
    if idx <= length(state.N)
        thickness_local, wet_mass_local, bulk_density_local, _ =
            _column_summary(state.N, state.mass, state.mass_w, state.density, idx, output.thickness)
        output.thickness[idx] = thickness_local
        output.wet_mass[idx] = wet_mass_local
        output.bulk_density[idx] = bulk_density_local
        output.base_mass[idx] = state.mass_base[idx]
    end
end

"""
    summarize_domain_state!(thickness, wet_mass, bulk_density, base_mass, smb_ice, liquid_water, runoff, melt, refreezing, vapor_mass, sublimation, latent_heat_flux_sum, albedo, state)

Fill preallocated summary arrays with one-column diagnostics from `state`.
Outputs are column-wise totals or aggregates in SI-like model units.
"""
function summarize_domain_state!(
    thickness::AbstractVector,
    wet_mass::AbstractVector,
    bulk_density::AbstractVector,
    base_mass::AbstractVector,
    smb_ice::AbstractVector,
    liquid_water::AbstractVector,
    runoff::AbstractVector,
    melt::AbstractVector,
    refreezing::AbstractVector,
    vapor_mass::AbstractVector,
    sublimation::AbstractVector,
    latent_heat_flux_sum::AbstractVector,
    albedo::AbstractVector,
    state,
)
    output = (; thickness, wet_mass, bulk_density, base_mass, smb_ice, liquid_water,
              runoff, melt, refreezing, vapor_mass, sublimation, latent_heat_flux_sum,
              albedo)
    return _launch_summary_kernel!(
        _summarize_domain_state_kernel!,
        state,
        output,
        state,
    )
end

"""
    summarize_domain_state(state)

Allocate and return a named tuple of per-column summary arrays for `state`.
This is a convenience wrapper around `summarize_domain_state!`.
"""
function summarize_domain_state(state)
    summary = _summary_buffers(state, _DOMAIN_SUMMARY_FIELDS)
    summarize_domain_state!(summary..., state)
    return summary
end

"""
    summarize_year_state!(thickness, wet_mass, bulk_density, base_mass, state)

Fill preallocated arrays with the smaller summary set used to compare
equilibrium years. Mutates the output arrays and returns `nothing`.
"""
function summarize_year_state!(
    thickness::AbstractVector,
    wet_mass::AbstractVector,
    bulk_density::AbstractVector,
        base_mass::AbstractVector,
        state,
)
    output = (; thickness, wet_mass, bulk_density, base_mass)
    return _launch_summary_kernel!(
        _summarize_year_state_kernel!,
        state,
        output,
        state,
    )
end
