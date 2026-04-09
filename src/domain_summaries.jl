"""
Domain-wide summary helpers.
"""

"""
    _summarize_domain_state_kernel!(...)

KernelAbstractions kernel that summarizes each column into thickness, wet
mass, bulk density, basal mass, ice SMB, liquid water, and runoff arrays. All
output arrays are mutated in-place.
"""
@kernel function _summarize_domain_state_kernel!(
    thickness,
    wet_mass,
    bulk_density,
    base_mass,
    smb_ice,
    liquid_water,
    runoff,
    N,
    mass,
    mass_w,
    density,
    mass_base_state,
    smb_ice_state,
    runoff_state,
)
    idx = @index(Global)
    if idx <= length(N)
        n = N[idx]
        thickness_local = zero(eltype(thickness))
        wet_mass_local = zero(eltype(wet_mass))
        solid_mass_local = zero(eltype(wet_mass))
        liquid_water_local = zero(eltype(liquid_water))
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
        thickness[idx] = thickness_local
        wet_mass[idx] = wet_mass_local
        bulk_density[idx] = thickness_local > EPS_TINY ? solid_mass_local / thickness_local : zero(eltype(bulk_density))
        base_mass[idx] = mass_base_state[idx]
        smb_ice[idx] = smb_ice_state[idx]
        liquid_water[idx] = liquid_water_local
        runoff[idx] = runoff_state[idx]
    end
end

"""
    _summarize_cycle_state_kernel!(...)

KernelAbstractions kernel that summarizes the subset of column diagnostics
needed for equilibrium-cycle tracking. Mutates the supplied output arrays
in-place.
"""
@kernel function _summarize_cycle_state_kernel!(
    thickness,
    wet_mass,
    bulk_density,
    base_mass,
    N,
    mass,
    mass_w,
    density,
    mass_base_state,
)
    idx = @index(Global)
    if idx <= length(N)
        n = N[idx]
        thickness_local = zero(eltype(thickness))
        wet_mass_local = zero(eltype(wet_mass))
        solid_mass_local = zero(eltype(wet_mass))
        for layer_index in 1:n
            solid = mass[layer_index, idx]
            liquid = mass_w[layer_index, idx]
            rho = density[layer_index, idx]
            solid_mass_local += solid
            wet_mass_local += solid + liquid
            if solid > zero(solid) && rho > EPS_TINY
                thickness_local += solid / rho
            end
        end
        thickness[idx] = thickness_local
        wet_mass[idx] = wet_mass_local
        bulk_density[idx] = thickness_local > EPS_TINY ? solid_mass_local / thickness_local : zero(eltype(bulk_density))
        base_mass[idx] = mass_base_state[idx]
    end
end

"""
    summarize_domain_state!(thickness, wet_mass, bulk_density, base_mass, smb_ice, liquid_water, runoff, domain; backend=:threads)

Fill preallocated summary arrays with one-column diagnostics from `domain`.
Outputs are column-wise totals or aggregates in SI-like model units, and the
chosen backend controls whether the work runs on Julia threads or a
KernelAbstractions kernel.
"""
function summarize_domain_state!(
    thickness::AbstractVector,
    wet_mass::AbstractVector,
    bulk_density::AbstractVector,
    base_mass::AbstractVector,
    smb_ice::AbstractVector,
    liquid_water::AbstractVector,
    runoff::AbstractVector,
    domain::AbstractSnowpackDomain;
    backend::Symbol=:threads,
)
    if backend == :threads
        Threads.@threads for idx in 1:column_count(domain)
            n = domain.N[idx]
            thickness_local = zero(eltype(thickness))
            wet_mass_local = zero(eltype(wet_mass))
            solid_mass_local = zero(eltype(wet_mass))
            liquid_water_local = zero(eltype(liquid_water))
            @inbounds for layer_index in 1:n
                solid = domain.mass[layer_index, idx]
                liquid = domain.mass_w[layer_index, idx]
                rho = domain.density[layer_index, idx]
                solid_mass_local += solid
                wet_mass_local += solid + liquid
                liquid_water_local += liquid
                if solid > zero(solid) && rho > EPS_TINY
                    thickness_local += solid / rho
                end
            end
            thickness[idx] = thickness_local
            wet_mass[idx] = wet_mass_local
            bulk_density[idx] = thickness_local > EPS_TINY ? solid_mass_local / thickness_local : zero(eltype(bulk_density))
            base_mass[idx] = domain.mass_base[idx]
            smb_ice[idx] = domain.smb_ice[idx]
            liquid_water[idx] = liquid_water_local
            runoff[idx] = domain.runoff[idx]
        end
        return nothing
    elseif backend == :kernelabstractions
        kernel! = _summarize_domain_state_kernel!(_ka_backend(domain.mass))
        event = kernel!(
            thickness,
            wet_mass,
            bulk_density,
            base_mass,
            smb_ice,
            liquid_water,
            runoff,
            domain.N,
            domain.mass,
            domain.mass_w,
            domain.density,
            domain.mass_base,
            domain.smb_ice,
            domain.runoff;
            ndrange=column_count(domain),
        )
        _wait_kernel(event)
        return nothing
    end
    error("Unsupported summary backend `$backend`.")
end

"""
    summarize_domain_state(domain; backend=:threads)

Allocate and return a named tuple of per-column summary arrays for `domain`.
This is a convenience wrapper around `summarize_domain_state!`.
"""
function summarize_domain_state(domain::AbstractSnowpackDomain; backend::Symbol=:threads)
    ncol = column_count(domain)
    NF = eltype(domain)
    allocate() = similar(domain.mass, NF, ncol)
    thickness = allocate()
    wet_mass = allocate()
    bulk_density = allocate()
    base_mass = allocate()
    smb_ice = allocate()
    liquid_water = allocate()
    runoff = allocate()
    summarize_domain_state!(
        thickness,
        wet_mass,
        bulk_density,
        base_mass,
        smb_ice,
        liquid_water,
        runoff,
        domain;
        backend=backend,
    )
    return (
        thickness=thickness,
        wet_mass=wet_mass,
        bulk_density=bulk_density,
        base_mass=base_mass,
        smb_ice=smb_ice,
        liquid_water=liquid_water,
        runoff=runoff,
    )
end

"""
    summarize_cycle_state!(thickness, wet_mass, bulk_density, base_mass, domain; backend=:threads)

Fill preallocated arrays with the smaller summary set used to compare
equilibrium cycles. Mutates the output arrays and returns `nothing`.
"""
function summarize_cycle_state!(
    thickness::AbstractVector,
    wet_mass::AbstractVector,
    bulk_density::AbstractVector,
    base_mass::AbstractVector,
    domain::AbstractSnowpackDomain;
    backend::Symbol=:threads,
)
    if backend == :threads
        Threads.@threads for idx in 1:column_count(domain)
            n = domain.N[idx]
            thickness_local = zero(eltype(thickness))
            wet_mass_local = zero(eltype(wet_mass))
            solid_mass_local = zero(eltype(wet_mass))
            @inbounds for layer_index in 1:n
                solid = domain.mass[layer_index, idx]
                liquid = domain.mass_w[layer_index, idx]
                rho = domain.density[layer_index, idx]
                solid_mass_local += solid
                wet_mass_local += solid + liquid
                if solid > zero(solid) && rho > EPS_TINY
                    thickness_local += solid / rho
                end
            end
            thickness[idx] = thickness_local
            wet_mass[idx] = wet_mass_local
            bulk_density[idx] = thickness_local > EPS_TINY ? solid_mass_local / thickness_local : zero(eltype(bulk_density))
            base_mass[idx] = domain.mass_base[idx]
        end
        return nothing
    elseif backend == :kernelabstractions
        kernel! = _summarize_cycle_state_kernel!(_ka_backend(domain.mass))
        event = kernel!(
            thickness,
            wet_mass,
            bulk_density,
            base_mass,
            domain.N,
            domain.mass,
            domain.mass_w,
            domain.density,
            domain.mass_base;
            ndrange=column_count(domain),
        )
        _wait_kernel(event)
        return nothing
    end
    error("Unsupported summary backend `$backend`.")
end

"""
    summarize_cycle_state(domain; backend=:threads)

Allocate and return a named tuple of cycle-level summary arrays for `domain`.
This is the allocating counterpart to `summarize_cycle_state!`.
"""
function summarize_cycle_state(domain::AbstractSnowpackDomain; backend::Symbol=:threads)
    ncol = column_count(domain)
    NF = eltype(domain)
    allocate() = similar(domain.mass, NF, ncol)
    thickness = allocate()
    wet_mass = allocate()
    bulk_density = allocate()
    base_mass = allocate()
    summarize_cycle_state!(thickness, wet_mass, bulk_density, base_mass, domain; backend=backend)
    return (
        thickness=thickness,
        wet_mass=wet_mass,
        bulk_density=bulk_density,
        base_mass=base_mass,
    )
end
