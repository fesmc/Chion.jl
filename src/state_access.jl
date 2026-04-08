"""
State accessors and formatted state output.
"""

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
    albedo_dynamic,
    idx::Int,
    c::SnowpackPhysicalConstants,
)
    snow_cover = _snow_cover_fraction(N_storage, mass, mass_w, density, idx)
    n = _n_active(N_storage, idx)
    if n == 0
        return Dict(
            "N" => 0,
            "n_active" => 0,
            "mass" => Float64[],
            "solid_mass" => Float64[],
            "mass_w" => Float64[],
            "liquid_water_mass" => Float64[],
            "density" => Float64[],
            "total_mass" => 0.0,
            "total_liquid_water" => 0.0,
            "total_wet_mass" => 0.0,
            "thickness" => Float64[],
            "total_thickness" => 0.0,
            "surface_temperature" => c.T0,
            "snow_cover" => snow_cover,
            "surface_albedo" => _get_scalar(albedo_dynamic, idx),
            "albedo_dynamic" => _get_scalar(albedo_dynamic, idx),
            "smb_ice" => _get_scalar(smb_ice, idx),
            "ice_sheet_smb" => _get_scalar(smb_ice, idx),
        )
    end

    active_solid_mass = [@inbounds _get_layer(mass, layer_index, idx) for layer_index in 1:n]
    active_liquid_water_mass = [@inbounds _get_layer(mass_w, layer_index, idx) for layer_index in 1:n]
    active_density = [@inbounds _get_layer(density, layer_index, idx) for layer_index in 1:n]
    thickness = active_solid_mass ./ active_density

    return Dict(
        "N" => n,
        "n_active" => n,
        "mass" => active_solid_mass,
        "solid_mass" => active_solid_mass,
        "mass_w" => active_liquid_water_mass,
        "liquid_water_mass" => active_liquid_water_mass,
        "density" => active_density,
        "total_mass" => sum(active_solid_mass),
        "total_liquid_water" => sum(active_liquid_water_mass),
        "total_wet_mass" => sum(active_solid_mass) + sum(active_liquid_water_mass),
        "thickness" => thickness,
        "total_thickness" => sum(thickness),
        "surface_temperature" => _get_layer(temperature, 1, idx),
        "snow_cover" => snow_cover,
        "surface_albedo" => _get_scalar(albedo_dynamic, idx),
        "albedo_dynamic" => _get_scalar(albedo_dynamic, idx),
        "smb_ice" => _get_scalar(smb_ice, idx),
        "ice_sheet_smb" => _get_scalar(smb_ice, idx),
    )
end

"""
    get_state(domain, idx=1)

Return a dictionary snapshot for column `idx` of `domain`. Values are copied
into plain Julia containers so callers can inspect state without mutating the
domain.
"""
function get_state(domain::AbstractSnowpackDomain, idx::Int=1)
    return _state_dict(
        domain.N,
        domain.mass,
        domain.mass_w,
        domain.density,
        domain.temperature,
        domain.smb_ice,
        domain.albedo_dynamic,
        idx,
        domain.c,
    )
end

"""
    print_state(domain, idx=1)

Print a short formatted summary of column `idx` to standard output. This is a
diagnostic convenience wrapper around `get_state`.
"""
function print_state(domain::AbstractSnowpackDomain, idx::Int=1)
    state = get_state(domain, idx)
    println("=" ^ 60)
    println("Snowpack Domain Column State")
    println("=" ^ 60)
    println("Column index: ", idx)
    println("Active layers: ", state["N"])
    println("Total mass: ", round(state["total_mass"], digits=2), " kg/m^2")
    println("Total thickness: ", round(state["total_thickness"], digits=3), " m")
    println("Snow cover: ", round(state["snow_cover"], digits=3))
    println("Surface albedo: ", round(state["surface_albedo"], digits=3))
    println()
end

"""
    compute_auxiliary!(domain, idx)

Recompute derived diagnostics for column `idx` in-place. Currently this
updates snow-cover fraction and surface albedo.
"""
function compute_auxiliary!(domain::AbstractSnowpackDomain, idx::Int)
    update_snow_cover!(domain, idx)
    update_surface_albedo!(domain, idx)
    return nothing
end

"""
    compute_auxiliary!(domain)

Recompute derived diagnostics for every column in `domain`. Mutates the
domain’s auxiliary fields in-place and returns `nothing`.
"""
function compute_auxiliary!(domain::AbstractSnowpackDomain)
    for idx in 1:column_count(domain)
        compute_auxiliary!(domain, idx)
    end
    return nothing
end
