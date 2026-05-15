"""State containers owned by `Simulation.ref` and `Simulation.now`."""

"""Abstract supertype for mutable model state owned by `Simulation`."""
abstract type AbstractSnowModelState end

"""Mutable state for `BESSIModel`."""
struct BESSIState{D <: SnowpackDomain} <: AbstractSnowModelState
    domain::D
end

function BESSIState(model::BESSIModel)
    return BESSIState(SnowpackDomain(;
        c=model.c,
        Ntot=model.Ntot,
        ncol=ncols(model.grid),
        mass_max=model.mass_max,
        mass_split=model.mass_split,
        mass_min=model.mass_min,
        rho_max=model.rho_max,
        density_init=model.density_init,
        temperature_init=model.temperature_init,
    ))
end

"""Mutable state for `PDDModel`."""
struct PDDState{ST <: AbstractVector{Float64}} <: AbstractSnowModelState
    snowpack_swe::ST
    smb_ice::ST
    runoff::ST
    pdd_sum::ST
end

function PDDState(model::PDDModel)
    ncol = ncols(model.grid)
    return PDDState(
        zeros(Float64, ncol),
        zeros(Float64, ncol),
        zeros(Float64, ncol),
        zeros(Float64, ncol),
    )
end

"""Placeholder state for `ITMModel` until its physics are implemented."""
struct ITMState <: AbstractSnowModelState end

initial_state(model::BESSIModel) = BESSIState(model)
initial_state(model::PDDModel) = PDDState(model)
initial_state(::ITMModel) = ITMState()

function _copy_domain_state!(dest::SnowpackDomain, src::SnowpackDomain)
    dest.Ntot == src.Ntot || error("Cannot copy domain state with different `Ntot`.")
    dest.ncol == src.ncol || error("Cannot copy domain state with different column count.")
    dest.N .= src.N
    dest.mass .= src.mass
    dest.mass_w .= src.mass_w
    dest.density .= src.density
    dest.temperature .= src.temperature
    dest.mass_base .= src.mass_base
    dest.smb_ice .= src.smb_ice
    dest.runoff .= src.runoff
    dest.melt .= src.melt
    dest.refreezing .= src.refreezing
    dest.Tsrf .= src.Tsrf
    dest.snow_cover .= src.snow_cover
    dest.albedo_dynamic .= src.albedo_dynamic
    return dest
end
