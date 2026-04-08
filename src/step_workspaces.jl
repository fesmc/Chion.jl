"""
Scratch storage for stepping and energy-flux solves.
"""

struct EnergyWorkspace{
    LT,
    DT,
    UT,
    RT,
    IT,
    PT,
    TT,
    KT,
}
    lower::LT
    diag::DT
    upper::UT
    rhs::RT
    interface_conductance::IT
    previous_temperature::PT
    layer_thickness::TT
    thermal_conductivity::KT
end

function EnergyWorkspace(::Type{NF}, Ntot::Int) where {NF <: AbstractFloat}
    return EnergyWorkspace(
        zeros(NF, Ntot),
        zeros(NF, Ntot),
        zeros(NF, Ntot),
        zeros(NF, Ntot),
        zeros(NF, Ntot),
        zeros(NF, Ntot),
        zeros(NF, Ntot),
        zeros(NF, Ntot),
    )
end

function EnergyWorkspace(::Type{NF}, Ntot::Int, ncol::Int) where {NF <: AbstractFloat}
    return EnergyWorkspace(
        zeros(NF, Ntot, ncol),
        zeros(NF, Ntot, ncol),
        zeros(NF, Ntot, ncol),
        zeros(NF, Ntot, ncol),
        zeros(NF, Ntot, ncol),
        zeros(NF, Ntot, ncol),
        zeros(NF, Ntot, ncol),
        zeros(NF, Ntot, ncol),
    )
end

function EnergyWorkspace(storage, ::Type{NF}, Ntot::Int) where {NF <: AbstractFloat}
    allocate() = similar(storage, NF, Ntot)
    return EnergyWorkspace(
        allocate(),
        allocate(),
        allocate(),
        allocate(),
        allocate(),
        allocate(),
        allocate(),
        allocate(),
    )
end

EnergyWorkspace(domain::AbstractSnowpackDomain) = EnergyWorkspace(domain.mass, number_type(domain.c), domain.Ntot)

function EnergyWorkspace(storage, ::Type{NF}, Ntot::Int, ncol::Int) where {NF <: AbstractFloat}
    allocate() = similar(storage, NF, Ntot, ncol)
    return EnergyWorkspace(
        allocate(),
        allocate(),
        allocate(),
        allocate(),
        allocate(),
        allocate(),
        allocate(),
        allocate(),
    )
end

struct StepWorkspace{LWT,ET}
    liquid_water_before_energy::LWT
    energy::ET
end

function StepWorkspace(::Type{NF}, Ntot::Int) where {NF <: AbstractFloat}
    return StepWorkspace(
        zeros(NF, Ntot),
        EnergyWorkspace(NF, Ntot),
    )
end

function StepWorkspace(domain::AbstractSnowpackDomain)
    NF = number_type(domain.c)
    return StepWorkspace(
        similar(domain.mass, NF, domain.Ntot),
        EnergyWorkspace(domain.mass, NF, domain.Ntot),
    )
end

function threaded_workspaces(domain::AbstractSnowpackDomain)
    return [StepWorkspace(domain) for _ in 1:Threads.maxthreadid()]
end

struct ColumnarStepWorkspace{LWT,ET}
    liquid_water_before_energy::LWT
    energy::ET
end

function ColumnarStepWorkspace(::Type{NF}, Ntot::Int, ncol::Int) where {NF <: AbstractFloat}
    return ColumnarStepWorkspace(
        zeros(NF, Ntot, ncol),
        EnergyWorkspace(NF, Ntot, ncol),
    )
end

function ColumnarStepWorkspace(domain::AbstractSnowpackDomain)
    NF = number_type(domain.c)
    ncol = column_count(domain)
    return ColumnarStepWorkspace(
        similar(domain.mass, NF, domain.Ntot, ncol),
        EnergyWorkspace(domain.mass, NF, domain.Ntot, ncol),
    )
end

@inline function column_workspace(workspace::ColumnarStepWorkspace, idx::Int)
    @views return StepWorkspace(
        workspace.liquid_water_before_energy[:, idx],
        EnergyWorkspace(
            workspace.energy.lower[:, idx],
            workspace.energy.diag[:, idx],
            workspace.energy.upper[:, idx],
            workspace.energy.rhs[:, idx],
            workspace.energy.interface_conductance[:, idx],
            workspace.energy.previous_temperature[:, idx],
            workspace.energy.layer_thickness[:, idx],
            workspace.energy.thermal_conductivity[:, idx],
        ),
    )
end

@inline function _launch_step_columns_kernel!(
    kernel!,
    domain::AbstractSnowpackDomain,
    workspace::ColumnarStepWorkspace,
    forcing::SnowpackStepFields,
    dt_value,
    trailing_args...,
)
    return kernel!(
        domain.N,
        domain.mass,
        domain.mass_w,
        domain.density,
        domain.temperature,
        domain.mass_base,
        domain.smb_ice,
        domain.runoff,
        domain.Tsrf,
        domain.snow_cover,
        domain.albedo_dynamic,
        domain.c,
        domain.Ntot,
        domain.mass_max,
        domain.mass_split,
        domain.mass_min,
        domain.f_base_max,
        workspace,
        forcing.air_temperature,
        forcing.snowfall_rate,
        forcing.rainfall_rate,
        forcing.shortwave_down,
        forcing.wind_speed,
        forcing.q_lw_down,
        forcing.has_q_lw_down,
        forcing.q_sh,
        forcing.has_q_sh,
        forcing.q_lh,
        forcing.has_q_lh,
        dt_value,
        trailing_args...;
        ndrange=column_count(domain),
    )
end

Adapt.@adapt_structure EnergyWorkspace
Adapt.@adapt_structure StepWorkspace
Adapt.@adapt_structure ColumnarStepWorkspace
