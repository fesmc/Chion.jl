"""Runtime container types shared by simulation, stepping, diagnostics, and output."""

struct IntegratorClocks
    run_wall_t0::Int
    simulation_wall_t0::Int
end

mutable struct ModelRuntime
    backend
    ncol::Int
    grid
    active::Vector{Bool}
    active_indices
end

mutable struct SimulationIntegrator
    sim
    options::RunOptions
    io::IO
    timings::StepTimingStats
    clocks::IntegratorClocks
    model_runtime::ModelRuntime
    diagnostics
    output
    time_index::Int
    completed_years::Int
    progress
    finalized::Bool
    result::Union{Nothing, SimulationResult}
end

@inline function Base.getproperty(integrator::SimulationIntegrator, name::Symbol)
    name === :forcing && return getfield(integrator, :sim).forcing
    return getfield(integrator, name)
end

function Base.propertynames(integrator::SimulationIntegrator, private::Bool=false)
    names = fieldnames(typeof(integrator))
    return private ? names : (names..., :forcing)
end
