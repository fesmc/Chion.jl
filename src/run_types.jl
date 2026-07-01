"""Runtime container types shared by simulation, stepping, diagnostics, and output."""

mutable struct ModelRuntime
    backend
    active::Vector{Bool}
    active_indices
end

mutable struct SimulationIntegrator
    sim
    io::IO
    timings::StepTimingStats
    wall_t0::Int
    model_runtime::ModelRuntime
    history::Vector{NamedTuple}
    netcdf_path::String
    time_index::Int
    completed_years::Int
    result::Union{Nothing, SimulationResult}
end
