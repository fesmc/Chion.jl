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

mutable struct NativeOutput <: AbstractOutput
    schedule
    netcdf::Union{Nothing, NetcdfOutput}
    nc_path::String
    monthly::MonthlyState
    step_vectors
    daily_vectors
    daily_written::Int
    steps_written::Int
end

mutable struct SimulationIntegrator
    sim
    options::RunOptions
    io::IO
    timings::StepTimingStats
    clocks::IntegratorClocks
    model_runtime::ModelRuntime
    diagnostics
    output::NativeOutput
    time_index::Int
    completed_years::Int
    current_forcing::SnowpackForcing
    progress
    finalized::Bool
    result::Union{Nothing, SimulationResult}
end
