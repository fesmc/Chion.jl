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

mutable struct DiagnosticsRuntime
    need_step_outputs::Bool
    need_monthly_outputs::Bool
    need_daily_outputs::Bool
    need_layer_outputs::Bool
    need_last_year_smb_delta::Bool
    need_step_diagnostics::Bool
    prev
    final
    backend_year_summary
    step_summary
    backend_step_summary
    deltas
    previous
    previous_year_smb_ice::Vector{Float64}
    history::Vector{NamedTuple}
end

mutable struct OutputRuntime
    schedule
    writer
    nc_path::String
    monthly_sums
    monthly_count::Vector{Int32}
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
    diagnostics::DiagnosticsRuntime
    output::OutputRuntime
    time_index::Int
    completed_years::Int
    current_forcing::SnowpackForcing
    progress
    finalized::Bool
    result::Union{Nothing, SimulationResult}
end
