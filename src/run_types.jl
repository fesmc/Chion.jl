"""Run configuration, result, and initialized runtime containers."""

const YearMetrics = @NamedTuple begin
    year::Int
    mean_thickness::Float64
    mean_wet_mass::Float64
    mean_bulk_density::Float64
    mean_base_mass::Float64
    mean_signed_delta_thickness::Float64
    mean_abs_delta_thickness::Float64
    max_abs_delta_thickness::Float64
    mean_signed_delta_wet_mass::Float64
    mean_abs_delta_wet_mass::Float64
    max_abs_delta_wet_mass::Float64
    mean_signed_delta_base_mass::Float64
    mean_abs_delta_base_mass::Float64
    max_abs_delta_base_mass::Float64
end

const YearHistory = Vector{YearMetrics}

struct RunOptions
    name::String
    input_label::String
    output_dir::String
    netcdf_path::String
    write_netcdf::Bool
    netcdf_variables::Vector{Symbol}
    years::Int
    backend::Symbol
    history_year_stride::Int
    compute_year_metrics::Bool
end

@inline function normalize_backend(backend)
    value = lowercase(strip(String(backend)))
    value == "cpu" && return :threads
    value in ("threads", "gpu") || error("Unsupported backend '$backend'. Use `threads`, `cpu`, or `gpu`.")
    return Symbol(value)
end

@inline normalize_history_year_stride(stride::Integer) =
    Int(stride) >= 0 ? Int(stride) : error("`history_year_stride` must be >= 0.")

@inline normalize_years(years::Integer) =
    Int(years) > 0 ? Int(years) : error("`years` must be positive.")

@inline should_record_year_metrics(year::Int, years::Int, stride::Int) =
    year == years || (stride > 0 && mod(year, stride) == 0)

function RunOptions(;
    name::AbstractString="chion_run",
    input_label::AbstractString="",
    output_dir::AbstractString="",
    netcdf_path::AbstractString="",
    write_netcdf::Bool=true,
    netcdf_variables=copy(NETCDF_VARIABLES),
    years::Integer=10,
    backend=:threads,
    history_year_stride::Integer=1,
    compute_year_metrics::Bool=true,
)
    resolved_name = String(name)
    resolved_output_dir = isempty(output_dir) ? _default_output_dir(resolved_name) : String(output_dir)
    return RunOptions(
        resolved_name,
        String(input_label),
        resolved_output_dir,
        String(netcdf_path),
        write_netcdf,
        normalize_netcdf_variables(netcdf_variables),
        normalize_years(years),
        normalize_backend(backend),
        normalize_history_year_stride(history_year_stride),
        compute_year_metrics,
    )
end

struct SimulationResult
    history::YearHistory
    status::Symbol
    years_completed::Int
    timings::StepTimingStats
    run_wall_sec::Float64
    netcdf_path::String
end

mutable struct ModelRuntime{D,A}
    data::D
    active::Vector{Bool}
    active_indices::A
end

struct ModelRuntimeData{S,F,W}
    state::S
    step_fields::F
    workspace::W
    is_gpu::Bool
end

mutable struct SimulationIntegrator{S,R<:ModelRuntime,I<:IO}
    sim::S
    io::I
    timings::StepTimingStats
    wall_t0::Int
    model_runtime::R
    history::YearHistory
    netcdf_path::String
    time_index::Int
    completed_years::Int
    result::Union{Nothing, SimulationResult}
end
