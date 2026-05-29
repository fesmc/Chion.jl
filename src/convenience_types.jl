"""Small option/result containers shared by simulation, I/O, and reporting."""

struct RunOptions
    name::String
    input_label::String
    output_dir::String
    netcdf_path::String
    write_outputs::Bool
    write_netcdf::Bool
    netcdf_variables::Vector{Symbol}
    years::Int
    backend::Symbol
    history_year_stride::Int
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
    write_outputs::Bool=true,
    write_netcdf::Bool=true,
    netcdf_variables=copy(NETCDF_VARIABLES),
    years::Integer=10,
    backend=:threads,
    history_year_stride::Integer=1,
)
    resolved_name = String(name)
    resolved_output_dir = isempty(output_dir) ? _default_output_dir(resolved_name) : String(output_dir)
    return RunOptions(
        resolved_name,
        String(input_label),
        resolved_output_dir,
        String(netcdf_path),
        write_outputs,
        write_netcdf,
        normalize_netcdf_variables(netcdf_variables),
        normalize_years(years),
        normalize_backend(backend),
        normalize_history_year_stride(history_year_stride),
    )
end

"""
    SimulationResult

Summary returned by `finalize!` and `run!`.
"""
struct SimulationResult
    history::Vector{NamedTuple}
    status::Symbol
    years_completed::Int
    timings::StepTimingStats
    simulation_wall_sec::Float64
    run_wall_sec::Float64
    netcdf_path::String
    summary_path::String
    history_csv_path::String
end
