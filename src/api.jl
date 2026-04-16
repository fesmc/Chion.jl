"""
Lightweight public run API.
"""

@inline function _normalize_run_save(save)
    save === nothing && return Symbol[]
    return normalize_netcdf_variables(save)
end

"""
    run!(domain, forcing; layout=nothing, save=[], backend=:threads, ...)

Run Chion directly from a prepared domain and forcing container.

`domain` already carries the physics constants, so the common flow is:

1. construct a `SnowpackDomain` with `c=physics(...)`
2. prepare a `ForcingData`
3. choose which fields to save via `save`
4. call `run!`

`save` accepts exact output symbols, output groups such as `:final` or
`:history`, or any mixture of both. When `save` is empty, no NetCDF file is
written.
"""
function run!(problem::LoadedProblem; kwargs...)
    return run!(problem.domain, problem.forcing; layout=problem.layout, kwargs...)
end

function run!(
    domain::SnowpackDomain,
    forcing::ForcingData;
    layout::Union{Nothing, GridLayout}=nothing,
    save=Symbol[],
    cycles::Integer=1,
    backend=:threads,
    write_outputs::Bool=false,
    output_dir::AbstractString="",
    netcdf_path::AbstractString="",
    history_stride::Integer=1,
    io::IO=stdout,
    timings::StepTimingStats=StepTimingStats(),
    run_wall_t0::Integer=time_ns(),
)
    save_variables = _normalize_run_save(save)
    options = RunOptions(
        name="chion_run",
        input_label="",
        cycles=Int(cycles),
        backend=normalize_backend(backend),
        write_outputs=write_outputs,
        output_dir=String(output_dir),
        netcdf_path=String(netcdf_path),
        write_netcdf=!isempty(save_variables),
        netcdf_variables=save_variables,
        history_stride=Int(history_stride),
    )
    return execute_run!(
        domain,
        forcing;
        layout=layout,
        options=options,
        io=io,
        timings=timings,
        run_wall_t0=run_wall_t0,
    )
end
