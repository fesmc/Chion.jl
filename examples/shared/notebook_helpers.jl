module ChionNotebookHelpers

using Dates
using CUDA
using Chion
const HAS_PLOTS = try
    @eval import Plots
    true
catch
    false
end

const SM = Chion.SnowpackModel

export pkg_root
export result_domain_cpu, summarize_column, output_dir_for
export cuda_preflight, device_report, run_case_capture
export plots_available, forcing_timeseries_plot, history_plot
export column_profile_plot, domain_metric_values, layout_heatmap_plot

pkg_root() = normpath(joinpath(@__DIR__, "..", ".."))

function output_dir_for(slug::AbstractString)
    out_dir = joinpath(pkg_root(), "plots", "pluto", slug)
    mkpath(out_dir)
    return out_dir
end

function run_case_capture(case::Chion.SnowpackCase; copy_domain::Bool=true)
    return mktemp() do _, io
        result = Chion.run_case(case; io=io, copy_domain=copy_domain)
        flush(io)
        seekstart(io)
        return (
            result=result,
            log=read(io, String),
        )
    end
end

plots_available() = HAS_PLOTS

function _plots_module()
    HAS_PLOTS || error("Plots.jl is not available. Install it in the environment that launches Pluto to enable plotting cells.")
    try
        Plots.default(fmt=:svg)
    catch
        nothing
    end
    return Plots
end

result_domain_cpu(result::Chion.RunResult) = result_domain_cpu(result.domain)

function result_domain_cpu(domain::SM.SnowpackDomain)
    return domain.mass isa Array ? domain : SM.cpu_domain(domain)
end

function summarize_column(result_or_domain, idx::Integer=1)
    domain = result_or_domain isa Chion.RunResult ? result_domain_cpu(result_or_domain) : result_domain_cpu(result_or_domain)
    state = Chion.get_state(domain, Int(idx))
    return (
        idx=Int(idx),
        n_active=state["n_active"],
        total_mass=state["total_mass"],
        total_liquid_water=state["total_liquid_water"],
        total_wet_mass=state["total_wet_mass"],
        total_thickness=state["total_thickness"],
        surface_temperature_c=state["surface_temperature"] - domain.c.T0,
        snow_cover=state["snow_cover"],
        surface_albedo=state["surface_albedo"],
        smb_ice=state["smb_ice"],
    )
end

function _finite_extrema(field)
    finite_values = Float64[]
    sizehint!(finite_values, length(field))
    for value in field
        if isfinite(value)
            push!(finite_values, Float64(value))
        end
    end
    isempty(finite_values) && return (-1.0, 1.0)
    vmin = minimum(finite_values)
    vmax = maximum(finite_values)
    vmin == vmax && return (vmin - 1.0, vmax + 1.0)
    return (vmin, vmax)
end

function _symmetric_clims(field)
    vmax = 0.0
    found = false
    for value in field
        if isfinite(value)
            vmax = max(vmax, abs(Float64(value)))
            found = true
        end
    end
    found || return (-1.0, 1.0)
    vmax > 0.0 || return (-1.0, 1.0)
    return (-vmax, vmax)
end

function _layout_grid(layout::Chion.GridLayout, values::AbstractVector{<:Real})
    length(values) == length(layout.js) || error("Value count must match the layout point count.")
    grid = fill(NaN, size(layout.mask))
    @inbounds for idx in eachindex(values)
        grid[layout.js[idx], layout.is[idx]] = Float64(values[idx])
    end
    return grid
end

function _domain_summary(result_or_domain)
    domain = result_domain_cpu(result_or_domain)
    return SM.summarize_domain_state(domain; backend=:threads)
end

function domain_metric_values(result_or_domain, metric::Symbol)
    summary = _domain_summary(result_or_domain)
    if metric == :thickness
        return Vector{Float64}(summary.thickness)
    elseif metric == :wet_mass
        return Vector{Float64}(summary.wet_mass)
    elseif metric == :bulk_density
        return Vector{Float64}(summary.bulk_density)
    elseif metric == :base_mass
        return Vector{Float64}(summary.base_mass)
    elseif metric == :smb_ice
        return Vector{Float64}(summary.smb_ice)
    elseif metric == :liquid_water
        return Vector{Float64}(summary.liquid_water)
    elseif metric == :runoff
        return Vector{Float64}(summary.runoff)
    end
    error("Unsupported metric '$metric'.")
end

function forcing_timeseries_plot(
    forcing::Chion.ForcingData;
    idx::Integer=1,
    title_prefix::AbstractString="Forcing overview",
)
    P = _plots_module()
    air_temperature_c = forcing.air_temperature[Int(idx), :] .- 273.15
    snowfall_mm_day = forcing.snowfall_rate[Int(idx), :] .* 86_400.0
    rainfall_mm_day = forcing.rainfall_rate[Int(idx), :] .* 86_400.0
    p1 = P.plot(
        forcing.time_values,
        air_temperature_c;
        lw=3,
        color=:steelblue,
        marker=:circle,
        xlabel="Time",
        ylabel="C",
        title="$(title_prefix): air temperature",
        legend=false,
        framestyle=:box,
    )
    p2 = P.plot(
        forcing.time_values,
        snowfall_mm_day;
        lw=3,
        color=:royalblue,
        marker=:circle,
        label="snow",
        xlabel="Time",
        ylabel="mmWE/day",
        title="$(title_prefix): snowfall and rainfall",
        framestyle=:box,
    )
    P.plot!(
        p2,
        forcing.time_values,
        rainfall_mm_day;
        lw=3,
        color=:firebrick,
        marker=:diamond,
        label="rain",
    )
    p3 = P.plot(
        forcing.time_values,
        forcing.shortwave_down[Int(idx), :];
        lw=3,
        color=:darkorange,
        marker=:circle,
        xlabel="Time",
        ylabel="W/m^2",
        title="$(title_prefix): shortwave down",
        legend=false,
        framestyle=:box,
    )
    p4 = P.plot(
        forcing.time_values,
        forcing.wind_speed[Int(idx), :];
        lw=3,
        color=:seagreen,
        marker=:circle,
        xlabel="Time",
        ylabel="m/s",
        title="$(title_prefix): wind speed",
        legend=false,
        framestyle=:box,
    )
    return P.plot(p1, p2, p3, p4; layout=(2, 2), size=(950, 650))
end

function history_plot(history::Vector{<:NamedTuple}; title::AbstractString="Cycle history")
    isempty(history) && error("History is empty. Run at least one cycle before plotting cycle history.")
    P = _plots_module()
    cycles = getproperty.(history, :cycle)
    mean_thickness = getproperty.(history, :mean_thickness)
    mean_wet_mass = getproperty.(history, :mean_wet_mass)
    mean_base_mass = getproperty.(history, :mean_base_mass)
    mean_abs_dth = getproperty.(history, :mean_abs_delta_thickness)
    mean_abs_dswe = getproperty.(history, :mean_abs_delta_wet_mass)
    mean_abs_dbase = getproperty.(history, :mean_abs_delta_base_mass)
    p1 = P.plot(cycles, mean_thickness; lw=3, marker=:circle, color=:steelblue, xlabel="Cycle", ylabel="m", title="Mean thickness", framestyle=:box, legend=false)
    p2 = P.plot(cycles, mean_wet_mass; lw=3, marker=:circle, color=:forestgreen, xlabel="Cycle", ylabel="mmWE", title="Mean wet mass", framestyle=:box, legend=false)
    p3 = P.plot(cycles, mean_base_mass; lw=3, marker=:circle, color=:purple, xlabel="Cycle", ylabel="mmWE", title="Mean base mass", framestyle=:box, legend=false)
    p4 = P.plot(cycles, mean_abs_dth; lw=3, marker=:circle, color=:firebrick, xlabel="Cycle", ylabel="m", title="Mean abs dThickness", framestyle=:box, legend=false)
    p5 = P.plot(cycles, mean_abs_dswe; lw=3, marker=:circle, color=:darkorange, xlabel="Cycle", ylabel="mmWE", title="Mean abs dSWE", framestyle=:box, legend=false)
    p6 = P.plot(cycles, mean_abs_dbase; lw=3, marker=:circle, color=:indigo, xlabel="Cycle", ylabel="mmWE", title="Mean abs dBase", framestyle=:box, legend=false)
    return P.plot(p1, p2, p3, p4, p5, p6; layout=(2, 3), size=(1150, 700), plot_title=title)
end

function column_profile_plot(result_or_domain, idx::Integer=1; title::AbstractString="Final column profile")
    P = _plots_module()
    domain = result_domain_cpu(result_or_domain)
    state = Chion.get_state(domain, Int(idx))
    layer_density = Float64.(state["density"])
    layer_mass = Float64.(state["mass"])
    n = length(layer_density)
    layers = collect(1:n)
    p1 = P.bar(
        layers,
        layer_mass;
        color=:steelblue,
        xlabel="Layer",
        ylabel="kg/m^2",
        title="$(title): layer mass",
        framestyle=:box,
        legend=false,
    )
    p2 = P.plot(
        layers,
        layer_density;
        lw=3,
        marker=:circle,
        color=:firebrick,
        xlabel="Layer",
        ylabel="kg/m^3",
        title="$(title): layer density",
        framestyle=:box,
        legend=false,
    )
    return P.plot(p1, p2; layout=(1, 2), size=(900, 350))
end

function layout_heatmap_plot(
    layout::Chion.GridLayout,
    values::AbstractVector{<:Real};
    title::AbstractString,
    unit::AbstractString="",
    color=:viridis,
    symmetric::Bool=false,
)
    P = _plots_module()
    grid = _layout_grid(layout, values)
    clims = symmetric ? _symmetric_clims(grid) : _finite_extrema(grid)
    plot_color = symmetric && color == :viridis ? P.cgrad([:navy, :white, :firebrick]) : color
    return P.heatmap(
        layout.x,
        layout.y,
        grid;
        title=title,
        xlabel="x",
        ylabel="y",
        aspect_ratio=:equal,
        color=plot_color,
        clims=clims,
        colorbar_title=unit,
        framestyle=:box,
    )
end

function cuda_preflight()
    functional = try
        SM.cuda_available()
    catch
        false
    end
    message = functional ?
        "CUDA.functional() is true. You can use `backend=:gpu`." :
        "CUDA.functional() is false. Request a GPU, load the CUDA module, and use a writable depot if precompilation fails."
    return (
        functional=functional,
        recommended_modules=("julia/1.12.2", "cuda/13.1.0"),
        depot_hint="export JULIA_DEPOT_PATH=/tmp/chion-pluto:\$HOME/.julia",
        message=message,
    )
end

function device_report()
    preflight = cuda_preflight()
    if !preflight.functional
        return (
            functional=false,
            device_name="unavailable",
            capability=nothing,
            total_memory_bytes=nothing,
            free_memory_bytes=nothing,
            pool_status=preflight.message,
        )
    end

    device = CUDA.device()
    device_name = try
        String(CUDA.name(device))
    catch
        sprint(show, device)
    end
    capability = try
        CUDA.capability(device)
    catch
        nothing
    end
    total_memory_bytes = try
        Int(CUDA.totalmem(device))
    catch
        nothing
    end
    free_memory_bytes = try
        Int(CUDA.available_memory())
    catch
        nothing
    end
    pool_status = try
        sprint(io -> CUDA.pool_status(io))
    catch err
        "CUDA.pool_status unavailable: $(sprint(showerror, err))"
    end

    return (
        functional=true,
        device_name=device_name,
        capability=capability,
        total_memory_bytes=total_memory_bytes,
        free_memory_bytes=free_memory_bytes,
        pool_status=pool_status,
    )
end

end
