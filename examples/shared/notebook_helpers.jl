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

const SM = Chion

export pkg_root
export result_domain_cpu, summarize_column, output_dir_for
export cuda_preflight, device_report
export plots_available, forcing_timeseries_plot, history_plot
export column_profile_plot, domain_metric_values, layout_heatmap_plot
export seed_domain!, regular_grid

pkg_root() = normpath(joinpath(@__DIR__, "..", ".."))

function output_dir_for(slug::AbstractString)
    out_dir = joinpath(pkg_root(), "plots", "pluto", slug)
    mkpath(out_dir)
    return out_dir
end

"""
    seed_domain!(model; surface_mass, density, temperature_c) -> BESSIState

Initialise the first layer of every column in a new `BESSIState` with the given
surface mass, density, and temperature.
"""
function seed_domain!(
    model::Chion.BESSIModel;
    surface_mass,
    density::Real=model.c.rho_s,
    temperature_c::Real=-10.0,
)
    state = Chion.initial_state(model)
    return seed_domain!(state; surface_mass=surface_mass, density=density, temperature_c=temperature_c)
end

function seed_domain!(
    state::Chion.BESSIState;
    surface_mass,
    density::Real=state.domain.c.rho_s,
    temperature_c::Real=-10.0,
)
    domain = state.domain
    fill!(domain.N, 0)
    fill!(domain.mass, 0.0)
    fill!(domain.mass_w, 0.0)
    fill!(domain.density, 0.0)
    fill!(domain.temperature, domain.c.T0)
    fill!(domain.mass_base, 0.0)
    fill!(domain.smb_ice, 0.0)
    fill!(domain.runoff, 0.0)
    fill!(domain.snow_cover, 0.0)
    fill!(domain.albedo_dynamic, domain.c.alpha_dry)

    if surface_mass isa AbstractVector
        length(surface_mass) == domain.ncol || error("`surface_mass` must match the domain column count.")
        @inbounds for idx in eachindex(surface_mass)
            if surface_mass[idx] > 0
                domain.N[idx] = 1
                domain.mass[1, idx] = Float64(surface_mass[idx])
                domain.density[1, idx] = Float64(density)
                domain.temperature[1, idx] = domain.c.T0 + Float64(temperature_c)
                domain.Tsrf[idx] = domain.temperature[1, idx]
            else
                domain.Tsrf[idx] = domain.c.T0
            end
        end
    else
        @inbounds for idx in 1:domain.ncol
            if surface_mass > 0
                domain.N[idx] = 1
                domain.mass[1, idx] = Float64(surface_mass)
                domain.density[1, idx] = Float64(density)
                domain.temperature[1, idx] = domain.c.T0 + Float64(temperature_c)
                domain.Tsrf[idx] = domain.temperature[1, idx]
            else
                domain.Tsrf[idx] = domain.c.T0
            end
        end
    end

    Chion.compute_auxiliary!(domain)
    return state
end

"""
    regular_grid(nx, ny; x, y) -> SnowpackGrid

Build a `SnowpackGrid` for a regular `nx × ny` column layout.
"""
function regular_grid(nx::Integer, ny::Integer; x=collect(1:nx), y=collect(1:ny))
    nx > 0 || error("`nx` must be positive.")
    ny > 0 || error("`ny` must be positive.")
    xvals = Float64.(collect(x))
    yvals = Float64.(collect(y))
    length(xvals) == nx || error("`x` must have length `nx`.")
    length(yvals) == ny || error("`y` must have length `ny`.")
    js = [j for j in 1:ny for _ in 1:nx]
    is = [i for _ in 1:ny for i in 1:nx]
    mask = ones(Float64, ny, nx)
    return Chion.SnowpackGrid(nx * ny; x=xvals, y=yvals, js=js, is=is, mask=mask)
end

plots_available() = HAS_PLOTS

function _plots_module()
    HAS_PLOTS || error("Plots.jl is not available.")
    try Plots.default(fmt=:svg) catch end
    return Plots
end

function result_domain_cpu(::Chion.SimulationResult)
    error("Access the domain via simulation.now instead of the result.")
end

function result_domain_cpu(domain::SM.SnowpackDomain)
    return domain.mass isa Array ? domain : Chion.cpu_domain(domain)
end

result_domain_cpu(state::Chion.BESSIState) = result_domain_cpu(state.domain)
result_domain_cpu(simulation::Chion.Simulation) = result_domain_cpu(simulation.now)

function summarize_column(model_or_domain, idx::Integer=1)
    model_or_domain isa Chion.BESSIModel && error("BESSIModel is configuration-only; pass a BESSIState or Simulation.")
    domain = result_domain_cpu(model_or_domain)
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
        isfinite(value) && push!(finite_values, Float64(value))
    end
    isempty(finite_values) && return (-1.0, 1.0)
    vmin, vmax = minimum(finite_values), maximum(finite_values)
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

function _layout_grid(grid::Chion.SnowpackGrid, values::AbstractVector{<:Real})
    Chion.has_spatial_coords(grid) || error("Grid has no spatial coordinates.")
    length(values) == length(grid.js) || error("Value count must match the grid point count.")
    out = fill(NaN, size(grid.mask))
    @inbounds for idx in eachindex(values)
        out[grid.js[idx], grid.is[idx]] = Float64(values[idx])
    end
    return out
end

function _domain_summary(model_or_domain)
    model_or_domain isa Chion.BESSIModel && error("BESSIModel is configuration-only; pass a BESSIState or Simulation.")
    domain = result_domain_cpu(model_or_domain)
    return SM.summarize_domain_state(domain)
end

function domain_metric_values(model_or_domain, metric::Symbol)
    summary = _domain_summary(model_or_domain)
    metric == :thickness    && return Vector{Float64}(summary.thickness)
    metric == :wet_mass     && return Vector{Float64}(summary.wet_mass)
    metric == :bulk_density && return Vector{Float64}(summary.bulk_density)
    metric == :base_mass    && return Vector{Float64}(summary.base_mass)
    metric == :smb_ice      && return Vector{Float64}(summary.smb_ice)
    metric == :liquid_water && return Vector{Float64}(summary.liquid_water)
    metric == :runoff       && return Vector{Float64}(summary.runoff)
    error("Unsupported metric '$metric'.")
end

function forcing_timeseries_plot(
    forcing::Chion.SnowpackForcing;
    idx::Integer=1,
    title_prefix::AbstractString="Forcing overview",
)
    P = _plots_module()
    air_temperature_c = forcing.air_temperature[Int(idx), :] .- 273.15
    snowfall_mm_day   = forcing.snowfall_rate[Int(idx), :] .* 86_400.0
    rainfall_mm_day   = forcing.rainfall_rate[Int(idx), :] .* 86_400.0
    p1 = P.plot(forcing.time_values, air_temperature_c; lw=3, color=:steelblue, marker=:circle,
        xlabel="Time", ylabel="C", title="$(title_prefix): air temperature", legend=false, framestyle=:box)
    p2 = P.plot(forcing.time_values, snowfall_mm_day; lw=3, color=:royalblue, marker=:circle,
        label="snow", xlabel="Time", ylabel="mmWE/day", title="$(title_prefix): snowfall and rainfall", framestyle=:box)
    P.plot!(p2, forcing.time_values, rainfall_mm_day; lw=3, color=:firebrick, marker=:diamond, label="rain")
    p3 = P.plot(forcing.time_values, forcing.shortwave_down[Int(idx), :]; lw=3, color=:darkorange, marker=:circle,
        xlabel="Time", ylabel="W/m^2", title="$(title_prefix): shortwave down", legend=false, framestyle=:box)
    p4 = P.plot(forcing.time_values, forcing.wind_speed[Int(idx), :]; lw=3, color=:seagreen, marker=:circle,
        xlabel="Time", ylabel="m/s", title="$(title_prefix): wind speed", legend=false, framestyle=:box)
    return P.plot(p1, p2, p3, p4; layout=(2, 2), size=(950, 650))
end

function history_plot(history::Vector{<:NamedTuple}; title::AbstractString="Year history")
    isempty(history) && error("History is empty.")
    P = _plots_module()
    years          = getproperty.(history, :year)
    mean_thickness = getproperty.(history, :mean_thickness)
    mean_wet_mass  = getproperty.(history, :mean_wet_mass)
    mean_base_mass = getproperty.(history, :mean_base_mass)
    mean_abs_dth   = getproperty.(history, :mean_abs_delta_thickness)
    mean_abs_dswe  = getproperty.(history, :mean_abs_delta_wet_mass)
    mean_abs_dbase = getproperty.(history, :mean_abs_delta_base_mass)
    p1 = P.plot(years, mean_thickness;  lw=3, marker=:circle, color=:steelblue,   xlabel="Year", ylabel="m",    title="Mean thickness",      framestyle=:box, legend=false)
    p2 = P.plot(years, mean_wet_mass;   lw=3, marker=:circle, color=:forestgreen, xlabel="Year", ylabel="mmWE", title="Mean wet mass",        framestyle=:box, legend=false)
    p3 = P.plot(years, mean_base_mass;  lw=3, marker=:circle, color=:purple,      xlabel="Year", ylabel="mmWE", title="Mean base mass",       framestyle=:box, legend=false)
    p4 = P.plot(years, mean_abs_dth;   lw=3, marker=:circle, color=:firebrick,   xlabel="Year", ylabel="m",    title="Mean abs dThickness",  framestyle=:box, legend=false)
    p5 = P.plot(years, mean_abs_dswe;  lw=3, marker=:circle, color=:darkorange,  xlabel="Year", ylabel="mmWE", title="Mean abs dSWE",        framestyle=:box, legend=false)
    p6 = P.plot(years, mean_abs_dbase; lw=3, marker=:circle, color=:indigo,      xlabel="Year", ylabel="mmWE", title="Mean abs dBase",       framestyle=:box, legend=false)
    return P.plot(p1, p2, p3, p4, p5, p6; layout=(2, 3), size=(1150, 700), plot_title=title)
end

function column_profile_plot(model_or_domain, idx::Integer=1; title::AbstractString="Final column profile")
    P = _plots_module()
    model_or_domain isa Chion.BESSIModel && error("BESSIModel is configuration-only; pass a BESSIState or Simulation.")
    domain = result_domain_cpu(model_or_domain)
    state = Chion.get_state(domain, Int(idx))
    layer_density = Float64.(state["density"])
    layer_mass    = Float64.(state["mass"])
    layers = collect(1:length(layer_density))
    p1 = P.bar(layers, layer_mass; color=:steelblue, xlabel="Layer", ylabel="kg/m^2",
        title="$(title): layer mass", framestyle=:box, legend=false)
    p2 = P.plot(layers, layer_density; lw=3, marker=:circle, color=:firebrick, xlabel="Layer",
        ylabel="kg/m^3", title="$(title): layer density", framestyle=:box, legend=false)
    return P.plot(p1, p2; layout=(1, 2), size=(900, 350))
end

function layout_heatmap_plot(
    grid::Chion.SnowpackGrid,
    values::AbstractVector{<:Real};
    title::AbstractString,
    unit::AbstractString="",
    color=:viridis,
    symmetric::Bool=false,
)
    P = _plots_module()
    g = _layout_grid(grid, values)
    clims = symmetric ? _symmetric_clims(g) : _finite_extrema(g)
    plot_color = symmetric && color == :viridis ? P.cgrad([:navy, :white, :firebrick]) : color
    return P.heatmap(grid.x, grid.y, g; title=title, xlabel="x", ylabel="y",
        aspect_ratio=:equal, color=plot_color, clims=clims,
        colorbar_title=unit, framestyle=:box)
end

function cuda_preflight()
    functional = try SM.cuda_available() catch; false end
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
    preflight.functional || return (
        functional=false, device_name="unavailable",
        capability=nothing, total_memory_bytes=nothing,
        free_memory_bytes=nothing, pool_status=preflight.message,
    )
    device = CUDA.device()
    return (
        functional=true,
        device_name=try String(CUDA.name(device)) catch; sprint(show, device) end,
        capability=try CUDA.capability(device) catch; nothing end,
        total_memory_bytes=try Int(CUDA.totalmem(device)) catch; nothing end,
        free_memory_bytes=try Int(CUDA.available_memory()) catch; nothing end,
        pool_status=try sprint(io -> CUDA.pool_status(io)) catch err
            "CUDA.pool_status unavailable: $(sprint(showerror, err))" end,
    )
end

end
