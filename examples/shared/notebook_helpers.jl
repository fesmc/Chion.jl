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
export build_physics, build_synthetic_case, build_run_options
export build_prescribed_case_data
export result_domain_cpu, summarize_column, output_dir_for
export cuda_preflight, device_report, run_case_capture
export plots_available, forcing_timeseries_plot, history_plot
export column_profile_plot, domain_metric_values, layout_heatmap_plot

pkg_root() = normpath(joinpath(@__DIR__, "..", ".."))

@inline _as_symbol(value::Symbol) = value
@inline _as_symbol(value) = Symbol(lowercase(strip(String(value))))

function build_physics(;
    albedo_scheme=:dynamic,
    densification=:bessi,
    fresh_snow_density=:constant,
)
    return Chion.SnowpackPhysicalConstants(
        Float64;
        albedo_scheme=_as_symbol(albedo_scheme),
        low_density_densification=_as_symbol(densification),
        fresh_snow_density_scheme=_as_symbol(fresh_snow_density),
    )
end

function build_synthetic_case(;
    variant::Symbol=:multi_column,
    physics::Chion.SnowpackPhysicalConstants{Float64}=build_physics(),
    ntot::Integer=5,
    ntime::Integer=12,
    nx::Union{Nothing, Integer}=nothing,
    ny::Union{Nothing, Integer}=nothing,
)
    data = Chion.load_forcing(
        Chion.synthetic_forcing(; variant=variant, ntime=ntime, nx=nx, ny=ny);
        physics=physics,
        ntot=ntot,
    )
    return deepcopy(data.domain), data.forcing, data.layout, data.metadata
end

function _regular_grid_layout(nx::Integer, ny::Integer)
    Int(nx) > 0 || error("`nx` must be positive.")
    Int(ny) > 0 || error("`ny` must be positive.")
    x = Float64[i - 1 for i in 1:Int(nx)]
    y = Float64[j - 1 for j in 1:Int(ny)]
    ncol = Int(nx) * Int(ny)
    js = Vector{Int}(undef, ncol)
    is = Vector{Int}(undef, ncol)
    idx = 1
    @inbounds for j in 1:Int(ny), i in 1:Int(nx)
        js[idx] = j
        is[idx] = i
        idx += 1
    end
    return Chion.EquilibriumGridLayout(x, y, js, is, ones(Float64, Int(ny), Int(nx)))
end

function _expand_column_vector(values, ncol::Int, name::AbstractString)
    if values isa Number
        return fill(Float64(values), ncol)
    end
    data = collect(values)
    if ndims(data) == 1
        length(data) == ncol || error("`$name` must be a scalar or a vector with length $ncol.")
        return Float64.(data)
    end
    error("`$name` must be a scalar or a vector with length $ncol.")
end

function _expand_numeric_timeseries(field, ncol::Int, ntime::Int, name::AbstractString)
    if field isa Number
        return fill(Float64(field), ncol, ntime)
    end
    data = collect(field)
    if ndims(data) == 1
        length(data) == ntime || error("`$name` must have length $ntime.")
        return repeat(reshape(Float64.(data), 1, ntime), ncol, 1)
    elseif ndims(data) == 2
        size(data) == (ncol, ntime) || error("`$name` must have size ($ncol, $ntime).")
        return Matrix{Float64}(data)
    end
    error("`$name` must be a scalar, a vector of length $ntime, or a matrix of size ($ncol, $ntime).")
end

function _expand_bool_timeseries(field, ncol::Int, ntime::Int, name::AbstractString)
    if field isa Bool
        return fill(field, ncol, ntime)
    end
    data = collect(field)
    if ndims(data) == 1
        length(data) == ntime || error("`$name` must have length $ntime.")
        return repeat(reshape(Bool.(data), 1, ntime), ncol, 1)
    elseif ndims(data) == 2
        size(data) == (ncol, ntime) || error("`$name` must have size ($ncol, $ntime).")
        return Bool.(data)
    end
    error("`$name` must be a Bool, a vector of length $ntime, or a matrix of size ($ncol, $ntime).")
end

function build_prescribed_case_data(;
    physics::Chion.SnowpackPhysicalConstants{Float64}=build_physics(),
    ntot::Integer=5,
    nx::Integer=1,
    ny::Integer=1,
    dt_days,
    air_temperature_c,
    snowfall_mm_day,
    rainfall_mm_day,
    shortwave_down,
    wind_speed=5.0,
    q_lw_down=nothing,
    has_q_lw_down=nothing,
    q_sh=nothing,
    has_q_sh=nothing,
    q_lh=nothing,
    has_q_lh=nothing,
    time_values=nothing,
    initial_surface_mass=250.0,
    initial_density=320.0,
    initial_temperature_c=-12.0,
    initial_albedo=physics.alpha_dry,
    forcing_label::AbstractString="prescribed_forcing",
)
    Int(ntot) > 0 || error("`ntot` must be positive.")
    dt_days_v = Float64.(collect(dt_days))
    isempty(dt_days_v) && error("`dt_days` must not be empty.")
    all(>(0.0), dt_days_v) || error("All `dt_days` entries must be positive.")
    ntime = length(dt_days_v)
    layout = _regular_grid_layout(nx, ny)
    ncol = length(layout.js)

    air_temperature = _expand_numeric_timeseries(air_temperature_c, ncol, ntime, "air_temperature_c") .+ physics.T0
    snowfall_rate = _expand_numeric_timeseries(snowfall_mm_day, ncol, ntime, "snowfall_mm_day") ./ 86_400.0
    rainfall_rate = _expand_numeric_timeseries(rainfall_mm_day, ncol, ntime, "rainfall_mm_day") ./ 86_400.0
    shortwave_down_m = _expand_numeric_timeseries(shortwave_down, ncol, ntime, "shortwave_down")
    wind_speed_m = _expand_numeric_timeseries(wind_speed, ncol, ntime, "wind_speed")

    q_lw_down_m = isnothing(q_lw_down) ? zeros(Float64, ncol, ntime) : _expand_numeric_timeseries(q_lw_down, ncol, ntime, "q_lw_down")
    has_q_lw_down_m = isnothing(q_lw_down) ? fill(false, ncol, ntime) :
        (isnothing(has_q_lw_down) ? fill(true, ncol, ntime) : _expand_bool_timeseries(has_q_lw_down, ncol, ntime, "has_q_lw_down"))
    q_sh_m = isnothing(q_sh) ? zeros(Float64, ncol, ntime) : _expand_numeric_timeseries(q_sh, ncol, ntime, "q_sh")
    has_q_sh_m = isnothing(q_sh) ? fill(false, ncol, ntime) :
        (isnothing(has_q_sh) ? fill(true, ncol, ntime) : _expand_bool_timeseries(has_q_sh, ncol, ntime, "has_q_sh"))
    q_lh_m = isnothing(q_lh) ? zeros(Float64, ncol, ntime) : _expand_numeric_timeseries(q_lh, ncol, ntime, "q_lh")
    has_q_lh_m = isnothing(q_lh) ? fill(false, ncol, ntime) :
        (isnothing(has_q_lh) ? fill(true, ncol, ntime) : _expand_bool_timeseries(has_q_lh, ncol, ntime, "has_q_lh"))

    surface_mass = _expand_column_vector(initial_surface_mass, ncol, "initial_surface_mass")
    surface_density = _expand_column_vector(initial_density, ncol, "initial_density")
    surface_temperature_c = _expand_column_vector(initial_temperature_c, ncol, "initial_temperature_c")
    surface_albedo = _expand_column_vector(initial_albedo, ncol, "initial_albedo")

    N = Int[m > 0.0 ? 1 : 0 for m in surface_mass]
    mass = zeros(Float64, Int(ntot), ncol)
    mass_w = zeros(Float64, Int(ntot), ncol)
    density = zeros(Float64, Int(ntot), ncol)
    temperature = fill(physics.T0, Int(ntot), ncol)
    Tsrf = physics.T0 .+ surface_temperature_c
    snow_cover = Float64.(N .> 0)

    @inbounds for col in 1:ncol
        if N[col] > 0
            mass[1, col] = surface_mass[col]
            density[1, col] = surface_density[col]
            temperature[1, col] = physics.T0 + surface_temperature_c[col]
        end
    end

    state = Chion.SnowpackStateFields(
        N,
        mass,
        mass_w,
        density,
        temperature;
        Tsrf=Tsrf,
        snow_cover=snow_cover,
        albedo_dynamic=surface_albedo,
        physics=physics,
    )
    forcing = Chion.EquilibriumForcing(
        time_values=time_values,
        dt_days=dt_days_v,
        air_temperature=air_temperature,
        snowfall_rate=snowfall_rate,
        rainfall_rate=rainfall_rate,
        shortwave_down=shortwave_down_m,
        wind_speed=wind_speed_m,
        q_lw_down=q_lw_down_m,
        has_q_lw_down=has_q_lw_down_m,
        q_sh=q_sh_m,
        has_q_sh=has_q_sh_m,
        q_lh=q_lh_m,
        has_q_lh=has_q_lh_m,
    )
    metadata = (
        format=:prescribed,
        nx=Int(nx),
        ny=Int(ny),
        grid_shape=(Int(ny), Int(nx)),
        ncol=ncol,
        ntime=ntime,
        ntot=Int(ntot),
    )
    return Chion.SnowpackCaseData(
        Chion.SnowpackDomain(state),
        forcing;
        layout=layout,
        forcing_label=String(forcing_label),
        metadata=metadata,
    )
end

function output_dir_for(slug::AbstractString)
    out_dir = joinpath(pkg_root(), "plots", "pluto", slug)
    mkpath(out_dir)
    return out_dir
end

function build_run_options(;
    name::AbstractString="pluto_run",
    forcing_label::AbstractString="",
    backend=:cpu,
    out_dir::AbstractString=output_dir_for(String(name)),
    out_nc::AbstractString="",
    write_outputs::Bool=true,
    write_netcdf::Bool=false,
    netcdf_variables="final,history",
    max_cycles::Integer=3,
    cycle_metrics_stride::Integer=1,
)
    return Chion.EquilibriumRunOptions(
        name=String(name),
        forcing_label=String(forcing_label),
        out_dir=String(out_dir),
        out_nc=String(out_nc),
        write_outputs=write_outputs,
        write_netcdf=write_netcdf,
        netcdf_variables=netcdf_variables,
        max_cycles=Int(max_cycles),
        cycle_metrics_stride=Int(cycle_metrics_stride),
        backend=_as_symbol(backend),
    )
end

function run_case_capture(domain::SM.SnowpackDomain, forcing, layout, options::Chion.EquilibriumRunOptions)
    return mktemp() do _, io
        result = Chion.run_equilibrium!(deepcopy(domain), forcing; layout=layout, options=options, io=io)
        flush(io)
        seekstart(io)
        return (
            result=result,
            log=read(io, String),
        )
    end
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
    return Plots
end

result_domain_cpu(result::Chion.EquilibriumResult) = result_domain_cpu(result.domain)

function result_domain_cpu(domain::SM.SnowpackDomain)
    return domain.mass isa Array ? domain : SM.cpu_domain(domain)
end

function summarize_column(result_or_domain, idx::Integer=1)
    domain = result_or_domain isa Chion.EquilibriumResult ? result_domain_cpu(result_or_domain) : result_domain_cpu(result_or_domain)
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

function _layout_grid(layout::Chion.EquilibriumGridLayout, values::AbstractVector{<:Real})
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
    forcing::Chion.EquilibriumForcing;
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
    layout::Chion.EquilibriumGridLayout,
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
