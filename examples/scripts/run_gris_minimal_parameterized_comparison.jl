#!/usr/bin/env julia

"""
Compare five clearly defined Greenland surface-forcing configurations with MAR.

Each configuration uses 15 snow layers and 1, 4, 12, or 24 diurnal substeps.
The model is spun up with the daily MAR climatology before one additional
diagnostic year is collected.

The forcing file uses its native EPSG:3413 polar-stereographic x/y coordinates
(km). `AL2` is a two-sector MAR field; sector one is used as the surface
albedo, consistently with Chion's forcing loader convention for two-level
fields.
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Chion
using Dates: day, month
import CairoMakie
using NCDatasets
using Plots
using Statistics: cor, mean, quantile

const CONFIG = (
    forcing_file=get(ENV, "CHION_FORCING_FILE", "/p/projects/ou/labs/ai/Nils/MAR3.14/MARv3.14.3-10km-daily-ERA5-1940-1980_daily_climatology.nc"),
    output_dir=get(
        ENV,
        "CHION_OUTPUT_DIR",
        joinpath(@__DIR__, "..", "plots", "gris_minimal_parameterized", "step_ntot_sweep"),
    ),
    mask_threshold=parse(Float64, get(ENV, "CHION_MASK_THRESHOLD", "50.0")),
    spinup_years=parse(Int, get(ENV, "CHION_SPINUP_YEARS", "200")),
    backend=Symbol(get(ENV, "CHION_BACKEND", "gpu")),
    # Default to a uniform wind.  Antarctic diagnostic runs can instead pass
    # CHION_USE_MAR_WIND_COMPONENTS=true to use sqrt(U2Z^2 + V2Z^2).
    wind_speed_m_s=5.0,
    use_mar_wind_components=get(ENV, "CHION_USE_MAR_WIND_COMPONENTS", "false") == "true",
    x_name=get(ENV, "CHION_X_NAME", "x"),
    y_name=get(ENV, "CHION_Y_NAME", "y"),
    force_aging_albedo=get(ENV, "CHION_FORCE_AGING_ALBEDO", "false") == "true",
    # Controlled hydrology experiment: use MAR's daily melt total, allocated
    # within a day by the reconstructed net-shortwave cycle. This intentionally
    # replaces, rather than closes, Chion's SEB-derived melt.
    prescribe_melt=get(ENV, "CHION_PRESCRIBE_MELT", "false") == "true",
    antarctica_mode=get(ENV, "CHION_ANTARCTICA_MODE", "false") == "true",
    surface_albedo_name="AL2",
    # Keep the scatter PDFs compact: vectorizing every grid-cell point makes
    # the six monthly budget panels unnecessarily large.  This deterministic
    # subsample is representative; panel statistics still use every valid
    # grid-cell value.
    scatter_max_points=20_000,
)

const COMPARISON_VARIABLES = (
    q_sw_net=("Net shortwave flux (W m⁻²)", "q_sw_net"),
    albedo=("Surface albedo (-)", "albedo"),
    q_lw_down=("Incoming longwave flux (W m⁻²)", "q_lw_down"),
    q_sh=("Sensible heat flux (W m⁻²)", "q_sh"),
    q_lh=("Latent heat flux (W m⁻²)", "q_lh"),
)

const BUDGET_VARIABLES = (
    melt="Melt",
    runoff="Runoff",
    smb_ice="Ice SMB",
    refreezing="Refreezing",
    vapor_mass="Vapor mass",
    surface_balance="SMB",
)

# Ice SMB remains in the mass-budget products, but is excluded from the
# MAR--Chion scatter figure at the requested 2 × 5 layout.
const BUDGET_SCATTER_VARIABLES = (:melt, :runoff, :refreezing, :vapor_mass, :surface_balance)

const CUMULATIVE_BUDGET_VARIABLES = (
    :melt,
    :runoff,
    :smb_ice,
    :refreezing,
    :vapor_mass,
    :sublimation,
)

const ENERGY_VARIABLES = (
    shortwave="Net shortwave",
    longwave="Net longwave",
    sensible_heat="Sensible heat",
    latent_heat="Latent heat",
)

# Fixed limits make a given field directly comparable among every member of
# the diurnal-substep × layer-count sweep. Values are annual mmWE.
const END_STATE_MAP_SPECS = (
    (name=:surface_smb, label="Annual surface SMB [mmWE]", colormap=:RdBu, clims=(-2000.0, 2000.0)),
    (name=:runoff, label="Annual runoff [mmWE]", colormap=:batlow, clims=(0.0, 2500.0)),
    (name=:melt, label="Annual melt [mmWE]", colormap=:batlow, clims=(0.0, 2500.0)),
    (name=:refreezing, label="Annual refreezing [mmWE]", colormap=:batlow, clims=(0.0, 1500.0)),
    (name=:sublimation, label="Annual sublimation [mmWE]", colormap=:batlow, clims=(0.0, 500.0)),
)

const END_STATE_DIFFERENCE_LIMITS = (
    surface_smb=(-500.0, 500.0),
    runoff=(-500.0, 500.0),
    melt=(-500.0, 500.0),
    refreezing=(-500.0, 500.0),
    sublimation=(-50.0, 50.0),
)

# Shared SEMIX surface-flux configuration for every case. Prescribed fields
# bypass their corresponding parameterization, so keeping one model setup
# makes differences among the five forcing families straightforward to audit.
const SURFACE_FLUX_PARAMETERS = (
    ϵ_air=parse(Float64, get(ENV, "CHION_EPSILON_AIR", "0.80")),
    longwave_scheme=Symbol(get(ENV, "CHION_LONGWAVE_SCHEME", "cloud_proxy")),
    # Keep snow emissivity explicit so surface-layer experiments cannot
    # accidentally conflate heat capacity with a longwave-energy change.
    ϵ_snow=0.98,
    D_sh=10.0,
    latent_heat_flux_ratio=1.0,
    semix_sensible_exchange_factor=parse(
        Float64,
        get(ENV, "CHION_SEMIX_SENSIBLE_EXCHANGE_FACTOR", "2.5"),
    ),
    semix_stable_coefficient=parse(
        Float64,
        get(ENV, "CHION_SEMIX_STABLE_COEFFICIENT", "40.0"),
    ),
)

"""Materialize device arrays before CPU-side diagnostic accumulation and plotting."""
to_host(values) = values isa Array ? values : Array(values)

const DIURNAL_SUBSTEPS = Tuple(parse.(Int, split(get(ENV, "CHION_DIURNAL_SUBSTEPS_LIST", "1,4,12,24"), ',')))
const SNOW_LAYERS = 15

const ALL_PRESCRIBED = (
    q_sw_net=true,
    albedo=true,
    q_lw_down=true,
    q_sh=true,
    q_lh=true,
)

const ALL_PARAMETERIZED = (
    q_sw_net=false,
    albedo=false,
    q_lw_down=false,
    q_sh=false,
    q_lh=false,
)

const PRESCRIBED_LONGWAVE = merge(ALL_PARAMETERIZED, (q_lw_down=true,))
const PRESCRIBED_SHORTWAVE = merge(ALL_PARAMETERIZED, (q_sw_net=true, albedo=true))
const PRESCRIBED_ALBEDO = merge(ALL_PARAMETERIZED, (albedo=true,))
const PRESCRIBED_TURBULENCE = merge(ALL_PARAMETERIZED, (q_sh=true, q_lh=true))
const PRESCRIBED_EXCEPT_SENSIBLE_HEAT = merge(ALL_PRESCRIBED, (q_sh=false,))

const CONSTANT_2K_AMPLITUDE = (
    temperature_amplitude_c=2.0,
    temperature_gradient_c_per_km=0.0,
    temperature_reference_height_m=0.0,
    temperature_max_c=2.0,
)

const HEIGHT_DEPENDENT_AMPLITUDE = (
    temperature_amplitude_c=2.0,
    temperature_gradient_c_per_km=1.0,
    temperature_reference_height_m=0.0,
    temperature_max_c=5.0,
)

const CASE_FAMILIES = (
    (
        name=:fully_prescribed,
        title="Fully prescribed",
        prescribed=ALL_PRESCRIBED,
        refreezing_correction=1.0,
        amplitude=HEIGHT_DEPENDENT_AMPLITUDE,
    ),
    (
        name=:parameterized_constant_2k,
        title="Fully parameterized, constant 5 m s⁻¹ wind and 2 K amplitude",
        prescribed=ALL_PARAMETERIZED,
        refreezing_correction=1.0,
        amplitude=CONSTANT_2K_AMPLITUDE,
    ),
    (
        name=:parameterized_constant_2k_refreezing_3,
        title="Fully parameterized, constant 2 K amplitude, refreezing factor 3",
        prescribed=ALL_PARAMETERIZED,
        refreezing_correction=3.0,
        amplitude=CONSTANT_2K_AMPLITUDE,
    ),
    (
        name=:parameterized_prescribed_lwd_refreezing_1,
        title="Prescribed incoming longwave, otherwise parameterized, refreezing factor 1",
        prescribed=PRESCRIBED_LONGWAVE,
        refreezing_correction=1.0,
        amplitude=CONSTANT_2K_AMPLITUDE,
    ),
    (
        name=:parameterized_mar_albedo,
        title="Parameterized SEB with MAR albedo",
        prescribed=PRESCRIBED_ALBEDO,
        refreezing_correction=1.0,
        amplitude=CONSTANT_2K_AMPLITUDE,
    ),
    (
        name=:parameterized_prescribed_sw,
        title="Prescribed MAR net shortwave and albedo",
        prescribed=PRESCRIBED_SHORTWAVE,
        refreezing_correction=1.0,
        amplitude=CONSTANT_2K_AMPLITUDE,
    ),
    (
        name=:parameterized_prescribed_turbulence,
        title="Prescribed MAR sensible and latent heat",
        prescribed=PRESCRIBED_TURBULENCE,
        refreezing_correction=1.0,
        amplitude=CONSTANT_2K_AMPLITUDE,
    ),
    (
        name=:prescribed_except_sensible_heat,
        title="All MAR fields prescribed except sensible heat",
        prescribed=PRESCRIBED_EXCEPT_SENSIBLE_HEAT,
        refreezing_correction=1.0,
        amplitude=CONSTANT_2K_AMPLITUDE,
    ),
    (
        name=:parameterized_height_dependent_amplitude,
        title="Fully parameterized, constant 5 m s⁻¹ wind and height-dependent amplitude",
        prescribed=ALL_PARAMETERIZED,
        refreezing_correction=1.0,
        amplitude=HEIGHT_DEPENDENT_AMPLITUDE,
    ),
)

const CASES = Tuple(
    merge(
        family,
        (
            name=Symbol("$(family.name)_steps_$(substeps)_ntot_$(SNOW_LAYERS)"),
            title="$(family.title) ($(substeps) substep$(substeps == 1 ? "" : "s"))",
            diurnal_substeps=substeps,
            ntot=SNOW_LAYERS,
        ),
        family.amplitude,
    )
    for family in CASE_FAMILIES for substeps in DIURNAL_SUBSTEPS
)

"""Read a MAR field, select its surface sector when needed, and retain grid columns."""
function read_mar_columns(path, name, grid, ntime)
    rows = grid.is .+ (grid.js .- 1) .* length(grid.x)
    return NCDataset(dataset -> read_forcing_columns(dataset, name, rows, ntime, length(grid.y), length(grid.x)), path)
end

"""Read the lowest MAR UVZ level (10 m) for a wind component."""
read_mar_wind_component(path, name, grid, ntime) = read_mar_columns(path, name, grid, ntime)

"""Read a static MAR grid field and retain the selected snowpack columns."""
function read_mar_static_columns(path, name, grid)
    field = NCDataset(path) do dataset
        raw, dim_names = Chion._read_variable_data(dataset, name)
        Chion._as_y_x(raw, dim_names, length(grid.y), length(grid.x), name)
    end
    return [field[grid.js[column], grid.is[column]] for column in eachindex(grid.js)]
end

has_mar_variable(path, name) = NCDataset(path) do dataset
    haskey(dataset, name)
end

"""Regular projected MAR grids store x/y in km, so their cell area is km²."""
function regular_grid_area_km2(grid)
    length(grid.x) > 1 && length(grid.y) > 1 || error("Cannot infer regular-grid cell area from a one-cell axis.")
    area = abs((grid.x[2] - grid.x[1]) * (grid.y[2] - grid.y[1]))
    return fill(area, length(grid.js))
end

"""Return a forcing object with the requested MAR products enabled or disabled."""
function forcing_for_case(forcing, reference, prescribed)
    return SnowpackForcing(
        time_values=forcing.time_values,
        dt_days=forcing.dt_days,
        air_temperature=forcing.air_temperature,
        snowfall_rate=forcing.snowfall_rate,
        rainfall_rate=forcing.rainfall_rate,
        shortwave_down=forcing.shortwave_down,
        wind_speed=forcing.wind_speed,
        q_sw_net=reference.q_sw_net,
        has_q_sw_net=prescribed.q_sw_net,
        q_lw_down=forcing.q_lw_down,
        has_q_lw_down=prescribed.q_lw_down,
        q_sh=forcing.q_sh,
        has_q_sh=prescribed.q_sh,
        q_lh=forcing.q_lh,
        has_q_lh=prescribed.q_lh,
        relative_humidity=forcing.relative_humidity,
        has_relative_humidity=forcing.has_relative_humidity,
        air_pressure=forcing.air_pressure,
        surface_height=forcing.surface_height,
        prescribed_albedo=reference.albedo,
        has_prescribed_albedo=CONFIG.force_aging_albedo ? false : prescribed.albedo,
        prescribed_melt=CONFIG.prescribe_melt ? reference.budgets.melt ./ 86_400.0 : nothing,
        has_prescribed_melt=CONFIG.prescribe_melt,
        latitude_deg=forcing.latitude_deg,
    )
end

"""Replace only the wind forcing while retaining the loaded atmospheric fields."""
function forcing_with_wind_speed(forcing, wind_speed)
    return SnowpackForcing(
        time_values=forcing.time_values,
        dt_days=forcing.dt_days,
        air_temperature=forcing.air_temperature,
        snowfall_rate=forcing.snowfall_rate,
        rainfall_rate=forcing.rainfall_rate,
        shortwave_down=forcing.shortwave_down,
        wind_speed=wind_speed,
        q_sw_net=forcing.q_sw_net,
        has_q_sw_net=forcing.has_q_sw_net,
        q_lw_down=forcing.q_lw_down,
        has_q_lw_down=forcing.has_q_lw_down,
        q_sh=forcing.q_sh,
        has_q_sh=forcing.has_q_sh,
        q_lh=forcing.q_lh,
        has_q_lh=forcing.has_q_lh,
        relative_humidity=forcing.relative_humidity,
        has_relative_humidity=forcing.has_relative_humidity,
        air_pressure=forcing.air_pressure,
        surface_height=forcing.surface_height,
        prescribed_albedo=forcing.prescribed_albedo,
        has_prescribed_albedo=forcing.has_prescribed_albedo,
        latitude_deg=forcing.latitude_deg,
    )
end

"""Build the common BESSI/SEMIX model with only case-specific options varied."""
function model_for_case(grid, case)
    alpha_dry = parse(Float64, get(ENV, "CHION_ALPHA_DRY", "0.81"))
    alpha_wet = parse(Float64, get(ENV, "CHION_ALPHA_WET", "0.70"))
    alpha_ice = parse(Float64, get(ENV, "CHION_ALPHA_ICE", "0.40"))
    max_lwc_albedo = parse(Float64, get(ENV, "CHION_MAX_LWC_ALBEDO", "0.10"))
    parameterized_albedo = Symbol(get(ENV, "CHION_ALBEDO_SCHEME", "dynamic"))
    aging_cold_timescale_days = parse(Float64, get(ENV, "CHION_AGING_COLD_DAYS", "20.0"))
    aging_melting_timescale_days = parse(Float64, get(ENV, "CHION_AGING_MELTING_DAYS", "2.0"))
    constant_diurnal_amplitude = get(ENV, "CHION_DIURNAL_TEMPERATURE_CONSTANT_AMPLITUDE_C", "")
    temperature_amplitude_c = isempty(constant_diurnal_amplitude) ?
                              case.temperature_amplitude_c : parse(Float64, constant_diurnal_amplitude)
    temperature_gradient_c_per_km = isempty(constant_diurnal_amplitude) ?
                                      case.temperature_gradient_c_per_km : 0.0
    temperature_max_c = isempty(constant_diurnal_amplitude) ?
                        case.temperature_max_c : temperature_amplitude_c
    near_surface_layer_max_thicknesses_m = let raw = get(ENV, "CHION_NEAR_SURFACE_LAYER_MAX_THICKNESSES_M", "")
        isempty(raw) ? (0.02, 0.05, 0.10, 0.30) : begin
            values = parse.(Float64, split(raw, ','))
            length(values) == 4 || error("CHION_NEAR_SURFACE_LAYER_MAX_THICKNESSES_M must provide four comma-separated values.")
            (values[1], values[2], values[3], values[4])
        end
    end
    model = BESSIModel(
        grid;
        Ntot=case.ntot,
        albedo=CONFIG.force_aging_albedo ? :aging : (case.prescribed.albedo ? :prescribed : parameterized_albedo),
        seb_scheme=Symbol(get(ENV, "CHION_SEB_SCHEME", "semix")),
        turbulent_flux_scheme=:semix,
        densification=:bessi,
        fresh_snow_density=:constant,
        alpha_dry=CONFIG.antarctica_mode ? 0.82 : alpha_dry,
        alpha_wet=alpha_wet,
        alpha_ice=CONFIG.antarctica_mode ? 0.40 : alpha_ice,
        max_lwc_albedo=max_lwc_albedo,
        aging_cold_timescale_days=aging_cold_timescale_days,
        aging_melting_timescale_days=aging_melting_timescale_days,
        refreezing_correction=case.refreezing_correction,
        near_surface_layer_max_thicknesses_m=near_surface_layer_max_thicknesses_m,
        ice_substrate_layers=parse(Int, get(ENV, "CHION_ICE_SUBSTRATE_LAYERS", "5")),
        ice_substrate_top_thickness_m=parse(Float64, get(ENV, "CHION_ICE_SUBSTRATE_TOP_THICKNESS_M", "0.05")),
        diurnal_shortwave_substeps=true,
        diurnal_shortwave_max_substeps=case.diurnal_substeps,
        # The reconstructed solar cycle can be assessed independently of an
        # assumed air-temperature cycle.  Keep the historical default, while
        # allowing controlled shortwave-only convergence experiments.
        diurnal_temperature_cycle=get(ENV, "CHION_DIURNAL_TEMPERATURE_CYCLE", "true") == "true",
        diurnal_temperature_amplitude_c=temperature_amplitude_c,
        diurnal_temperature_amplitude_gradient_c_per_km=temperature_gradient_c_per_km,
        diurnal_temperature_amplitude_reference_height_m=case.temperature_reference_height_m,
        diurnal_temperature_amplitude_max_c=temperature_max_c,
        SURFACE_FLUX_PARAMETERS...,
    )
    @info "BESSI case runtime options" diurnal_substeps=case.diurnal_substeps diurnal_temperature_cycle=model.diurnal_temperature_cycle near_surface_layer_max_thicknesses_m=model.parameters.near_surface_layer_max_thicknesses_m ice_substrate_layers=model.parameters.ice_substrate_layers
    return model
end

function comparison_buffers(ncol)
    make_buffer() = Dict(name => zeros(Float64, 12, ncol) for name in keys(COMPARISON_VARIABLES))
    return (truth=make_buffer(), configured=make_buffer(), month_days=zeros(Float64, 12))
end

function budget_buffers(ncol)
    make_buffer() = Dict(name => zeros(Float64, 12, ncol) for name in keys(BUDGET_VARIABLES))
    return (truth=make_buffer(), configured=make_buffer(), month_days=zeros(Float64, 12))
end

function energy_buffers(ncol)
    make_buffer() = Dict(name => zeros(Float64, 12, ncol) for name in keys(ENERGY_VARIABLES))
    return (truth=make_buffer(), configured=make_buffer())
end

"""Buffers for diagnosing the turbulent vapour-pressure closure itself."""
function latent_diagnostic_buffers(ncol)
    return (
        vapor_pressure_deficit=zeros(Float64, 12, ncol),
        # dQ_LH / d(e_a - e_si), evaluated at the diagnosed surface state.
        # This makes an effective-vapour-pressure fit explicit and avoids
        # inferring a transfer coefficient by dividing by a near-zero deficit.
        vapor_conductance=zeros(Float64, 12, ncol),
        chion_latent_heat=zeros(Float64, 12, ncol),
        mar_latent_heat=zeros(Float64, 12, ncol),
        month_days=zeros(Float64, 12),
    )
end

function accumulate_latent_diagnostics!(buffers, time, dt_days, vapor_pressure_deficit, vapor_conductance, chion_latent_heat, mar_latent_heat)
    month_index = month(time)
    @views buffers.vapor_pressure_deficit[month_index, :] .+= vapor_pressure_deficit .* dt_days
    @views buffers.vapor_conductance[month_index, :] .+= vapor_conductance .* dt_days
    @views buffers.chion_latent_heat[month_index, :] .+= chion_latent_heat .* dt_days
    @views buffers.mar_latent_heat[month_index, :] .+= mar_latent_heat .* dt_days
    buffers.month_days[month_index] += dt_days
    return nothing
end

function finalize_latent_diagnostics(buffers)
    for month_index in 1:12
        @views buffers.vapor_pressure_deficit[month_index, :] ./= buffers.month_days[month_index]
        @views buffers.vapor_conductance[month_index, :] ./= buffers.month_days[month_index]
        @views buffers.chion_latent_heat[month_index, :] ./= buffers.month_days[month_index]
        @views buffers.mar_latent_heat[month_index, :] ./= buffers.month_days[month_index]
    end
    annual_mean(values) = vec(sum(values .* reshape(buffers.month_days, :, 1); dims=1) ./ sum(buffers.month_days))
    return (; buffers...,
        annual_vapor_pressure_deficit=annual_mean(buffers.vapor_pressure_deficit),
        annual_chion_latent_heat=annual_mean(buffers.chion_latent_heat),
        annual_mar_latent_heat=annual_mean(buffers.mar_latent_heat),
    )
end

"""Write the monthly fields required to fit an effective air-vapour pressure."""
function write_latent_flux_fit_input(path, diagnostics, grid, surface_height, area_km2)
    NCDataset(path, "c") do dataset
        ncolumn = length(grid.js)
        defDim(dataset, "month", 12)
        defDim(dataset, "column", ncolumn)
        defVar(dataset, "month", Int32, ("month",))[:] = Int32.(1:12)
        defVar(dataset, "month_days", Float64, ("month",); attrib=Dict("units" => "days"))[:] = diagnostics.month_days
        defVar(dataset, "vapor_pressure_deficit", Float64, ("month", "column"); attrib=Dict("units" => "Pa", "long_name" => "e_a minus e_si at Chion surface temperature"))[:, :] = diagnostics.vapor_pressure_deficit
        defVar(dataset, "vapor_conductance", Float64, ("month", "column"); attrib=Dict("units" => "W m-2 Pa-1", "long_name" => "bulk latent-flux sensitivity to e_a minus e_si"))[:, :] = diagnostics.vapor_conductance
        defVar(dataset, "chion_latent_heat", Float64, ("month", "column"); attrib=Dict("units" => "W m-2"))[:, :] = diagnostics.chion_latent_heat
        defVar(dataset, "mar_latent_heat", Float64, ("month", "column"); attrib=Dict("units" => "W m-2"))[:, :] = diagnostics.mar_latent_heat
        defVar(dataset, "surface_height", Float64, ("column",); attrib=Dict("units" => "m"))[:] = surface_height
        defVar(dataset, "area", Float64, ("column",); attrib=Dict("units" => "km2"))[:] = area_km2
        defVar(dataset, "x", Float64, ("column",); attrib=Dict("units" => "km"))[:] = grid.x[grid.is]
        defVar(dataset, "y", Float64, ("column",); attrib=Dict("units" => "km"))[:] = grid.y[grid.js]
        dataset.attrib["title"] = "Input for effective air-vapour-pressure fit; diagnostic-only, not a model forcing"
        dataset.attrib["formula"] = "Q_LH_corrected = Q_LH_chion - vapor_conductance * delta_e_a"
    end
    return nothing
end

"""Accumulate a daily set of values into calendar-month means."""
function accumulate_comparison!(buffers, time, dt_days, truth, configured)
    month_index = month(time)
    for name in keys(COMPARISON_VARIABLES)
        @views buffers.truth[name][month_index, :] .+= truth[name] .* dt_days
        @views buffers.configured[name][month_index, :] .+= configured[name] .* dt_days
    end
    buffers.month_days[month_index] += dt_days
    return nothing
end

function finalize_comparison(buffers)
    for name in keys(COMPARISON_VARIABLES), month_index in 1:12
        @views buffers.truth[name][month_index, :] ./= buffers.month_days[month_index]
        @views buffers.configured[name][month_index, :] ./= buffers.month_days[month_index]
    end
    yearly_truth = Dict(
        name => vec(sum(buffers.truth[name] .* reshape(buffers.month_days, :, 1); dims=1) ./ sum(buffers.month_days))
        for name in keys(COMPARISON_VARIABLES)
    )
    yearly_configured = Dict(
        name => vec(sum(buffers.configured[name] .* reshape(buffers.month_days, :, 1); dims=1) ./ sum(buffers.month_days))
        for name in keys(COMPARISON_VARIABLES)
    )
    return (; buffers..., yearly_truth, yearly_configured)
end

"""Accumulate interval mass changes into calendar-month totals."""
function accumulate_budget!(buffers, time, dt_days, truth, configured)
    month_index = month(time)
    for name in keys(BUDGET_VARIABLES)
        @views buffers.truth[name][month_index, :] .+= truth[name]
        @views buffers.configured[name][month_index, :] .+= configured[name]
    end
    buffers.month_days[month_index] += dt_days
    return nothing
end

function finalize_budget(buffers)
    yearly_truth = Dict(name => vec(sum(buffers.truth[name]; dims=1)) for name in keys(BUDGET_VARIABLES))
    yearly_configured = Dict(name => vec(sum(buffers.configured[name]; dims=1)) for name in keys(BUDGET_VARIABLES))
    return (; buffers..., yearly_truth, yearly_configured)
end

"""Accumulate interval energy changes (J m⁻²) into calendar-month totals."""
function accumulate_energy!(buffers, time, truth, configured)
    month_index = month(time)
    for name in keys(ENERGY_VARIABLES)
        @views buffers.truth[name][month_index, :] .+= truth[name]
        @views buffers.configured[name][month_index, :] .+= configured[name]
    end
    return nothing
end

function selected_scatter_points(x, y, max_points)
    valid = findall(isfinite.(x) .& isfinite.(y))
    isempty(valid) && return Float64[], Float64[]
    indices = length(valid) <= max_points ? valid : valid[round.(Int, range(1, length(valid); length=max_points))]
    return x[indices], y[indices]
end

function scatter_panel(x, y, label, period, max_points; unit=nothing, panel_label=nothing)
    valid = isfinite.(x) .& isfinite.(y)
    xmetrics, ymetrics = x[valid], y[valid]
    xplot, yplot = selected_scatter_points(x, y, max_points)
    isempty(xplot) && return plot(title="$label — $period (no valid data)")
    # Plot a representative subset, but calculate skill from every valid
    # grid-cell value so file-size optimization cannot change the diagnostics.
    rmse = sqrt(mean((ymetrics .- xmetrics) .^ 2))
    correlation = length(xmetrics) > 1 ? cor(xmetrics, ymetrics) : NaN
    r_squared = correlation^2
    lower = min(minimum(xplot), minimum(yplot))
    upper = max(maximum(xplot), maximum(yplot))
    pad = max((upper - lower) * 0.04, eps(Float64))
    limits = (lower - pad, upper + pad)
    panel = scatter(
        xplot,
        yplot;
        markersize=1.0,
        markerstrokewidth=0,
        markeralpha=0.16,
        label="Chion",
        xlabel="$period $label (MAR)$(isnothing(unit) ? "" : " [$unit]")",
        ylabel="$period $label (Chion)$(isnothing(unit) ? "" : " [$unit]")",
        title=isnothing(panel_label) ? "" : panel_label,
        titlelocation=:left,
        titlefont=Plots.font(11),
        legend=:bottomright,
        aspect_ratio=:equal,
        xlims=limits,
        ylims=limits,
    )
    plot!(panel, [limits[1], limits[2]], [limits[1], limits[2]]; color=:black, linewidth=1.2, linestyle=:dash, label="1:1")
    stats = length(xplot) > 1 ? "R²=$(round(r_squared; digits=3))\nr=$(round(correlation; digits=3))\nRMSE=$(round(rmse; digits=3))" : "One annual total"
    annotate!(panel, limits[1] + 0.04 * (limits[2] - limits[1]), limits[2] - 0.04 * (limits[2] - limits[1]),
        Plots.text(stats, 8, :left, :top))
    return panel
end

function comparison_figure(comparison, case)
    panels = Any[]
    for (name, (label, _)) in pairs(COMPARISON_VARIABLES)
        case.prescribed[name] && continue
        push!(panels, scatter_panel(vec(comparison.truth[name]), vec(comparison.configured[name]), label, "monthly", CONFIG.scatter_max_points))
        push!(panels, scatter_panel(comparison.yearly_truth[name], comparison.yearly_configured[name], label, "yearly", CONFIG.scatter_max_points))
    end
    isempty(panels) && return plot(title="$(case.title): all variables are prescribed")
    return plot(
        panels...;
        layout=(length(panels) ÷ 2, 2),
        size=(1500, 440 * (length(panels) ÷ 2)),
        left_margin=8Plots.mm,
        bottom_margin=6Plots.mm,
        top_margin=8Plots.mm,
    )
end

"""Add one rasterized budget scatter layer to a CairoMakie axis."""
function budget_scatter_panel!(position, x, y, label, period; unit, panel_label)
    valid = isfinite.(x) .& isfinite.(y)
    xmetrics, ymetrics = x[valid], y[valid]
    xplot, yplot = selected_scatter_points(x, y, CONFIG.scatter_max_points)
    axis = CairoMakie.Axis(
        position;
        xlabel="$period $label (MAR) [$unit]",
        ylabel="$period $label (Chion) [$unit]",
        title="",
        aspect=CairoMakie.DataAspect(),
    )
    isempty(xplot) && return axis

    rmse = sqrt(mean((ymetrics .- xmetrics) .^ 2))
    correlation = length(xmetrics) > 1 ? cor(xmetrics, ymetrics) : NaN
    lower = min(minimum(xplot), minimum(yplot))
    upper = max(maximum(xplot), maximum(yplot))
    pad = max((upper - lower) * 0.04, eps(Float64))
    limits = (lower - pad, upper + pad)
    CairoMakie.xlims!(axis, limits...)
    CairoMakie.ylims!(axis, limits...)
    # Only the dense marker series is rasterized in the PDF. Axes, labels,
    # annotations, and the 1:1 reference line remain vector graphics.
    CairoMakie.scatter!(
        axis,
        xplot,
        yplot;
        color=(:dodgerblue, 0.16),
        markersize=3.0,
        rasterize=2,
    )
    CairoMakie.lines!(axis, collect(limits), collect(limits); color=:black, linewidth=1.5, linestyle=:dash)
    stats = "R²=$(round(correlation^2; digits=3))\nr=$(round(correlation; digits=3))\nRMSE=$(round(rmse; digits=3))"
    CairoMakie.text!(
        axis,
        stats;
        position=(limits[1] + 0.04 * (limits[2] - limits[1]), limits[2] - 0.04 * (limits[2] - limits[1])),
        align=(:left, :top),
        fontsize=14,
    )
    return axis
end

"""Budget scatter figure with rasterized point layers and vector annotations."""
function budget_figure(comparison, area_km2, case)
    # 1 mmWE over 1 km² is 10⁻³ Mt. The monthly panel has 12 × ncell points.
    monthly_cell_total(values) = values .* reshape(area_km2, 1, :) .* 1e-3
    # Near-square cells preserve the 1:1 scatter aspect without unused
    # horizontal space; minimal layout gaps keep the 2 × 5 grid compact.
    figure = CairoMakie.Figure(size=(3000, 1200), fontsize=18, figure_padding=(2, 2, 2, 2))
    CairoMakie.colgap!(figure.layout, 2)
    CairoMakie.rowgap!(figure.layout, 4)
    for (index, name) in enumerate(BUDGET_SCATTER_VARIABLES)
        label = BUDGET_VARIABLES[name]
        monthly_panel_label = "($(Char('a' + index - 1)))"
        annual_panel_label = "($(Char('a' + length(BUDGET_SCATTER_VARIABLES) + index - 1)))"
        monthly_mar = monthly_cell_total(comparison.truth[name])
        monthly_chion = monthly_cell_total(comparison.configured[name])
        yearly_mar = vec(sum(monthly_mar; dims=1))
        yearly_chion = vec(sum(monthly_chion; dims=1))
        budget_scatter_panel!(
            figure[1, index], vec(monthly_mar), vec(monthly_chion), label, "Monthly";
            unit="Mt/month", panel_label=monthly_panel_label,
        )
        CairoMakie.Label(
            figure[1, index, CairoMakie.TopLeft()], monthly_panel_label;
            fontsize=20, halign=:left, valign=:bottom, padding=(0, 0, 4, 0),
        )
        budget_scatter_panel!(
            figure[2, index], yearly_mar, yearly_chion, label, "Annual";
            unit="Mt/yr", panel_label=annual_panel_label,
        )
        CairoMakie.Label(
            figure[2, index, CairoMakie.TopLeft()], annual_panel_label;
            fontsize=20, halign=:left, valign=:bottom, padding=(0, 0, 4, 0),
        )
    end
    return figure
end

"""Return annual per-grid-cell budget-scatter skill metrics in Mt/yr."""
function yearly_budget_scatter_metrics(comparison, area_km2, case)
    # 1 mmWE over 1 km² is 10⁻³ Mt.
    yearly_cell_total(values) = vec(sum(values .* reshape(area_km2, 1, :); dims=1)) .* 1e-3
    rows = NamedTuple[]
    for (name, _) in pairs(BUDGET_VARIABLES)
        mar = yearly_cell_total(comparison.truth[name])
        chion = yearly_cell_total(comparison.configured[name])
        selected = isfinite.(mar) .& isfinite.(chion)
        x, y = mar[selected], chion[selected]
        isempty(x) && continue
        push!(rows, (
            variable=name,
            diurnal_steps=case.diurnal_substeps,
            ntot=case.ntot,
            rmse=sqrt(mean((y .- x) .^ 2)),
            correlation=length(x) > 1 ? cor(x, y) : NaN,
        ))
    end
    return rows
end

function write_yearly_budget_scatter_metrics(path, rows)
    open(path, "w") do io
        println(io, "variable,diurnal_steps,ntot,rmse_mt_per_year,correlation")
        for row in rows
            println(io, "$(row.variable),$(row.diurnal_steps),$(row.ntot),$(row.rmse),$(row.correlation)")
        end
    end
    return nothing
end

"""Place an annual-total annotation above, rather than over, plotted series."""
function annotate_integrated_total!(panel, values, annotation)
    lower, upper = extrema(values)
    span = max(upper - lower, max(abs(lower), abs(upper), 1.0) * 0.1)
    padding = 0.18 * span
    ylims!(panel, (lower - 0.03 * span, upper + padding))
    annotate!(panel, 6.5, upper + 0.92 * padding, Plots.text(annotation, 9, :center, :top))
    return panel
end

"""Plot monthly area-integrated budgets in Gt/month, with annual totals in Gt/yr."""
function monthly_integrated_budget_figure(comparison, area_km2, case)
    # 1 mmWE over 1 km² is 10⁶ kg = 10⁻⁶ Gt.
    integrate_monthly(values) = vec(sum(values .* reshape(area_km2, 1, :); dims=2)) .* 1e-6
    panels = Any[]
    for (name, label) in pairs(BUDGET_VARIABLES)
        panel_index = length(panels) + 1
        panel_label = "($(Char('a' + panel_index - 1)))"
        mar = integrate_monthly(comparison.truth[name])
        chion = integrate_monthly(comparison.configured[name])
        mar_yearly = sum(mar)
        chion_yearly = sum(chion)
        labels = ["MAR" "Chion"]
        annotation = "MAR=$(round(mar_yearly; digits=1)) Gt/yr\nChion=$(round(chion_yearly; digits=1)) Gt/yr"
        if name == :smb_ice
            delta_firn = integrate_monthly(comparison.delta_firn)
            delta_firn_yearly = sum(delta_firn)
            chion_surface_smb = chion_yearly + delta_firn_yearly
            annotation = "MAR=$(round(mar_yearly; digits=1)) Gt/yr\nChion=$(round(chion_yearly; digits=1)) Gt/yr\nΔM firn=$(round(delta_firn_yearly; digits=1)) Gt/yr\nChion + ΔM=$(round(chion_surface_smb; digits=1)) Gt/yr"
            panel = plot(
                1:12,
                mar;
                label="MAR surface SMB",
                marker=:circle,
                markerstrokewidth=0,
                linewidth=2.8,
                xlabel="Month",
                ylabel="$label [Gt/month]",
                title=panel_label,
                titlelocation=:left,
                titlefont=Plots.font(11),
                xticks=1:12,
                framestyle=:box,
            )
            plot!(
                panel,
                1:12,
                chion;
                label="Chion ice SMB",
                marker=:square,
                markerstrokewidth=0,
                linewidth=2.8,
            )
            plot!(
                panel,
                1:12,
                chion .+ delta_firn;
                label="Chion ice SMB + ΔM firn",
                marker=:diamond,
                markerstrokewidth=0,
                linewidth=2.8,
            )
            annotate_integrated_total!(panel, vcat(mar, chion, chion .+ delta_firn), annotation)
            push!(panels, panel)
            continue
        end
        push!(
            panels,
            plot(
                1:12,
                hcat(mar, chion);
                label=labels,
                marker=[:circle :square],
                markerstrokewidth=0,
                linewidth=2.8,
                xlabel="Month",
                ylabel="$label [Gt/month]",
                title=panel_label,
                titlelocation=:left,
                titlefont=Plots.font(11),
                xticks=1:12,
                framestyle=:box,
            ),
        )
        annotate_integrated_total!(panels[end], vcat(mar, chion), annotation)
    end
    return plot(
        panels...;
        layout=(2, 3),
        size=(1800, 1150),
        left_margin=8Plots.mm,
        bottom_margin=6Plots.mm,
        top_margin=8Plots.mm,
    )
end

"""Plot final-year monthly integrated surface-energy components in PJ."""
function monthly_integrated_energy_figure(comparison, area_km2, case)
    # 1 J m⁻² over 1 km² is 10⁻⁹ PJ.
    integrate_monthly(values) = vec(sum(values .* reshape(area_km2, 1, :); dims=2)) .* 1e-9
    panels = Any[]
    for (name, label) in pairs(ENERGY_VARIABLES)
        panel_index = length(panels) + 1
        panel_label = "($(Char('a' + panel_index - 1)))"
        mar = integrate_monthly(comparison.truth[name])
        chion = integrate_monthly(comparison.configured[name])
        annotation = "MAR=$(round(sum(mar); digits=1)) PJ\nChion=$(round(sum(chion); digits=1)) PJ"
        panel = plot(
            1:12,
            hcat(mar, chion);
            label=["MAR" "Chion"],
            marker=[:circle :square],
            markerstrokewidth=0,
            linewidth=2.8,
            xlabel="Month",
            ylabel="$label [PJ]",
            title=panel_label,
            titlelocation=:left,
            titlefont=Plots.font(11),
            xticks=1:12,
            framestyle=:box,
        )
        annotate_integrated_total!(panel, vcat(mar, chion), annotation)
        push!(panels, panel)
    end
    return plot(
        panels...;
        layout=(2, 2),
        size=(1400, 1050),
        left_margin=8Plots.mm,
        bottom_margin=6Plots.mm,
        top_margin=8Plots.mm,
    )
end

"""Return NaN-separated coastline segments from the MAR land-mask cell edges."""
function coastline_segments(grid; threshold)
    mask = grid.mask .>= threshold
    ny, nx = size(mask)
    dx = length(grid.x) > 1 ? abs(grid.x[2] - grid.x[1]) : 1.0
    dy = length(grid.y) > 1 ? abs(grid.y[2] - grid.y[1]) : 1.0
    xs, ys = Float64[], Float64[]
    function append_edge!(x1, y1, x2, y2)
        append!(xs, (x1, x2, NaN))
        append!(ys, (y1, y2, NaN))
        return nothing
    end
    for j in 1:ny, i in 1:nx
        mask[j, i] || continue
        x, y = grid.x[i], grid.y[j]
        left = x - dx / 2
        right = x + dx / 2
        bottom = y - dy / 2
        top = y + dy / 2
        (i == 1 || !mask[j, i - 1]) && append_edge!(left, bottom, left, top)
        (i == nx || !mask[j, i + 1]) && append_edge!(right, bottom, right, top)
        (j == 1 || !mask[j - 1, i]) && append_edge!(left, bottom, right, bottom)
        (j == ny || !mask[j + 1, i]) && append_edge!(left, top, right, top)
    end
    return xs, ys
end

function map_panel(
    grid,
    values,
    title;
    coastline_level,
    colorbar_label=title,
    clims=nothing,
    clip_quantiles=nothing,
    colormap=nothing,
    show_title=true,
    panel_label=nothing,
)
    field = Chion.scatter_to_grid(Float64.(to_host(values)), grid.js, grid.is, size(grid.mask))
    finite_values = field[isfinite.(field)]
    isempty(finite_values) && return plot(
        title="$title (no finite values)",
        aspect_ratio=:equal,
        axis=false,
        grid=false,
        framestyle=:none,
    )
    nonnegative = all(>=(zero(eltype(finite_values))), finite_values)
    resolved_colormap = isnothing(colormap) ? (nonnegative ? :batlow : :RdBu) : colormap
    limits = if !isnothing(clims)
        clims
    elseif !isnothing(clip_quantiles)
        lower_quantile, upper_quantile = clip_quantiles
        if nonnegative
            (quantile(finite_values, lower_quantile), quantile(finite_values, upper_quantile))
        else
            extent = quantile(abs.(finite_values), upper_quantile)
            (-extent, extent)
        end
    elseif nonnegative
        extrema(finite_values)
    else
        extent = maximum(abs, finite_values)
        (-extent, extent)
    end
    if limits[1] == limits[2]
        padding = max(abs(limits[1]) * 0.05, one(limits[1]))
        limits = (limits[1] - padding, limits[2] + padding)
    end
    panel = heatmap(
        grid.x,
        grid.y,
        field;
        title=show_title ? title : "",
        aspect_ratio=:equal,
        c=resolved_colormap,
        clims=limits,
        colorbar=true,
        colorbar_title=colorbar_label,
        colorbar_tickfontsize=8,
        colorbar_titlefontsize=11,
        axis=false,
        grid=false,
        framestyle=:none,
        margin=0Plots.mm,
    )
    coastline_x, coastline_y = coastline_segments(grid; threshold=coastline_level)
    plot!(panel, coastline_x, coastline_y; color=:black, linewidth=0.35, label=false, primary=false)
    if !isnothing(panel_label)
        dx = length(grid.x) > 1 ? abs(grid.x[2] - grid.x[1]) : 1.0
        dy = length(grid.y) > 1 ? abs(grid.y[2] - grid.y[1]) : 1.0
        x_min, x_max = extrema(grid.x)
        y_min, y_max = extrema(grid.y)
        map_left = x_min - dx / 2
        map_top = y_max + dy / 2
        plot!(panel; xlims=(map_left - 2dx, x_max + dx / 2), ylims=(y_min - dy / 2, map_top + 1.5dy))
        annotate!(panel, map_left, map_top, Plots.text(panel_label, 13, :right, :bottom))
    end
    return panel
end

function end_state_figure(grid, state, annual, case)
    panels = [
        map_panel(
            grid,
            getfield(annual, spec.name),
            spec.label;
            coastline_level=CONFIG.mask_threshold,
            colorbar_label=spec.label,
            colormap=spec.colormap,
            clims=spec.clims,
            show_title=false,
            panel_label="($(Char('a' + index - 1)))",
        )
        for (index, spec) in enumerate(END_STATE_MAP_SPECS)
    ]
    return plot(
        panels...;
        # Five annual surface-budget maps; ice SMB is intentionally omitted.
        layout=(2, 3),
        size=(1800, 1050),
        left_margin=4Plots.mm,
        right_margin=1Plots.mm,
        bottom_margin=0Plots.mm,
        top_margin=2Plots.mm,
    )
end

"""Annual MAR budget fields comparable with Chion's end-state diagnostics."""
function annual_mar_end_state(reference, dt_days)
    annual_sum(values) = vec(sum(values .* reshape(dt_days, 1, :); dims=2))
    return (
        smb_ice=annual_sum(reference.budgets.smb_ice),
        surface_smb=annual_sum(reference.budgets.surface_balance),
        runoff=annual_sum(reference.budgets.runoff),
        melt=annual_sum(reference.budgets.melt),
        refreezing=annual_sum(reference.budgets.refreezing),
        # MAR SU is positive for sublimation; the reference vapor field uses
        # the opposite sign to match Chion's deposition-positive convention.
        sublimation=annual_sum(-reference.budgets.vapor_mass),
    )
end

"""Plot Chion minus MAR annual budget maps using Crameri's `bam` colormap."""
function end_state_difference_figure(grid, mar_annual, annual, case)
    panels = [
        map_panel(
            grid,
            getfield(annual, spec.name) .- getfield(mar_annual, spec.name),
            "Chion − MAR $(spec.label)";
            coastline_level=CONFIG.mask_threshold,
            colorbar_label="$(replace(spec.label, " [mmWE]" => "")) Difference [mmWE]",
            colormap=:bam,
            clims=getfield(END_STATE_DIFFERENCE_LIMITS, spec.name),
            show_title=false,
            panel_label="($(Char('a' + index - 1)))",
        )
        for (index, spec) in enumerate(END_STATE_MAP_SPECS)
    ]
    return plot(
        panels...;
        layout=(2, 3),
        size=(1800, 1050),
        left_margin=4Plots.mm,
        right_margin=1Plots.mm,
        bottom_margin=0Plots.mm,
        top_margin=2Plots.mm,
    )
end

"""Save map figures with narrow GR colourbars, without changing other figures."""
function save_map_figure(path, figure)
    original_width = Plots.gr_cbar_width[]
    Plots.gr_cbar_width[] = 0.009
    try
        savefig(figure, path)
    finally
        Plots.gr_cbar_width[] = original_width
    end
    return nothing
end

"""Map the thermodynamic driver of Chion's latent flux and its MAR mismatch."""
function latent_flux_diagnostic_figure(grid, diagnostics, case)
    annual_difference = diagnostics.annual_chion_latent_heat .- diagnostics.annual_mar_latent_heat
    winter_difference = vec(mean((diagnostics.chion_latent_heat .- diagnostics.mar_latent_heat)[[1, 2, 3, 12], :]; dims=1))
    summer_difference = vec(mean((diagnostics.chion_latent_heat .- diagnostics.mar_latent_heat)[6:8, :]; dims=1))
    panels = (
        map_panel(grid, diagnostics.annual_vapor_pressure_deficit,
            "Annual mean eₐ − eₛᵢ(Tₛ) [Pa]";
            coastline_level=CONFIG.mask_threshold, colormap=:bam, clims=(-150.0, 150.0), show_title=false, panel_label="(a)"),
        map_panel(grid, annual_difference,
            "Annual Chion − MAR latent heat [W m⁻²]";
            coastline_level=CONFIG.mask_threshold, colormap=:bam, clims=(-30.0, 30.0), show_title=false, panel_label="(b)"),
        map_panel(grid, winter_difference,
            "DJFM Chion − MAR latent heat [W m⁻²]";
            coastline_level=CONFIG.mask_threshold, colormap=:bam, clims=(-30.0, 30.0), show_title=false, panel_label="(c)"),
        map_panel(grid, summer_difference,
            "JJAS Chion − MAR latent heat [W m⁻²]";
            coastline_level=CONFIG.mask_threshold, colormap=:bam, clims=(-30.0, 30.0), show_title=false, panel_label="(d)"),
    )
    return plot(panels...; layout=(2, 2), size=(1400, 1100), title="$(case.title): latent-flux closure diagnosis")
end

function latent_flux_elevation_summary(diagnostics, surface_height, area_km2)
    edges = collect(0.0:500.0:3500.0)
    difference = diagnostics.annual_chion_latent_heat .- diagnostics.annual_mar_latent_heat
    rows = NamedTuple[]
    for (lower, upper) in zip(edges[1:end-1], edges[2:end])
        selected = (surface_height .>= lower) .& (surface_height .< upper)
        any(selected) || continue
        weights = area_km2[selected]
        weighted_mean(values) = sum(values[selected] .* weights) / sum(weights)
        push!(rows, (
            lower=lower,
            upper=upper,
            vapor_pressure_deficit_pa=weighted_mean(diagnostics.annual_vapor_pressure_deficit),
            chion_latent_heat_w_m2=weighted_mean(diagnostics.annual_chion_latent_heat),
            mar_latent_heat_w_m2=weighted_mean(diagnostics.annual_mar_latent_heat),
            difference_w_m2=weighted_mean(difference),
        ))
    end
    return rows
end

function write_latent_flux_elevation_summary(path, rows)
    open(path, "w") do io
        println(io, "elevation_lower_m,elevation_upper_m,vapor_pressure_deficit_pa,chion_latent_heat_w_m2,mar_latent_heat_w_m2,chion_minus_mar_w_m2")
        for row in rows
            println(io, "$(row.lower),$(row.upper),$(row.vapor_pressure_deficit_pa),$(row.chion_latent_heat_w_m2),$(row.mar_latent_heat_w_m2),$(row.difference_w_m2)")
        end
    end
    return nothing
end

"""Summarize annual MAR--Chion refreezing by surface-elevation band in Gt."""
function refreezing_elevation_summary(mar_refreezing, chion_refreezing, mar_melt, chion_melt, surface_height, area_km2)
    edges = collect(0.0:500.0:3500.0)
    rows = NamedTuple[]
    for (lower, upper) in zip(edges[1:end-1], edges[2:end])
        selected = (surface_height .>= lower) .& (surface_height .< upper)
        isempty(findall(selected)) && continue
        mar = sum(mar_refreezing[selected] .* area_km2[selected]) * 1e-6
        chion = sum(chion_refreezing[selected] .* area_km2[selected]) * 1e-6
        push!(rows, (
            lower=lower,
            upper=upper,
            mar=mar,
            chion=chion,
            deficit=mar - chion,
            mar_melt=sum(mar_melt[selected] .* area_km2[selected]) * 1e-6,
            chion_melt=sum(chion_melt[selected] .* area_km2[selected]) * 1e-6,
        ))
    end
    return rows
end

function write_refreezing_elevation_summary(path, rows)
    open(path, "w") do io
        println(io, "elevation_lower_m,elevation_upper_m,mar_refreezing_gt,chion_refreezing_gt,mar_minus_chion_gt,mar_melt_gt,chion_melt_gt")
        for row in rows
            println(io, "$(row.lower),$(row.upper),$(row.mar),$(row.chion),$(row.deficit),$(row.mar_melt),$(row.chion_melt)")
        end
    end
    return nothing
end

function refreezing_difference_figure(grid, mar_refreezing, chion_refreezing, case)
    shared_upper = max(maximum(mar_refreezing), maximum(chion_refreezing))
    shared_limits = (0.0, shared_upper)
    difference = chion_refreezing .- mar_refreezing
    difference_extent = maximum(abs, difference)
    panels = (
        map_panel(grid, mar_refreezing, "MAR annual refreezing";
            coastline_level=CONFIG.mask_threshold, colorbar_label="Annual refreezing [mmWE]", clims=shared_limits),
        map_panel(grid, chion_refreezing, "Chion annual refreezing";
            coastline_level=CONFIG.mask_threshold, colorbar_label="Annual refreezing [mmWE]", clims=shared_limits),
        map_panel(grid, difference, "Chion − MAR annual refreezing";
            coastline_level=CONFIG.mask_threshold, colorbar_label="Chion − MAR annual refreezing [mmWE]", clims=(-difference_extent, difference_extent)),
    )
    return plot(panels...; layout=(1, 3), size=(1800, 600), title="$(case.title): refreezing diagnosis")
end

function annual_budget_start(state)
    return NamedTuple(name => copy(to_host(getfield(state, name))) for name in CUMULATIVE_BUDGET_VARIABLES)
end

function annual_budget_change(state, start)
    return NamedTuple(name => to_host(getfield(state, name)) .- getfield(start, name) for name in keys(start))
end

"""Total solid plus liquid firn mass per column (mmWE)."""
function firn_mass(state)
    # `TransposedLayerMatrix` is a GPU storage wrapper whose elementwise
    # broadcast falls back to scalar indexing. Transfer each logical field
    # first, then form this CPU-side diagnostic.
    return vec(sum(to_host(state.mass) .+ to_host(state.mass_w); dims=1))
end

"""Diagnose the firn state immediately before the 1 May forcing interval."""
function pre_melt_firn_state(state)
    N = to_host(state.N)
    mass = to_host(state.mass)
    mass_w = to_host(state.mass_w)
    density = to_host(state.density)
    temperature = to_host(state.temperature)
    ncol = length(N)
    cold_content = zeros(Float64, ncol)
    pore_capacity = zeros(Float64, ncol)
    liquid_water = zeros(Float64, ncol)
    depth_to_dense_layer = fill(NaN, ncol)
    c = state.c
    dense_threshold = 810.0
    for column in 1:ncol
        depth = 0.0
        for layer in 1:N[column]
            layer_mass = mass[layer, column]
            layer_density = density[layer, column]
            if layer_mass <= 0 || layer_density <= 0
                continue
            end
            cold_content[column] += max(c.T0 - temperature[layer, column], 0.0) * c.ci * layer_mass / c.Lm
            pore_capacity[column] += max(layer_mass / layer_density - layer_mass / c.rho_i, 0.0) * c.rho_w
            liquid_water[column] += mass_w[layer, column]
            if isnan(depth_to_dense_layer[column]) && layer_density >= dense_threshold
                depth_to_dense_layer[column] = depth
            end
            depth += layer_mass / layer_density
        end
    end
    return (
        cold_content=cold_content,
        pore_capacity=pore_capacity,
        liquid_water=liquid_water,
        depth_to_dense_layer=depth_to_dense_layer,
    )
end

function firn_state_elevation_summary(firn, refreezing, surface_height, area_km2)
    edges = collect(0.0:500.0:3500.0)
    rows = NamedTuple[]
    for (lower, upper) in zip(edges[1:end-1], edges[2:end])
        selected = (surface_height .>= lower) .& (surface_height .< upper)
        !any(selected) && continue
        area = area_km2[selected]
        dense_depth = firn.depth_to_dense_layer[selected]
        finite_dense = isfinite.(dense_depth)
        push!(rows, (
            lower=lower,
            upper=upper,
            refreezing_gt=sum(refreezing[selected] .* area) * 1e-6,
            cold_content_gt=sum(firn.cold_content[selected] .* area) * 1e-6,
            pore_capacity_gt=sum(firn.pore_capacity[selected] .* area) * 1e-6,
            liquid_water_gt=sum(firn.liquid_water[selected] .* area) * 1e-6,
            dense_layer_area_fraction=sum(area[finite_dense]) / sum(area),
            dense_layer_depth_m=any(finite_dense) ? sum(dense_depth[finite_dense] .* area[finite_dense]) / sum(area[finite_dense]) : NaN,
        ))
    end
    return rows
end

function write_firn_state_summary(path, rows)
    open(path, "w") do io
        println(io, "elevation_lower_m,elevation_upper_m,chion_refreezing_gt,cold_content_capacity_gt,pore_capacity_gt,liquid_water_gt,dense_layer_area_fraction,mean_dense_layer_depth_m")
        for row in rows
            println(io, "$(row.lower),$(row.upper),$(row.refreezing_gt),$(row.cold_content_gt),$(row.pore_capacity_gt),$(row.liquid_water_gt),$(row.dense_layer_area_fraction),$(row.dense_layer_depth_m)")
        end
    end
    return nothing
end

function pre_melt_firn_figure(grid, firn, case)
    panels = (
        map_panel(grid, firn.cold_content, "1 May cold-content capacity (mmWE)"; coastline_level=CONFIG.mask_threshold),
        map_panel(grid, firn.pore_capacity, "1 May pore capacity (mmWE)"; coastline_level=CONFIG.mask_threshold),
        map_panel(grid, firn.liquid_water, "1 May liquid water (mmWE)"; coastline_level=CONFIG.mask_threshold),
        map_panel(grid, firn.depth_to_dense_layer, "1 May depth to density ≥810 kg m⁻³ (m)"; coastline_level=CONFIG.mask_threshold),
    )
    return plot(panels...; layout=(2, 2), size=(1200, 1100), title="$(case.title): pre-melt firn state")
end

"""Allocate column-wise, diagnostic-only hydrology pathway accounting."""
function hydrology_pathway_buffers(ncol)
    z() = zeros(Float64, ncol)
    return (
        liquid_input=z(), refreezing=z(), runoff=z(),
        runoff_preexisting_bare_ice=z(), runoff_snow_covered=z(),
        runoff_on_snow_disappearance_day=z(),
        input_with_cold_content=z(), input_without_cold_content=z(),
        input_with_pore_capacity=z(), snow_disappearance_count=z(),
    )
end

function monthly_hydrology_buffers(ncol)
    z() = zeros(Float64, 12, ncol)
    return (liquid_input=z(), refreezing=z(), runoff=z(), cold_overlap=z(), no_cold_overlap=z())
end

function accumulate_monthly_hydrology!(buffers, time, before, melt, rainfall, refreezing, runoff)
    month_index = month(time)
    liquid_input = max.(melt .+ rainfall, 0.0)
    cold_overlap = min.(liquid_input, before.cold_content)
    @views buffers.liquid_input[month_index, :] .+= liquid_input
    @views buffers.refreezing[month_index, :] .+= max.(refreezing, 0.0)
    @views buffers.runoff[month_index, :] .+= max.(runoff, 0.0)
    @views buffers.cold_overlap[month_index, :] .+= cold_overlap
    @views buffers.no_cold_overlap[month_index, :] .+= liquid_input .- cold_overlap
    return nothing
end

function write_monthly_hydrology_elevation_summary(path, hydrology, budgets, surface_height, area_km2)
    edges = collect(0.0:500.0:3500.0)
    open(path, "w") do io
        println(io, "month,elevation_lower_m,elevation_upper_m,mar_melt_gt,chion_melt_gt,mar_refreezing_gt,chion_refreezing_gt,liquid_input_gt,cold_overlap_gt,no_cold_overlap_gt,runoff_gt")
        for month_index in 1:12, (lower, upper) in zip(edges[1:end-1], edges[2:end])
            selected = (surface_height .>= lower) .& (surface_height .< upper)
            !any(selected) && continue
            area = area_km2[selected]
            gt(values) = sum(values[selected] .* area) * 1e-6
            println(io, join((
                month_index, lower, upper,
                gt(@view budgets.truth[:melt][month_index, :]),
                gt(@view budgets.configured[:melt][month_index, :]),
                gt(@view budgets.truth[:refreezing][month_index, :]),
                gt(@view budgets.configured[:refreezing][month_index, :]),
                gt(@view hydrology.liquid_input[month_index, :]),
                gt(@view hydrology.cold_overlap[month_index, :]),
                gt(@view hydrology.no_cold_overlap[month_index, :]),
                gt(@view hydrology.runoff[month_index, :]),
            ), ','))
        end
    end
    return nothing
end

function write_spatial_melt_cold_overlap(path, firn, mar_melt, chion_melt, mar_refreezing, chion_refreezing, surface_height, area_km2)
    edges = collect(0.0:500.0:3500.0)
    open(path, "w") do io
        println(io, "elevation_lower_m,elevation_upper_m,mar_melt_gt,chion_melt_gt,may_cold_content_gt,mar_melt_may_cold_potential_gt,chion_melt_may_cold_potential_gt,mar_melt_cold_correlation,chion_melt_cold_correlation,area_mar_melt_without_chion_fraction")
        for (lower, upper) in zip(edges[1:end-1], edges[2:end])
            selected = (surface_height .>= lower) .& (surface_height .< upper)
            !any(selected) && continue
            area = area_km2[selected]
            mm_to_gt(values) = sum(values .* area) * 1e-6
            cold = firn.cold_content[selected]
            mar = mar_melt[selected]
            chion = chion_melt[selected]
            missed = (mar .> 1.0) .& (chion .<= 1.0)
            println(io, join((
                lower, upper, mm_to_gt(mar), mm_to_gt(chion), mm_to_gt(cold),
                mm_to_gt(min.(mar, cold)), mm_to_gt(min.(chion, cold)),
                length(mar) > 1 ? cor(mar, cold) : NaN,
                length(chion) > 1 ? cor(chion, cold) : NaN,
                sum(area[missed]) / sum(area),
            ), ','))
        end
    end
    return nothing
end

function write_fine_elevation_melt_refreezing(path, mar_melt, chion_melt, mar_refreezing, chion_refreezing, surface_height, area_km2)
    edges = collect(0.0:100.0:3500.0)
    open(path, "w") do io
        println(io, "elevation_lower_m,elevation_upper_m,mar_melt_gt,chion_melt_gt,melt_bias_gt,mar_refreezing_gt,chion_refreezing_gt,refreezing_bias_gt")
        for (lower, upper) in zip(edges[1:end-1], edges[2:end])
            selected = (surface_height .>= lower) .& (surface_height .< upper)
            !any(selected) && continue
            area = area_km2[selected]
            gt(values) = sum(values[selected] .* area) * 1e-6
            mar_melt_gt, chion_melt_gt = gt(mar_melt), gt(chion_melt)
            mar_refreezing_gt, chion_refreezing_gt = gt(mar_refreezing), gt(chion_refreezing)
            println(io, join((lower, upper, mar_melt_gt, chion_melt_gt,
                chion_melt_gt - mar_melt_gt, mar_refreezing_gt, chion_refreezing_gt,
                chion_refreezing_gt - mar_refreezing_gt), ','))
        end
    end
    return nothing
end

function write_energy_component_elevation_summary(path, comparison, energy, surface_height, area_km2)
    edges = collect(0.0:500.0:3500.0)
    open(path, "w") do io
        println(io, "elevation_lower_m,elevation_upper_m,mar_lw_down_wm2,chion_lw_down_wm2,mar_sensible_wm2,chion_sensible_wm2,mar_net_longwave_pj,chion_net_longwave_pj,mar_sensible_pj,chion_sensible_pj")
        month_days = comparison.month_days
        for (lower, upper) in zip(edges[1:end-1], edges[2:end])
            selected = (surface_height .>= lower) .& (surface_height .< upper)
            !any(selected) && continue
            area = area_km2[selected]
            annual_flux(values) = sum(values[:, selected] .* reshape(month_days, :, 1) .* reshape(area, 1, :)) /
                                  (sum(month_days) * sum(area))
            pj(values) = sum(values[:, selected] .* reshape(area, 1, :)) * 1e-9
            println(io, join((
                lower, upper,
                annual_flux(comparison.truth[:q_lw_down]),
                annual_flux(comparison.configured[:q_lw_down]),
                annual_flux(comparison.truth[:q_sh]),
                annual_flux(comparison.configured[:q_sh]),
                pj(energy.truth[:longwave]), pj(energy.configured[:longwave]),
                pj(energy.truth[:sensible_heat]), pj(energy.configured[:sensible_heat]),
            ), ','))
        end
    end
    return nothing
end

const MELT_EVENT_FIELDS = (
    :mar_melt_rate, :chion_melt_rate,
    :mar_shortwave, :chion_shortwave,
    :mar_lw_down, :chion_lw_down,
    :mar_longwave, :chion_longwave,
    :mar_sensible, :chion_sensible,
    :mar_latent, :chion_latent,
    :mar_albedo, :chion_albedo,
    :air_temperature, :relative_humidity, :air_vapor_pressure,
    :wind_speed, :mar_surface_temperature, :chion_surface_temperature,
    :mar_air_surface_temperature_difference,
    :chion_air_surface_temperature_difference,
    :mar_heating_chion_cooling_fraction,
    :pre_event_cold_capacity,
)

function melt_event_energy_buffers()
    edges = collect(1500.0:100.0:3500.0)
    nbins = length(edges) - 1
    allocate() = Dict(field => zeros(Float64, nbins) for field in MELT_EVENT_FIELDS)
    return (
        edges=edges,
        event_weight=zeros(Float64, nbins),
        missed_weight=zeros(Float64, nbins),
        event=allocate(),
        missed=allocate(),
    )
end

function accumulate_melt_event_energy!(
    diagnostic, surface_height, area_km2, dt_days, mar_melt, chion_melt,
    mar_flux, chion_flux, mar_lw_down, chion_lw_down, mar_albedo, chion_albedo,
    air_temperature, chion_surface_temperature, pre_event_cold_capacity,
    relative_humidity, wind_speed,
)
    threshold = 0.01
    mar_longwave_up = max.(mar_lw_down .- mar_flux.longwave, 0.0)
    mar_surface_temperature = (mar_longwave_up ./ (0.98 * 5.670373e-8)) .^ 0.25
    values = (
        mar_melt_rate=mar_melt ./ dt_days,
        chion_melt_rate=chion_melt ./ dt_days,
        mar_shortwave=mar_flux.shortwave,
        chion_shortwave=chion_flux.shortwave,
        mar_lw_down=mar_lw_down,
        chion_lw_down=chion_lw_down,
        mar_longwave=mar_flux.longwave,
        chion_longwave=chion_flux.longwave,
        mar_sensible=mar_flux.sensible_heat,
        chion_sensible=chion_flux.sensible_heat,
        mar_latent=mar_flux.latent_heat,
        chion_latent=chion_flux.latent_heat,
        mar_albedo=mar_albedo,
        chion_albedo=chion_albedo,
        air_temperature=air_temperature,
        relative_humidity=relative_humidity,
        air_vapor_pressure=[
            Chion._bessi_air_vapor_pressure(
                air_temperature[column], relative_humidity[column], 273.15,
            ) for column in eachindex(air_temperature)
        ],
        wind_speed=wind_speed,
        mar_surface_temperature=mar_surface_temperature,
        chion_surface_temperature=chion_surface_temperature,
        mar_air_surface_temperature_difference=air_temperature .- mar_surface_temperature,
        chion_air_surface_temperature_difference=air_temperature .- chion_surface_temperature,
        mar_heating_chion_cooling_fraction=
            Float64.((mar_flux.sensible_heat .> 0.0) .& (chion_flux.sensible_heat .< 0.0)),
        pre_event_cold_capacity=pre_event_cold_capacity,
    )
    for bin in eachindex(diagnostic.event_weight)
        lower, upper = diagnostic.edges[bin], diagnostic.edges[bin + 1]
        in_bin = (surface_height .>= lower) .& (surface_height .< upper)
        event = in_bin .& (mar_melt .> threshold)
        missed = event .& (chion_melt .<= threshold)
        for (selected, weight_sum, sums) in (
            (event, diagnostic.event_weight, diagnostic.event),
            (missed, diagnostic.missed_weight, diagnostic.missed),
        )
            any(selected) || continue
            weights = area_km2[selected] .* dt_days
            weight_sum[bin] += sum(weights)
            for field in MELT_EVENT_FIELDS
                sums[field][bin] += sum(getproperty(values, field)[selected] .* weights)
            end
        end
    end
    return nothing
end

function write_melt_event_energy_summary(path, diagnostic)
    open(path, "w") do io
        columns = String["elevation_lower_m", "elevation_upper_m", "category", "area_days_km2"]
        append!(columns, String.(MELT_EVENT_FIELDS))
        println(io, join(columns, ','))
        for bin in eachindex(diagnostic.event_weight)
            for (category, weight, sums) in (
                ("mar_melt_event", diagnostic.event_weight[bin], diagnostic.event),
                ("missed_by_chion", diagnostic.missed_weight[bin], diagnostic.missed),
            )
                weight <= 0 && continue
                means = (sums[field][bin] / weight for field in MELT_EVENT_FIELDS)
                println(io, join((
                    diagnostic.edges[bin], diagnostic.edges[bin + 1], category, weight, means...,
                ), ','))
            end
        end
    end
    return nothing
end

function melt_cold_overlap_figure(grid, mar_melt, chion_melt, firn, case)
    difference = chion_melt .- mar_melt
    difference_extent = maximum(abs, difference)
    melt_upper = max(maximum(mar_melt), maximum(chion_melt))
    panels = (
        map_panel(grid, mar_melt, "MAR annual melt"; coastline_level=CONFIG.mask_threshold, colorbar_label="mmWE", clims=(0.0, melt_upper)),
        map_panel(grid, chion_melt, "Chion annual melt"; coastline_level=CONFIG.mask_threshold, colorbar_label="mmWE", clims=(0.0, melt_upper)),
        map_panel(grid, difference, "Chion − MAR melt"; coastline_level=CONFIG.mask_threshold, colorbar_label="mmWE", clims=(-difference_extent, difference_extent)),
        map_panel(grid, firn.cold_content, "Chion 1 May cold capacity"; coastline_level=CONFIG.mask_threshold, colorbar_label="mmWE"),
    )
    return plot(panels...; layout=(2, 2), size=(1200, 1100), title="$(case.title): melt–cold-content overlap")
end

function accumulate_hydrology_pathways!(diagnostic, before, n_before, n_after, melt, rainfall, refreezing, runoff)
    liquid_input = max.(melt .+ rainfall, 0.0)
    diagnostic.liquid_input .+= liquid_input
    diagnostic.refreezing .+= max.(refreezing, 0.0)
    diagnostic.runoff .+= max.(runoff, 0.0)

    preexisting_bare = n_before .== 0
    disappeared = (n_before .> 0) .& (n_after .== 0)
    remained_snow = (n_before .> 0) .& .!disappeared
    diagnostic.runoff_preexisting_bare_ice[preexisting_bare] .+= max.(runoff[preexisting_bare], 0.0)
    diagnostic.runoff_on_snow_disappearance_day[disappeared] .+= max.(runoff[disappeared], 0.0)
    diagnostic.runoff_snow_covered[remained_snow] .+= max.(runoff[remained_snow], 0.0)
    diagnostic.snow_disappearance_count[disappeared] .+= 1

    # These are conservative daily overlap tests, not claims about the exact
    # substep pathway: how much liquid input entered a column that had enough
    # cold content / pore volume at the beginning of that forcing interval.
    cold_overlap = min.(liquid_input, before.cold_content)
    diagnostic.input_with_cold_content .+= cold_overlap
    diagnostic.input_without_cold_content .+= liquid_input .- cold_overlap
    diagnostic.input_with_pore_capacity .+= min.(liquid_input, before.pore_capacity)
    return nothing
end

function hydrology_pathway_elevation_summary(diagnostic, surface_height, area_km2)
    edges = collect(0.0:500.0:3500.0)
    rows = NamedTuple[]
    for (lower, upper) in zip(edges[1:end-1], edges[2:end])
        selected = (surface_height .>= lower) .& (surface_height .< upper)
        !any(selected) && continue
        area = area_km2[selected]
        gt(field) = sum(field[selected] .* area) * 1e-6
        push!(rows, (
            lower=lower, upper=upper,
            liquid_input_gt=gt(diagnostic.liquid_input),
            refreezing_gt=gt(diagnostic.refreezing), runoff_gt=gt(diagnostic.runoff),
            bare_ice_runoff_gt=gt(diagnostic.runoff_preexisting_bare_ice),
            snow_covered_runoff_gt=gt(diagnostic.runoff_snow_covered),
            disappearance_day_runoff_gt=gt(diagnostic.runoff_on_snow_disappearance_day),
            cold_overlap_gt=gt(diagnostic.input_with_cold_content),
            no_cold_overlap_gt=gt(diagnostic.input_without_cold_content),
            pore_overlap_gt=gt(diagnostic.input_with_pore_capacity),
            disappearance_area_fraction=sum(area[diagnostic.snow_disappearance_count[selected] .> 0]) / sum(area),
        ))
    end
    return rows
end

function write_hydrology_pathway_summary(path, rows)
    names = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(names, ','))
        for row in rows
            println(io, join((getproperty(row, name) for name in names), ','))
        end
    end
    return nothing
end

"""Run a spinup, then collect one full year of diagnostic comparisons."""
function run_case(grid, base_forcing, reference, case)
    forcing = forcing_for_case(base_forcing, reference, case.prescribed)
    model = model_for_case(grid, case)
    spinup = Simulation(model; forcing, years=CONFIG.spinup_years, backend=CONFIG.backend, write_netcdf=false)
    run!(spinup)

    diagnostic = Simulation(model; forcing, state=spinup.now, years=1, backend=CONFIG.backend, write_netcdf=false)
    integrator = init_integrator(diagnostic; io=devnull)
    runtime = integrator.model_runtime.data
    buffers = comparison_buffers(length(grid.js))
    budgets = budget_buffers(length(grid.js))
    energy = energy_buffers(length(grid.js))
    latent_diagnostics = latent_diagnostic_buffers(length(grid.js))
    initial_budget = annual_budget_start(runtime.state)
    previous_budget = annual_budget_start(runtime.state)
    # The latent flux is evaluated within every diurnal substep.  Retain its
    # accumulated W m⁻² day diagnostic so the comparison products use the
    # same time-integrated flux that drives vapor exchange in the model.
    previous_latent_heat_flux_sum = copy(to_host(runtime.state.latent_heat_flux_sum))
    previous_firn_mass = firn_mass(runtime.state)
    previous_net_longwave_energy = zeros(Float64, length(grid.js))
    previous_absorbed_shortwave_energy = zeros(Float64, length(grid.js))
    previous_sensible_heat_energy = zeros(Float64, length(grid.js))
    previous_rain_heat_energy = zeros(Float64, length(grid.js))
    monthly_delta_firn = zeros(Float64, 12, length(grid.js))
    hydrology_pathways = hydrology_pathway_buffers(length(grid.js))
    monthly_hydrology = monthly_hydrology_buffers(length(grid.js))
    melt_event_energy = melt_event_energy_buffers()
    firn_pre_melt = nothing
    c = model.c
    for time_index in eachindex(forcing.time_values)
        # Diagnose against the state entering this forcing interval. Using the
        # post-step temperature would compare MAR's flux at t with Chion's
        # surface state at t + Δt and artificially degrade agreement.
        state = runtime.state
        n_before = to_host(state.N)
        firn_before = pre_melt_firn_state(state)
        if month(forcing.time_values[time_index]) == 5 && day(forcing.time_values[time_index]) == 1
            firn_pre_melt = pre_melt_firn_state(state)
        end
        air_temperature = to_host(forcing.air_temperature[:, time_index])
        surface_temperature = to_host(state.Tsrf)
        surface_albedo = to_host(state.albedo)
        shortwave_down = to_host(forcing.shortwave_down[:, time_index])
        wind_speed = to_host(forcing.wind_speed[:, time_index])
        relative_humidity = to_host(forcing.relative_humidity[:, time_index])
        has_relative_humidity = to_host(forcing.has_relative_humidity[:, time_index])
        air_pressure = to_host(forcing.air_pressure[:, time_index])
        parameterized = (
            q_sw_net=(one(eltype(surface_albedo)) .- surface_albedo) .* shortwave_down,
            albedo=surface_albedo,
            q_lw_down=c.σ .* c.ϵ_air .* air_temperature .^ 4,
            # Replaced after stepping by the surface-solver-integrated flux.
            q_sh=zeros(eltype(surface_temperature), length(surface_temperature)),
            q_lh=[Chion._resolved_turbulent_latent_heat_flux(
                c,
                surface_temperature[column],
                air_temperature[column],
                false,
                zero(eltype(surface_temperature)),
                has_relative_humidity[column],
                relative_humidity[column],
                air_pressure[column],
                wind_speed[column],
            ) for column in eachindex(surface_temperature)],
        )
        truth = (
            q_sw_net=reference.q_sw_net[:, time_index],
            albedo=reference.albedo[:, time_index],
            q_lw_down=to_host(forcing.q_lw_down[:, time_index]),
            q_sh=to_host(forcing.q_sh[:, time_index]),
            q_lh=to_host(forcing.q_lh[:, time_index]),
        )
        dt_seconds = forcing.dt_days[time_index] * c.seconds_per_day
        mar_energy_flux = (
            shortwave=reference.energy.shortwave[:, time_index],
            longwave=reference.energy.longwave[:, time_index],
            sensible_heat=truth.q_sh,
            latent_heat=truth.q_lh,
            rain_heat=to_host(forcing.rainfall_rate[:, time_index]) .* c.cw .* (air_temperature .- c.T0),
        )
        Chion.step_model!(model, diagnostic.now, integrator.model_runtime, runtime.step_fields, time_index)
        # Unlike the other forcing diagnostics, turbulent latent heat is not
        # a single daily flux under diurnal stepping.  Convert the actual
        # substep-integrated model diagnostic back to its daily-mean flux.
        interval_latent_heat_energy =
            (to_host(runtime.state.latent_heat_flux_sum) .-
             previous_latent_heat_flux_sum) .* c.seconds_per_day
        actual_latent_heat_flux = interval_latent_heat_energy ./ dt_seconds
        actual_surface_temperature = to_host(runtime.state.Tsrf)
        vapor_pressure_deficit = [
            Chion._bessi_air_vapor_pressure(
                air_temperature[column], relative_humidity[column], c.T0,
            ) - Chion._bessi_ice_saturation_vapor_pressure(
                actual_surface_temperature[column], c.T0,
            ) for column in eachindex(actual_surface_temperature)
        ]
        vapor_conductance = [
            Chion._surface_vapor_latent_heat(actual_surface_temperature[column], c) *
            Chion._bessi_vapor_exchange_coefficient(c) /
            max(air_pressure[column], eps(eltype(air_pressure)))
            for column in eachindex(actual_surface_temperature)
        ]
        net_longwave_energy = to_host(runtime.workspace.net_longwave_energy)
        interval_net_longwave_energy = net_longwave_energy .- previous_net_longwave_energy
        copyto!(previous_net_longwave_energy, net_longwave_energy)
        shortwave_energy = to_host(runtime.workspace.absorbed_shortwave_energy)
        interval_shortwave_energy = shortwave_energy .- previous_absorbed_shortwave_energy
        copyto!(previous_absorbed_shortwave_energy, shortwave_energy)
        sensible_heat_energy = to_host(runtime.workspace.sensible_heat_energy)
        interval_sensible_heat_energy = sensible_heat_energy .- previous_sensible_heat_energy
        copyto!(previous_sensible_heat_energy, sensible_heat_energy)
        rain_heat_energy = to_host(runtime.workspace.rain_heat_energy)
        interval_rain_heat_energy = rain_heat_energy .- previous_rain_heat_energy
        copyto!(previous_rain_heat_energy, rain_heat_energy)
        actual_shortwave_flux = interval_shortwave_energy ./ dt_seconds
        actual_sensible_heat_flux = interval_sensible_heat_energy ./ dt_seconds
        actual_rain_heat_flux = interval_rain_heat_energy ./ dt_seconds
        configured = (
            q_sw_net=case.prescribed.q_sw_net ? truth.q_sw_net : actual_shortwave_flux,
            albedo=(CONFIG.force_aging_albedo || !case.prescribed.albedo) ? parameterized.albedo : truth.albedo,
            q_lw_down=case.prescribed.q_lw_down ? truth.q_lw_down : parameterized.q_lw_down,
            q_sh=case.prescribed.q_sh ? truth.q_sh : actual_sensible_heat_flux,
            q_lh=case.prescribed.q_lh ? truth.q_lh : actual_latent_heat_flux,
        )
        accumulate_comparison!(buffers, forcing.time_values[time_index], forcing.dt_days[time_index], truth, configured)
        accumulate_latent_diagnostics!(
            latent_diagnostics,
            forcing.time_values[time_index],
            forcing.dt_days[time_index],
            vapor_pressure_deficit,
            vapor_conductance,
            actual_latent_heat_flux,
            truth.q_lh,
        )
        chion_energy_flux = (
            shortwave=actual_shortwave_flux,
            longwave=interval_net_longwave_energy ./ dt_seconds,
            sensible_heat=actual_sensible_heat_flux,
            latent_heat=actual_latent_heat_flux,
            rain_heat=actual_rain_heat_flux,
        )
        cumulative_budget = NamedTuple(
            name => to_host(getfield(runtime.state, name)) .- getfield(previous_budget, name)
            for name in CUMULATIVE_BUDGET_VARIABLES
        )
        rainfall_mass = to_host(forcing.rainfall_rate[:, time_index]) .* dt_seconds
        mar_melt_interval = reference.budgets.melt[:, time_index] .* forcing.dt_days[time_index]
        accumulate_melt_event_energy!(
            melt_event_energy,
            base_forcing.surface_height[:, 1],
            reference.area_km2,
            forcing.dt_days[time_index],
            mar_melt_interval,
            cumulative_budget.melt,
            mar_energy_flux,
            chion_energy_flux,
            truth.q_lw_down,
            configured.q_lw_down,
            truth.albedo,
            configured.albedo,
            air_temperature,
            actual_surface_temperature,
            firn_before.cold_content,
            relative_humidity,
            to_host(forcing.wind_speed[:, time_index]),
        )
        accumulate_hydrology_pathways!(
            hydrology_pathways, firn_before, n_before, to_host(runtime.state.N),
            cumulative_budget.melt, rainfall_mass,
            cumulative_budget.refreezing, cumulative_budget.runoff,
        )
        accumulate_monthly_hydrology!(
            monthly_hydrology, forcing.time_values[time_index], firn_before,
            cumulative_budget.melt, rainfall_mass,
            cumulative_budget.refreezing, cumulative_budget.runoff,
        )
        model_budget = (
            cumulative_budget...,
            surface_balance=(to_host(forcing.snowfall_rate[:, time_index]) .+
                             to_host(forcing.rainfall_rate[:, time_index])) .*
                            (forcing.dt_days[time_index] * model.c.seconds_per_day) .-
                            cumulative_budget.runoff .+
                            cumulative_budget.vapor_mass,
        )
        mar_budget = NamedTuple(
            name => reference.budgets[name][:, time_index] .* forcing.dt_days[time_index]
            for name in keys(BUDGET_VARIABLES)
        )
        accumulate_budget!(budgets, forcing.time_values[time_index], forcing.dt_days[time_index], mar_budget, model_budget)
        mar_energy = (
            mar_energy_flux...,
            net_surface_energy=mar_energy_flux.shortwave .+
                               mar_energy_flux.longwave .+
                               mar_energy_flux.sensible_heat .+
                               mar_energy_flux.latent_heat .+
                               mar_energy_flux.rain_heat,
            melt_energy=reference.budgets.melt[:, time_index] .* forcing.dt_days[time_index] .* c.Lm,
        )
        chion_energy = (
            chion_energy_flux...,
            net_surface_energy=chion_energy_flux.shortwave .+
                               chion_energy_flux.longwave .+
                               chion_energy_flux.sensible_heat .+
                               chion_energy_flux.latent_heat .+
                               chion_energy_flux.rain_heat,
            melt_energy=cumulative_budget.melt .* c.Lm,
        )
        interval_energy_truth = NamedTuple(name => mar_energy[name] .* dt_seconds for name in keys(ENERGY_VARIABLES))
        interval_energy_chion = NamedTuple(
            name => name == :latent_heat ? interval_latent_heat_energy : chion_energy[name] .* dt_seconds
            for name in keys(ENERGY_VARIABLES)
        )
        interval_energy_truth = (; interval_energy_truth..., melt_energy=mar_energy.melt_energy)
        interval_energy_chion = (; interval_energy_chion..., melt_energy=chion_energy.melt_energy)
        accumulate_energy!(energy, forcing.time_values[time_index], interval_energy_truth, interval_energy_chion)
        firn_mass_change = firn_mass(runtime.state) .- previous_firn_mass
        @views monthly_delta_firn[month(forcing.time_values[time_index]), :] .+= firn_mass_change
        for name in CUMULATIVE_BUDGET_VARIABLES
            copyto!(getfield(previous_budget, name), getfield(runtime.state, name))
        end
        copyto!(previous_firn_mass, firn_mass(runtime.state))
        copyto!(previous_latent_heat_flux_sum, to_host(runtime.state.latent_heat_flux_sum))
    end
    update_diagnostics!(diagnostic.now)
    isnothing(firn_pre_melt) && error("No 1 May pre-melt diagnostic snapshot was found in the forcing calendar.")
    annual = annual_budget_change(runtime.state, initial_budget)
    surface_smb = annual.smb_ice .+ vec(sum(monthly_delta_firn; dims=1))
    return (
        comparison=finalize_comparison(buffers),
        budgets=(; finalize_budget(budgets)..., delta_firn=monthly_delta_firn),
        energy=energy,
        latent_diagnostics=finalize_latent_diagnostics(latent_diagnostics),
        state=diagnostic.now,
        annual=(; annual..., surface_smb),
        firn_pre_melt=firn_pre_melt,
        hydrology_pathways=hydrology_pathways,
        monthly_hydrology=monthly_hydrology,
        melt_event_energy=melt_event_energy,
    )
end

function main()
    mkpath(CONFIG.output_dir)
    requested_case_names = filter(!isempty, split(get(ENV, "CHION_CASES", ""), ','))
    selected_cases = isempty(requested_case_names) ? CASES : filter(CASES) do case
        String(case.name) in requested_case_names
    end
    isempty(selected_cases) && error("`CHION_CASES` did not select any known case.")
    loaded = load_forcing_file(
        CONFIG.forcing_file;
        x_name=CONFIG.x_name,
        y_name=CONFIG.y_name,
        time_name="TIME",
        air_temperature_name="TTZ",
        wind_speed_name=nothing,
        wind_default=CONFIG.wind_speed_m_s,
        mask_name="MSK",
        mask_threshold=CONFIG.mask_threshold,
    )
    grid, base_forcing = loaded.grid, loaded.forcing
    if CONFIG.use_mar_wind_components
        ntime = length(base_forcing.time_values)
        u2z = read_mar_wind_component(CONFIG.forcing_file, "U2Z", grid, ntime)
        v2z = read_mar_wind_component(CONFIG.forcing_file, "V2Z", grid, ntime)
        base_forcing = forcing_with_wind_speed(base_forcing, hypot.(u2z, v2z))
        @info "Using prescribed MAR 2 m wind-speed magnitude" mean_wind_speed_m_s=mean(base_forcing.wind_speed)
    end
    albedo = read_mar_columns(CONFIG.forcing_file, CONFIG.surface_albedo_name, grid, length(base_forcing.time_values))
    albedo .= clamp.(albedo, 0.0, 1.0)
    ntime = length(base_forcing.time_values)
    mar_shortwave_up = read_mar_columns(CONFIG.forcing_file, "SWU", grid, ntime)
    # Antarctic MAR exposes surface temperature rather than LWU. ST2 is in
    # degrees Celsius, so convert it to Kelvin before applying Stefan–Boltzmann
    # to infer LWU for a consistent net-longwave reference.
    mar_longwave_up = has_mar_variable(CONFIG.forcing_file, "LWU") ?
                      read_mar_columns(CONFIG.forcing_file, "LWU", grid, ntime) :
                      0.98 .* 5.670373e-8 .* (read_mar_columns(CONFIG.forcing_file, "ST2", grid, ntime) .+ 273.15) .^ 4
    mar_refreezing = has_mar_variable(CONFIG.forcing_file, "RZ") ?
                     read_mar_columns(CONFIG.forcing_file, "RZ", grid, ntime) :
                     zeros(size(base_forcing.shortwave_down))
    area_km2 = has_mar_variable(CONFIG.forcing_file, "AREA") ?
               read_mar_static_columns(CONFIG.forcing_file, "AREA", grid) :
               regular_grid_area_km2(grid)
    reference = (
        albedo=albedo,
        # Use MAR's diagnosed net shortwave when it is prescribed. AL2 × SWD
        # is retained only as the independent albedo comparison target.
        q_sw_net=base_forcing.shortwave_down .- mar_shortwave_up,
        budgets=(
            melt=read_mar_columns(CONFIG.forcing_file, "ME", grid, ntime),
            runoff=read_mar_columns(CONFIG.forcing_file, "RU", grid, ntime),
            smb_ice=read_mar_columns(CONFIG.forcing_file, "SMB", grid, ntime),
            refreezing=mar_refreezing,
            vapor_mass=-read_mar_columns(CONFIG.forcing_file, "SU", grid, ntime),
            surface_balance=read_mar_columns(CONFIG.forcing_file, "SMB", grid, ntime),
        ),
        energy=(
            shortwave=base_forcing.shortwave_down .- mar_shortwave_up,
            longwave=base_forcing.q_lw_down .- mar_longwave_up,
        ),
        area_km2=area_km2,
    )
    mar_annual = annual_mar_end_state(reference, base_forcing.dt_days)

    for case in selected_cases
        @info "Running MAR--Chion comparison" case=case.name spinup_years=CONFIG.spinup_years antarctica=CONFIG.antarctica_mode
        result = run_case(grid, base_forcing, reference, case)
        annual_integrated_path = joinpath(CONFIG.output_dir, "$(case.name)_annual_integrated_budgets.csv")
        open(annual_integrated_path, "w") do io
            println(io, "variable,MAR_Gt,Chion_Gt,Chion_minus_MAR_Gt")
            for name in (:melt, :runoff, :refreezing)
                mar_total = sum(getfield(mar_annual, name) .* reference.area_km2) * 1e-6
                chion_total = sum(getfield(result.annual, name) .* reference.area_km2) * 1e-6
                println(io, "$(name),$(mar_total),$(chion_total),$(chion_total-mar_total)")
            end
        end
        scatter_path = joinpath(CONFIG.output_dir, "$(case.name)_flux_scatter.pdf")
        budget_path = joinpath(CONFIG.output_dir, "$(case.name)_budget_scatter.pdf")
        integrated_budget_path = joinpath(CONFIG.output_dir, "$(case.name)_monthly_integrated_budgets.pdf")
        energy_path = joinpath(CONFIG.output_dir, "$(case.name)_monthly_integrated_energy.pdf")
        latent_diagnostic_path = joinpath(CONFIG.output_dir, "$(case.name)_latent_flux_diagnostic.pdf")
        state_path = joinpath(CONFIG.output_dir, "$(case.name)_end_state.pdf")
        state_difference_path = joinpath(CONFIG.output_dir, "$(case.name)_end_state_difference.pdf")
        refreezing_path = joinpath(CONFIG.output_dir, "$(case.name)_refreezing_difference.pdf")
        refreezing_summary_path = joinpath(CONFIG.output_dir, "$(case.name)_refreezing_by_elevation.csv")
        firn_state_path = joinpath(CONFIG.output_dir, "$(case.name)_pre_melt_firn_state.pdf")
        firn_summary_path = joinpath(CONFIG.output_dir, "$(case.name)_pre_melt_firn_by_elevation.csv")
        hydrology_pathway_path = joinpath(CONFIG.output_dir, "$(case.name)_hydrology_pathways_by_elevation.csv")
        monthly_hydrology_path = joinpath(CONFIG.output_dir, "$(case.name)_monthly_hydrology_by_elevation.csv")
        spatial_overlap_path = joinpath(CONFIG.output_dir, "$(case.name)_spatial_melt_cold_overlap.csv")
        fine_elevation_path = joinpath(CONFIG.output_dir, "$(case.name)_melt_refreezing_100m_bins.csv")
        melt_cold_map_path = joinpath(CONFIG.output_dir, "$(case.name)_melt_cold_overlap_maps.pdf")
        energy_component_path = joinpath(CONFIG.output_dir, "$(case.name)_energy_components_by_elevation.csv")
        melt_event_energy_path = joinpath(CONFIG.output_dir, "$(case.name)_melt_event_energy_above_1500m.csv")
        yearly_metrics_path = joinpath(CONFIG.output_dir, "$(case.name)_budget_scatter_yearly_metrics.csv")
        latent_summary_path = joinpath(CONFIG.output_dir, "$(case.name)_latent_flux_by_elevation.csv")
        latent_fit_input_path = joinpath(CONFIG.output_dir, "$(case.name)_latent_flux_fit_input.nc")
        mar_refreezing = vec(sum(reference.budgets.refreezing .* reshape(base_forcing.dt_days, 1, :); dims=2))
        mar_melt = vec(sum(reference.budgets.melt .* reshape(base_forcing.dt_days, 1, :); dims=2))
        refreezing_rows = refreezing_elevation_summary(
            mar_refreezing,
            result.annual.refreezing,
            mar_melt,
            result.annual.melt,
            base_forcing.surface_height[:, 1],
            reference.area_km2,
        )
        firn_rows = firn_state_elevation_summary(
            result.firn_pre_melt,
            result.annual.refreezing,
            base_forcing.surface_height[:, 1],
            reference.area_km2,
        )
        hydrology_rows = hydrology_pathway_elevation_summary(
            result.hydrology_pathways,
            base_forcing.surface_height[:, 1],
            reference.area_km2,
        )
        yearly_metrics = yearly_budget_scatter_metrics(result.budgets, reference.area_km2, case)
        latent_rows = latent_flux_elevation_summary(
            result.latent_diagnostics,
            base_forcing.surface_height[:, 1],
            reference.area_km2,
        )
        savefig(comparison_figure(result.comparison, case), scatter_path)
        CairoMakie.save(budget_path, budget_figure(result.budgets, reference.area_km2, case))
        savefig(monthly_integrated_budget_figure(result.budgets, reference.area_km2, case), integrated_budget_path)
        savefig(monthly_integrated_energy_figure(result.energy, reference.area_km2, case), energy_path)
        save_map_figure(latent_diagnostic_path, latent_flux_diagnostic_figure(grid, result.latent_diagnostics, case))
        save_map_figure(state_path, end_state_figure(grid, result.state, result.annual, case))
        save_map_figure(state_difference_path, end_state_difference_figure(grid, mar_annual, result.annual, case))
        savefig(refreezing_difference_figure(grid, mar_refreezing, result.annual.refreezing, case), refreezing_path)
        savefig(pre_melt_firn_figure(grid, result.firn_pre_melt, case), firn_state_path)
        write_refreezing_elevation_summary(refreezing_summary_path, refreezing_rows)
        write_firn_state_summary(firn_summary_path, firn_rows)
        write_hydrology_pathway_summary(hydrology_pathway_path, hydrology_rows)
        write_monthly_hydrology_elevation_summary(
            monthly_hydrology_path, result.monthly_hydrology, result.budgets,
            base_forcing.surface_height[:, 1], reference.area_km2,
        )
        write_spatial_melt_cold_overlap(
            spatial_overlap_path, result.firn_pre_melt,
            mar_melt, result.annual.melt,
            mar_refreezing, result.annual.refreezing,
            base_forcing.surface_height[:, 1], reference.area_km2,
        )
        write_fine_elevation_melt_refreezing(
            fine_elevation_path, mar_melt, result.annual.melt,
            mar_refreezing, result.annual.refreezing,
            base_forcing.surface_height[:, 1], reference.area_km2,
        )
        savefig(melt_cold_overlap_figure(
            grid, mar_melt, result.annual.melt, result.firn_pre_melt, case,
        ), melt_cold_map_path)
        write_energy_component_elevation_summary(
            energy_component_path, result.comparison, result.energy,
            base_forcing.surface_height[:, 1], reference.area_km2,
        )
        write_melt_event_energy_summary(melt_event_energy_path, result.melt_event_energy)
        write_yearly_budget_scatter_metrics(yearly_metrics_path, yearly_metrics)
        write_latent_flux_elevation_summary(latent_summary_path, latent_rows)
        write_latent_flux_fit_input(
            latent_fit_input_path,
            result.latent_diagnostics,
            grid,
            base_forcing.surface_height[:, 1],
            reference.area_km2,
        )
        @info "Refreezing by elevation" rows=refreezing_rows
        @info "Pre-melt firn state by elevation" rows=firn_rows
        @info "Hydrology pathways by elevation" rows=hydrology_rows
        @info "Latent flux by elevation" rows=latent_rows
        @info "Wrote comparison figures" scatter_path budget_path integrated_budget_path energy_path latent_diagnostic_path latent_fit_input_path state_path state_difference_path refreezing_path refreezing_summary_path firn_state_path firn_summary_path hydrology_pathway_path yearly_metrics_path latent_summary_path
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
