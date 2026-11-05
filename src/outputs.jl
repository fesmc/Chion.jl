"""
Output specifications for `Experiment` runs: which variables are written, at which frequency,
how they are aggregated in space, and where.

    Output("smb_monthly.nc"; variables=[:surface_smb, :melt], frequency=:monthly)
    Output("budget.csv"; variables=:budget, frequency=:yearly, aggregate=:area_sum)
    Output("firn.nc"; variables=[:density_profile], frequency=Month(6), depths=0:0.5:20)

Fluxes (melt, runoff, …) are period totals, rates (latent heat flux) are period means, and state
variables are period means (`kind=:mean`, default) or values at the end of the period
(`kind=:snapshot`). Layered and depth-profile variables are always snapshots.
"""

# ---------------------------------------------------------------------------
# Variable catalogue
# ---------------------------------------------------------------------------

"""
`class`:    `:flux` (period total of a cumulative counter or of the forcing), `:rate` (period mean
            of an integrated counter), or `:state`.
`quantity`: what an area sum means: `:mass` (kg m⁻² → Gt), `:fraction` (→ km²), `:length`
            (m → km³), or `:intensive` (area-weighted mean instead of a sum).
`shape`:    `:column`, `:layer` (snow layers), `:substrate` (ice-substrate layers), or `:depth`
            (profile on the output depths).
"""
struct OutputVariable
    name::Symbol
    class::Symbol
    quantity::Symbol
    shape::Symbol
    units::String
    long_name::String
end

const _BESSI_OUTPUT_VARIABLES = Dict{Symbol, OutputVariable}(v.name => v for v in (
    OutputVariable(:surface_smb, :flux, :mass, :column, "kg m-2", "Surface mass balance (precipitation − runoff + vapour exchange)"),
    OutputVariable(:smb_ice, :flux, :mass, :column, "kg m-2", "Net mass forcing to the ice sheet"),
    OutputVariable(:melt, :flux, :mass, :column, "kg m-2", "Melt"),
    OutputVariable(:runoff, :flux, :mass, :column, "kg m-2", "Runoff"),
    OutputVariable(:refreezing, :flux, :mass, :column, "kg m-2", "Refreezing"),
    OutputVariable(:sublimation, :flux, :mass, :column, "kg m-2", "Sublimation"),
    OutputVariable(:vapor_mass, :flux, :mass, :column, "kg m-2", "Net vapour exchange (deposition − sublimation − evaporation)"),
    OutputVariable(:precipitation, :flux, :mass, :column, "kg m-2", "Precipitation"),
    OutputVariable(:snowfall, :flux, :mass, :column, "kg m-2", "Snowfall"),
    OutputVariable(:rainfall, :flux, :mass, :column, "kg m-2", "Rainfall"),
    OutputVariable(:latent_heat_flux, :rate, :intensive, :column, "W m-2", "Turbulent latent heat flux"),
    OutputVariable(:albedo, :state, :intensive, :column, "1", "Surface albedo"),
    OutputVariable(:Tsrf, :state, :intensive, :column, "K", "Surface temperature"),
    OutputVariable(:snow_age_days, :state, :intensive, :column, "day", "Time since the latest snowfall"),
    OutputVariable(:swe, :state, :mass, :column, "kg m-2", "Snow water equivalent (solid and liquid)"),
    OutputVariable(:snow_mass, :state, :mass, :column, "kg m-2", "Solid snow and firn mass"),
    OutputVariable(:liquid_water, :state, :mass, :column, "kg m-2", "Liquid water in the snowpack"),
    OutputVariable(:thickness, :state, :intensive, :column, "m", "Snow and firn depth"),
    OutputVariable(:bulk_density, :state, :intensive, :column, "kg m-3", "Bulk density of snow and firn"),
    OutputVariable(:firn_air_content, :state, :length, :column, "m", "Firn air content"),
    OutputVariable(:snow_cover, :state, :fraction, :column, "1", "Snow cover (SWE at or above the threshold)"),
    OutputVariable(:mass_base, :state, :mass, :column, "kg m-2", "Firn mass exported to the ice below"),
    OutputVariable(:N, :state, :intensive, :column, "1", "Number of snow layers"),
    OutputVariable(:mass, :state, :mass, :layer, "kg m-2", "Layer solid mass"),
    OutputVariable(:mass_w, :state, :mass, :layer, "kg m-2", "Layer liquid water"),
    OutputVariable(:density, :state, :intensive, :layer, "kg m-3", "Layer density"),
    OutputVariable(:temperature, :state, :intensive, :layer, "K", "Layer temperature"),
    OutputVariable(:ice_temperature, :state, :intensive, :substrate, "K", "Ice-substrate layer temperature"),
    OutputVariable(:density_profile, :state, :intensive, :depth, "kg m-3", "Density on depth bins (ice below the firn)"),
    OutputVariable(:temperature_profile, :state, :intensive, :depth, "K", "Temperature at the output depths"),
))

"""Variable selections accepted in place of a list."""
const _OUTPUT_VARIABLE_GROUPS = Dict{Symbol, Vector{Symbol}}(
    :budget => [:surface_smb, :precipitation, :melt, :runoff, :refreezing, :vapor_mass],
    :monthly => [:surface_smb, :smb_ice, :runoff, :melt, :refreezing, :sublimation, :latent_heat_flux, :albedo],
    :default => [:surface_smb, :melt, :runoff, :refreezing, :swe, :albedo],
    :snow => [:swe, :thickness, :snow_cover],
    :firn => [:density_profile, :temperature_profile],
)

# PDD and ITM: their state vectors, with the cumulative ones as fluxes.
const _CUMULATIVE_VECTOR_FIELDS = (:smb_ice, :runoff, :pdd_sum, :smb_cum, :melt_cum, :runoff_cum, :refreezing_cum)
function _vector_model_variable(model, name::Symbol)
    state_fields = output_variables(model)
    name in state_fields || error("`$name` is not an output variable of $(nameof(typeof(model))); " *
        "available: $(join(string.(state_fields), ", ")).")
    meta = get(NETCDF_METADATA, name, OutputMetadata(String(name), String(name), ""))
    class = name in _CUMULATIVE_VECTOR_FIELDS ? :flux : :state
    quantity = meta.units in ("mmWE", "kg m-2") ? :mass : :intensive
    return OutputVariable(name, class, quantity, :column, meta.units, meta.long_name)
end

output_variable(::BESSIModel, name::Symbol) = get(_BESSI_OUTPUT_VARIABLES, name) do
    error("Unknown output variable `$name`. Available: $(join(sort(string.(collect(keys(_BESSI_OUTPUT_VARIABLES)))), ", ")) " *
          "or a group ($(join(string.(collect(keys(_OUTPUT_VARIABLE_GROUPS))), ", "))).")
end
# Forcing totals and the surface mass balance (precipitation − runoff) for every model; ITM's
# instantaneous runoff, melt and refreezing are written as period totals of its cumulative fields.
const _ANY_MODEL_FLUXES = (:precipitation, :snowfall, :rainfall, :surface_smb)
function output_variable(model::AbstractSnowModel, name::Symbol)
    name in _ANY_MODEL_FLUXES && return _BESSI_OUTPUT_VARIABLES[name]
    model isa ITMModel && name in (:runoff, :melt, :refreezing) && return _BESSI_OUTPUT_VARIABLES[name]
    model isa PDDModel && name == :runoff && return _BESSI_OUTPUT_VARIABLES[name]
    return _vector_model_variable(model, name)
end

"""All output variables of `model`, for documentation and error messages."""
available_output_variables(::BESSIModel) = sort(collect(keys(_BESSI_OUTPUT_VARIABLES)))
available_output_variables(model::AbstractSnowModel) = unique([_ANY_MODEL_FLUXES..., output_variables(model)...])

function _expand_variables(spec)
    names = spec isa Symbol || spec isa AbstractString ? [Symbol(spec)] : Symbol.(collect(spec))
    out = Symbol[]
    for n in names
        append!(out, get(_OUTPUT_VARIABLE_GROUPS, n, [n]))
    end
    return unique(out)
end

# ---------------------------------------------------------------------------
# Output specification
# ---------------------------------------------------------------------------

struct Output
    path::String
    variables::Vector{Symbol}
    frequency::Any            # :step, :daily, :monthly, :yearly, or a Dates.Period
    aggregate::Symbol         # :columns, :area_sum, :area_mean
    kind::Symbol              # :mean or :snapshot (state variables)
    period::Union{Nothing, Tuple{DateTime, DateTime}}
    regions::Any              # nothing or NamedTuple of masks
    area::Union{Nothing, Vector{Float64}}
    depths::Vector{Float64}
    snow_cover_threshold::Float64
    columns::Union{Nothing, UnitRange{Int}}   # nothing (all) or the recorded columns
    grid::Union{Nothing, SnowpackGrid}        # grid of `columns`
end

"""
    Output(path; variables=:default, frequency=:monthly, aggregate=:columns, kind=:mean,
           period=nothing, regions=nothing, area=nothing, depths=0:0.5:20, snow_cover_threshold=5.0,
           columns=nothing, grid=nothing)

- `path`: `.nc` (NetCDF) or `.csv` (only for area-aggregated output or a single column).
- `variables`: a list of names or a group (`:budget`, `:monthly`, `:default`, `:snow`, `:firn`);
  see `available_output_variables(model)`.
- `frequency`: `:step`, `:daily`, `:monthly`, `:yearly`, or a period such as `Day(5)`, `Month(6)`.
  Periods are aligned to the calendar (months from January, days from 1970-01-01).
- `aggregate`: `:columns` (one value per column, on the x–y grid when the grid has coordinates),
  `:area_sum` (mass fluxes and stores in Gt, snow cover in km², firn air content in km³,
  intensive variables as area-weighted means) or `:area_mean` (all area-weighted means).
- `kind`: `:mean` (period mean) or `:snapshot` (end of period) for state variables.
- `period`: `(first, last)` years or dates; records outside are not written.
- `regions`: for area aggregation, a NamedTuple of masks (column vectors or y–x matrices),
  e.g. `(eurasia=mask1, north_america=mask2)`; default all columns.
- `area`: column areas in km² (default: from the forcing, e.g. MAR `AREA`).
- `depths`: bin edges (m) of `:density_profile` and depths of `:temperature_profile`.
- `columns`: record only these columns of the domain (a unit range), e.g. one member of an
  ensemble run side by side; `grid` is their `SnowpackGrid` (default: no spatial coordinates),
  and `regions` and `area` refer to them.
"""
function Output(
    path::AbstractString;
    variables=:default,
    frequency=:monthly,
    aggregate::Symbol=:columns,
    kind::Symbol=:mean,
    period=nothing,
    regions=nothing,
    area=nothing,
    depths=0.0:0.5:20.0,
    snow_cover_threshold::Real=5.0,
    columns=nothing,
    grid=nothing,
)
    ext = lowercase(splitext(path)[2])
    ext in (".nc", ".csv") || error("Output path must end in .nc or .csv: $path")
    aggregate in (:columns, :area_sum, :area_mean) || error("`aggregate` must be :columns, :area_sum or :area_mean.")
    kind in (:mean, :snapshot) || error("`kind` must be :mean or :snapshot.")
    _validate_frequency(frequency)
    !isnothing(regions) && aggregate == :columns && error("`regions` requires area aggregation.")
    cols = isnothing(columns) ? nothing : UnitRange{Int}(columns)
    isnothing(grid) || isnothing(cols) || ncols(grid) == length(cols) ||
        error("`grid` has $(ncols(grid)) columns for $(length(cols)) `columns`.")
    isnothing(grid) || !isnothing(cols) || error("`grid` requires `columns`.")
    return Output(String(path), _expand_variables(variables), frequency, aggregate, kind,
        _normalize_period(period), regions, isnothing(area) ? nothing : Float64.(collect(area)),
        Float64.(collect(depths)), Float64(snow_cover_threshold), cols, grid)
end

_validate_frequency(f::Dates.Period) = (Dates.value(f) > 0 || error("Output frequency must be positive."); nothing)
_validate_frequency(f::Symbol) = (f in (:step, :daily, :monthly, :yearly) ||
    error("Output frequency must be :step, :daily, :monthly, :yearly or a Dates.Period, got :$f."); nothing)
_validate_frequency(f::AbstractString) = _validate_frequency(Symbol(f))
_validate_frequency(f) = error("Unsupported output frequency $f.")

_normalize_period(::Nothing) = nothing
_normalize_period(p::AbstractRange{<:Integer}) = (DateTime(first(p)), DateTime(last(p) + 1) - Millisecond(1))
_normalize_period(p::Tuple{<:Integer, <:Integer}) = (DateTime(p[1]), DateTime(p[2] + 1) - Millisecond(1))
_normalize_period(p::Tuple) = (DateTime(p[1]), DateTime(p[2]))

"""Calendar bucket of time `t` for `frequency`; a period ends where the bucket changes."""
_bucket(f::Symbol, t::DateTime, k::Int) =
    f === :step ? k :
    f === :daily ? Dates.value(Date(t)) :
    f === :monthly ? 12 * year(t) + month(t) :
    year(t)
_bucket(f::AbstractString, t, k) = _bucket(Symbol(f), t, k)
_bucket(f::Dates.Month, t::DateTime, k::Int) = (12 * year(t) + month(t) - 1) ÷ Dates.value(f)
_bucket(f::Dates.Year, t::DateTime, k::Int) = fld(year(t), Dates.value(f))
_bucket(f::Dates.Day, t::DateTime, k::Int) = fld(Dates.value(Date(t) - Date(1970)), Dates.value(f))
_bucket(f::Dates.Week, t::DateTime, k::Int) = fld(Dates.value(Date(t) - Date(1970)), 7 * Dates.value(f))
_bucket(f::Dates.Hour, t::DateTime, k::Int) = fld(Dates.value(t - DateTime(1970)), 3_600_000 * Dates.value(f))

# ---------------------------------------------------------------------------
# Accumulation on the backend
# ---------------------------------------------------------------------------

const _COUNTER_OF = Dict(:smb_ice => :smb_ice, :melt => :melt, :runoff => :runoff, :refreezing => :refreezing,
    :sublimation => :sublimation, :vapor_mass => :vapor_mass, :latent_heat_flux => :latent_heat_flux_sum)
const _FORCING_FLUXES = (:surface_smb, :precipitation, :snowfall, :rainfall)

"""The cumulative state field behind the flux `name` of `model`."""
_counter(::BESSIModel, name::Symbol) = get(_COUNTER_OF, name, name)
_counter(::ITMModel, name::Symbol) = name in (:runoff, :melt, :refreezing) ? Symbol(name, :_cum) : name
_counter(::AbstractSnowModel, name::Symbol) = name

"""Cumulative vapour exchange, or zeros for models without it."""
_vapor_counter(state) = hasproperty(state, :vapor_mass) ? state.vapor_mass : zero(_runoff_counter_value(state))
_runoff_counter_value(state) = hasproperty(state, :runoff_cum) ? state.runoff_cum : state.runoff
const _DERIVED_COLUMN = (:swe, :snow_mass, :liquid_water, :thickness, :bulk_density, :firn_air_content, :snow_cover)

"""Per-column derived diagnostics of a BESSI state, written into `out` (a NamedTuple of vectors)."""
@kernel function _derived_column_kernel!(out, fields, rho_ice, threshold)
    idx = @index(Global)
    @inbounds begin
        n = fields.N[idx]
        sample = out.swe[idx]
        solid = zero(sample); liquid = zero(sample); wet = zero(sample); thick = zero(sample); fac = zero(sample)
        for layer_index in 1:n
            m = fields.mass[layer_index, idx]
            w = fields.mass_w[layer_index, idx]
            rho = fields.density[layer_index, idx]
            solid += m
            liquid += w
            wet += m + w                      # summed per layer, as in the state diagnostics
            valid = (m > zero(m)) & (rho > EPS_TINY)
            safe = ifelse(valid, rho, one(rho))
            thick += ifelse(valid, m / safe, zero(m))
            fac += ifelse(valid, m / safe - m / rho_ice, zero(m))
        end
        out.swe[idx] = wet
        out.snow_mass[idx] = solid
        out.liquid_water[idx] = liquid
        out.thickness[idx] = thick
        out.bulk_density[idx] = ifelse(thick > EPS_TINY, solid / ifelse(thick > EPS_TINY, thick, one(thick)), zero(thick))
        out.firn_air_content[idx] = max(fac, zero(fac))
        out.snow_cover[idx] = ifelse(wet >= threshold, one(wet), zero(wet))
    end
end

mutable struct OutputRecorder
    spec::Output
    vars::Vector{OutputVariable}
    model::Any
    flux_prev::Dict{Symbol, Any}       # backend vectors: counters at the period start
    sums::Dict{Symbol, Any}            # backend vectors: Σ value·dt for :mean state variables
    rate_prev::Dict{Symbol, Any}
    precip::Dict{Symbol, Vector{Float64}}   # host: period totals of forcing fluxes (kg m⁻²)
    derived::Any                       # backend NamedTuple of derived-column buffers, or nothing
    days::Float64
    period_start::Union{Nothing, DateTime}
    last_time::Union{Nothing, DateTime}
    needs_every_step::Bool
    needs_host_state::Bool
    writer::Any
    records::Int
    area::Vector{Float64}
    regions::Vector{Pair{String, Vector{Bool}}}
    columns::UnitRange{Int}                 # recorded columns of the domain
end

_needs_derived(rec::OutputRecorder) = any(v -> v.name in _DERIVED_COLUMN, rec.vars)

function _derived_buffers(state, n)
    NF = eltype(state.runoff)
    return NamedTuple{_DERIVED_COLUMN}(ntuple(_ -> fill!(similar(state.runoff, NF, n), zero(NF)), length(_DERIVED_COLUMN)))
end

function _update_derived!(rec::OutputRecorder, state)
    rec.derived === nothing && return nothing
    # Allocate on the device of `state` (the buffers are first made from the host state).
    typeof(rec.derived.swe) == typeof(similar(state.runoff)) || (rec.derived = _derived_buffers(state, ncols(state)))
    k! = _derived_column_kernel!(_ka_backend(state.mass))
    k!(rec.derived, get_fields(state), state.c.rho_i, convert(eltype(state.runoff), rec.spec.snow_cover_threshold);
        ndrange=ncols(state))
    return nothing
end

_column_source(rec::OutputRecorder, state, name::Symbol) =
    name in _DERIVED_COLUMN ? getfield(rec.derived, name) : getproperty(state, name)

function OutputRecorder(spec::Output, model, state, domain_grid, forcing_area, nc_meta)
    columns = something(spec.columns, 1:ncols(domain_grid))
    (first(columns) >= 1 && last(columns) <= ncols(domain_grid)) ||
        error("Output columns $columns lie outside the $(ncols(domain_grid)) columns of the domain.")
    grid = spec.columns === nothing ? domain_grid : something(spec.grid, SnowpackGrid(length(columns)))
    forcing_area = forcing_area === nothing || spec.columns === nothing ? forcing_area : forcing_area[columns]
    vars = [output_variable(model, v) for v in spec.variables]
    model isa BESSIModel || any(v -> v.shape != :column, vars) && error("Only column variables for $(nameof(typeof(model))).")
    if spec.aggregate != :columns
        bad = filter(v -> v.shape != :column, vars)
        isempty(bad) || error("Layer and profile variables cannot be area-aggregated: $(join(string.(getfield.(bad, :name)), ", ")).")
    end
    ext = lowercase(splitext(spec.path)[2])
    ext == ".csv" && spec.aggregate == :columns && ncols(grid) > 1 &&
        error("CSV output needs area aggregation or a single column; use .nc for $(ncols(grid)) columns: $(spec.path)")
    ncol = ncols(grid)
    area = spec.area !== nothing ? spec.area : forcing_area !== nothing ? forcing_area :
        spec.aggregate == :columns ? ones(ncol) :
        error("Area aggregation needs column areas: pass `area` to the Output or provide them in the Forcing.")
    length(area) == ncol || error("`area` has $(length(area)) values for $ncol columns.")
    regions = _region_masks(spec.regions, grid, ncol)
    rec = OutputRecorder(spec, vars, model, Dict(), Dict(), Dict(), Dict(), nothing, 0.0, nothing, nothing,
        false, false, nothing, 0, area, regions, columns)
    for v in vars
        if v.class == :state && v.shape == :column && spec.kind == :mean
            rec.needs_every_step = true
        elseif v.shape in (:layer, :substrate, :depth)
            rec.needs_host_state = true
        end
    end
    model isa BESSIModel && _needs_derived(rec) && (rec.derived = _derived_buffers(state, ncols(state)))   # whole domain
    rec.writer = _open_writer(spec, vars, grid, nc_meta, model, state, length(regions))
    return rec
end

function _region_masks(regions, grid, ncol)
    regions === nothing && return ["all" => trues(ncol)]
    out = Pair{String, Vector{Bool}}[]
    for (name, mask) in pairs(regions)
        push!(out, String(name) => _active_mask_vector(mask, grid, ncol))
    end
    return out
end

"""Start a new period at time `t` from the backend `state`."""
function _start_period!(rec::OutputRecorder, state, t::DateTime)
    for v in rec.vars
        if v.class == :flux && !(v.name in _FORCING_FLUXES)
            rec.flux_prev[v.name] = copy(getproperty(state, _counter(rec.model, v.name)))
        elseif v.class == :rate
            rec.rate_prev[v.name] = copy(getproperty(state, _COUNTER_OF[v.name]))
        end
        if v.name == :surface_smb
            rec.flux_prev[:runoff_for_smb] = copy(_runoff_counter_value(state))
            rec.flux_prev[:vapor_for_smb] = copy(_vapor_counter(state))
        end
        if v.class == :state && v.shape == :column && rec.spec.kind == :mean
            rec.sums[v.name] = fill!(similar(state.runoff, Float64, ncols(state)), 0.0)
        end
    end
    for v in rec.vars
        v.name in _FORCING_FLUXES && (rec.precip[v.name] = zeros(length(rec.columns)))
    end
    rec.days = 0.0
    rec.period_start = t
    return nothing
end

"""Add the forcing fluxes of the host forcing steps `ks` to the period totals."""
function _accumulate_forcing!(rec::OutputRecorder, forcing, ks, seconds_per_day)
    isempty(rec.precip) && return nothing
    sf, rf = forcing.snowfall_rate, forcing.rainfall_rate
    for k in ks
        dt_seconds = Float64(forcing.dt_days[k]) * seconds_per_day
        for (name, acc) in rec.precip
            @inbounds for (i, c) in enumerate(rec.columns)
                s = sf[c, k]; r = rf[c, k]
                acc[i] += (name === :snowfall ? s : name === :rainfall ? r : s + r) * dt_seconds
            end
        end
    end
    return nothing
end

"""After one step of length `dt_days`: add the state means."""
function _accumulate_step!(rec::OutputRecorder, state, dt_days::Float64)
    rec.days += dt_days
    isempty(rec.sums) && return nothing
    _update_derived!(rec, state)
    for (name, acc) in rec.sums
        src = _column_source(rec, state, name)
        acc .+= src .* dt_days
    end
    return nothing
end

"""The recorded columns of a host column vector or `(level, column)` matrix."""
_recorded(rec::OutputRecorder, x::AbstractVector) = length(rec.columns) == length(x) ? x : x[rec.columns]
_recorded(rec::OutputRecorder, x::AbstractMatrix) = length(rec.columns) == size(x, 2) ? x : x[:, rec.columns]

"""Values of every variable for the period that ends now, as host `Float64` vectors (or matrices)
over the recorded columns."""
function _period_values(rec::OutputRecorder, state, host_state)
    out = Dict{Symbol, Any}()
    rec.needs_every_step || _update_derived!(rec, state)
    for v in rec.vars
        name = v.name
        if name == :surface_smb
            acc = rec.precip[:surface_smb]
            Δrunoff = _host_vector(_runoff_counter_value(state) .- rec.flux_prev[:runoff_for_smb])
            Δvapor = _host_vector(_vapor_counter(state) .- rec.flux_prev[:vapor_for_smb])
            out[name] = acc .- _recorded(rec, Δrunoff) .+ _recorded(rec, Δvapor)
            continue
        elseif name in (:precipitation, :snowfall, :rainfall)
            out[name] = copy(rec.precip[name])
            continue
        elseif v.class == :flux
            counter = _counter(rec.model, name)
            out[name] = _host_vector(getproperty(state, counter) .- rec.flux_prev[name])
        elseif v.class == :rate
            Δ = _host_vector(getproperty(state, _COUNTER_OF[name]) .- rec.rate_prev[name])
            out[name] = rec.days > 0 ? Δ ./ rec.days : zero(Δ)
        elseif v.shape == :column && rec.spec.kind == :mean
            out[name] = rec.days > 0 ? _host_vector(rec.sums[name]) ./ rec.days : _host_vector(_column_source(rec, state, name))
        elseif v.shape == :column
            out[name] = _host_vector(_column_source(rec, state, name); copy_array=true)
        elseif v.shape in (:layer, :substrate)
            out[name] = Float64.(Array(getproperty(host_state, name)))           # (layer, column)
        elseif v.shape == :depth
            out[name] = _depth_profile(host_state, name, rec.spec.depths)          # (depth, column)
        end
        out[name] = _recorded(rec, out[name])
    end
    return out
end

"""Density in the depth bins (edges `depths`) or temperature at the depths `depths`, per column."""
function _depth_profile(state::BESSIState, name::Symbol, depths::Vector{Float64})
    N, m, ρ, T, Ti = Array(state.N), Array(state.mass), Array(state.density), Array(state.temperature), Array(state.ice_temperature)
    p = state.parameters
    rho_ice = Float64(p.c.rho_i)
    nsub = size(Ti, 1)
    sub_dz = [p.ice_substrate_top_thickness_m * 2.0^(k - 1) for k in 1:nsub]
    ncol = length(N)
    nout = name == :density_profile ? length(depths) - 1 : length(depths)
    out = fill(NaN, nout, ncol)
    Threads.@threads for c in 1:ncol
        tops = Float64[]; dz = Float64[]; dens = Float64[]; tmp = Float64[]
        z = 0.0
        for l in 1:N[c]
            (m[l, c] > 0 && ρ[l, c] > 0) || continue
            h = m[l, c] / ρ[l, c]
            push!(tops, z); push!(dz, h); push!(dens, ρ[l, c]); push!(tmp, T[l, c]); z += h
        end
        for k in 1:nsub
            push!(tops, z); push!(dz, sub_dz[k]); push!(dens, rho_ice); push!(tmp, Ti[k, c]); z += sub_dz[k]
        end
        if name == :density_profile
            for b in 1:nout
                lo, hi = depths[b], depths[b + 1]
                num = 0.0; den = 0.0
                for i in eachindex(tops)
                    ov = max(0.0, min(hi, tops[i] + dz[i]) - max(lo, tops[i]))
                    num += ov * dens[i]; den += ov
                end
                out[b, c] = den >= hi - lo - 1e-9 ? num / den : den > 0 ? (num + (hi - lo - den) * rho_ice) / (hi - lo) : rho_ice
            end
        else
            isempty(tops) && continue
            mids = tops .+ dz ./ 2
            for (k, d) in enumerate(depths)
                d > last(tops) + last(dz) && continue
                i = searchsortedlast(mids, d)
                out[k, c] = i == 0 ? tmp[1] : i == length(mids) ? tmp[end] :
                    tmp[i] + (tmp[i + 1] - tmp[i]) * (d - mids[i]) / (mids[i + 1] - mids[i])
            end
        end
    end
    return out
end

"""Area aggregation of one column vector over every region."""
function _aggregate(rec::OutputRecorder, v::OutputVariable, values::Vector{Float64})
    out = Float64[]
    for (_, mask) in rec.regions
        if rec.spec.aggregate == :area_sum && v.quantity != :intensive
            scale = v.quantity == :mass ? 1e-6 : v.quantity == :length ? 1e-3 : 1.0   # Gt, km³, km²
            push!(out, sum(values[mask] .* rec.area[mask]) * scale)
        else
            a = rec.area[mask]
            push!(out, sum(values[mask] .* a) / sum(a))
        end
    end
    return out
end

function _aggregate_units(rec::OutputRecorder, v::OutputVariable)
    rec.spec.aggregate == :area_sum || return v.units
    v.quantity == :mass && return "Gt"
    v.quantity == :fraction && return "km2"
    v.quantity == :length && return "km3"
    return v.units
end

"""Close the period that ends at `t_end` (time of the last step) and write one record."""
function _close_period!(rec::OutputRecorder, state, host_state, t_end::DateTime)
    t0 = rec.period_start
    inside = rec.spec.period === nothing || (rec.spec.period[1] <= t0 <= rec.spec.period[2])
    if inside
        values = _period_values(rec, state, host_state)
        rec.records += 1
        _write_record!(rec.writer, rec, values, t0, t_end, rec.records)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

mutable struct CSVWriter
    io::IO
    columns::Vector{Tuple{Symbol, Int}}    # (variable, region index)
end

mutable struct NetCDFWriter
    ds::NCDataset
    handles::Dict{Symbol, Any}
    spatial::Bool
    is::Vector{Int}
    js::Vector{Int}
    nx::Int
    ny::Int
    aggregated::Bool
end

_column_label(rec_regions, v, r, units) =
    (length(rec_regions) == 1 ? String(v) : "$(v)_$(first(rec_regions[r]))") * " ($units)"

function _open_writer(spec::Output, vars, grid, nc_meta, model, state, nregions)
    mkpath(dirname(abspath(spec.path)))
    ext = lowercase(splitext(spec.path)[2])
    if ext == ".csv"
        io = open(spec.path, "w")
        return CSVWriter(io, Tuple{Symbol, Int}[])     # header on the first record (needs region names)
    end
    ds = NCDataset(spec.path, "c")
    defDim(ds, "time", Inf)
    tv = defVar(ds, "time", Float64, ("time",); attrib=Dict("standard_name" => "time", "units" => NETCDF_TIME_UNITS,
        "calendar" => "proleptic_gregorian", "axis" => "T", "bounds" => "time_bnds", "long_name" => "start of the period"))
    defDim(ds, "nv", 2)
    defVar(ds, "time_bnds", Float64, ("nv", "time"); attrib=Dict("units" => NETCDF_TIME_UNITS))
    for (k, v) in nc_meta
        ds.attrib[k] = v
    end
    ds.attrib["output_frequency"] = string(spec.frequency)
    ds.attrib["output_aggregate"] = string(spec.aggregate)
    ds.attrib["state_kind"] = string(spec.kind)
    spatial = spec.aggregate == :columns && has_spatial_coords(grid)
    nx = spatial ? length(grid.x) : 0
    ny = spatial ? length(grid.y) : 0
    if spec.aggregate == :columns
        if spatial
            defDim(ds, "x", nx); defDim(ds, "y", ny)
            defVar(ds, "x", Float64.(grid.x), ("x",); attrib=Dict("long_name" => "x coordinate"))
            defVar(ds, "y", Float64.(grid.y), ("y",); attrib=Dict("long_name" => "y coordinate"))
            cdims = ("x", "y")
        else
            defDim(ds, "column", ncols(grid))
            cdims = ("column",)
        end
    else
        defDim(ds, "region", nregions)
        cdims = ("region",)
    end
    handles = Dict{Symbol, Any}()
    for v in vars
        dims = v.shape == :column ? (cdims..., "time") :
            v.shape == :layer ? (_ensure_dim!(ds, "layer", state.Ntot); (cdims..., "layer", "time")) :
            v.shape == :substrate ? (_ensure_dim!(ds, "substrate_layer", size(state.ice_temperature, 1)); (cdims..., "substrate_layer", "time")) :
            _depth_dims!(ds, v, spec.depths, cdims)
        units = spec.aggregate == :area_sum && v.quantity != :intensive ?
            (v.quantity == :mass ? "Gt" : v.quantity == :fraction ? "km2" : "km3") : v.units
        handles[v.name] = defVar(ds, String(v.name), Float32, dims; fillvalue=NaN32,
            attrib=Dict("units" => units, "long_name" => v.long_name,
                "cell_methods" => "time: " * (v.class == :flux ? "sum" : v.class == :rate ? "mean" :
                    (v.shape == :column && spec.kind == :mean) ? "mean" : "point")))
    end
    return NetCDFWriter(ds, handles, spatial, spatial ? grid.is : Int[], spatial ? grid.js : Int[], nx, ny,
        spec.aggregate != :columns)
end

_ensure_dim!(ds, name, n) = haskey(ds.dim, name) || defDim(ds, name, n)

function _depth_dims!(ds, v, depths, cdims)
    if v.name == :density_profile
        if !haskey(ds.dim, "depth_bin")
            defDim(ds, "depth_bin", length(depths) - 1)
            defVar(ds, "depth_bin", (depths[1:end-1] .+ depths[2:end]) ./ 2, ("depth_bin",);
                attrib=Dict("units" => "m", "long_name" => "centre of the depth bin"))
        end
        return (cdims..., "depth_bin", "time")
    end
    if !haskey(ds.dim, "depth")
        defDim(ds, "depth", length(depths))
        defVar(ds, "depth", depths, ("depth",); attrib=Dict("units" => "m"))
    end
    return (cdims..., "depth", "time")
end

function _write_record!(w::CSVWriter, rec::OutputRecorder, values, t0, t1, r)
    if r == 1
        for v in rec.vars, k in eachindex(rec.regions)
            push!(w.columns, (v.name, k))
        end
        labels = [_column_label(rec.regions, v.name, k, _aggregate_units(rec, v)) for v in rec.vars for k in eachindex(rec.regions)]
        rec.spec.aggregate == :columns && (labels = [String(v.name) * " ($(v.units))" for v in rec.vars])
        println(w.io, join(vcat(["period_start", "period_end"], labels), ','))
    end
    fields = String[string(t0), string(t1)]
    for v in rec.vars
        x = values[v.name]
        vals = rec.spec.aggregate == :columns ? [x[1]] : _aggregate(rec, v, x)
        append!(fields, repr.(vals))
    end
    println(w.io, join(fields, ','))
    flush(w.io)
    return nothing
end

function _write_record!(w::NetCDFWriter, rec::OutputRecorder, values, t0, t1, r)
    w.ds["time"][r] = _netcdf_time_days(t0)
    w.ds["time_bnds"][:, r] = [_netcdf_time_days(t0), _netcdf_time_days(t1)]
    for v in rec.vars
        x = values[v.name]
        h = w.handles[v.name]
        if w.aggregated
            h[:, r] = Float32.(_aggregate(rec, v, x))
        elseif x isa AbstractVector
            if w.spatial
                grid = fill(NaN32, w.nx, w.ny)
                @inbounds for c in eachindex(w.is)
                    grid[w.is[c], w.js[c]] = Float32(x[c])
                end
                h[:, :, r] = grid
            else
                h[:, r] = Float32.(x)
            end
        else                                    # (level, column)
            nl = size(x, 1)
            if w.spatial
                grid = fill(NaN32, w.nx, w.ny, nl)
                @inbounds for c in eachindex(w.is), l in 1:nl
                    grid[w.is[c], w.js[c], l] = Float32(x[l, c])
                end
                h[:, :, :, r] = grid
            else
                h[:, :, r] = permutedims(Float32.(x), (2, 1))
            end
        end
    end
    return nothing
end

function _close_writer!(w::CSVWriter, status)
    close(w.io)
    return nothing
end

function _close_writer!(w::NetCDFWriter, status)
    w.ds.attrib["status"] = string(status)
    close(w.ds)
    return nothing
end
