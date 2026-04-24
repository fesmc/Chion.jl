const OUTPUT_GROUPS = (
    final=(
        :final_thickness,
        :final_wet_mass,
        :final_bulk_density,
        :final_base_mass,
        :final_ice_sheet_smb,
        :final_runoff,
        :last_cycle_delta_thickness,
        :last_cycle_delta_wet_mass,
        :last_cycle_delta_base_mass,
        :last_cycle_delta_ice_sheet_smb,
    ),
    layers=(
        :n_active,
        :layer_density,
        :layer_thickness,
        :layer_snow_mass,
        :layer_liquid_mass,
        :layer_temperature_c,
    ),
    history=(
        :history_mean_thickness,
        :history_mean_wet_mass,
        :history_mean_bulk_density,
        :history_mean_base_mass,
        :history_mean_abs_delta_thickness,
        :history_mean_abs_delta_wet_mass,
        :history_mean_abs_delta_base_mass,
    ),
    monthly=(
        :monthly_mean_thickness,
        :monthly_mean_wet_mass,
        :monthly_mean_bulk_density,
        :monthly_mean_base_mass,
        :monthly_mean_ice_sheet_smb,
        :monthly_export_to_ice,
        :monthly_net_ice_sheet_forcing,
        :monthly_runoff,
    ),
    step=(:step_export_to_ice, :step_ice_sheet_smb),
)
const NETCDF_VARIABLES = unique(Symbol[var for group in values(OUTPUT_GROUPS) for var in group])
const FINAL_GRID_KEYS = OUTPUT_GROUPS.final
const LAYER_GRID_KEYS = OUTPUT_GROUPS.layers
const MONTHLY_GRID_KEYS = OUTPUT_GROUPS.monthly
@inline _grid_shape(layout) = size(layout.mask)

function scatter_to_grid(values::Vector{Float64}, js::Vector{Int}, is::Vector{Int}, grid_shape::Tuple{Int, Int})
    out = fill(NaN, grid_shape)
    @inbounds for idx in eachindex(values)
        out[js[idx], is[idx]] = values[idx]
    end
    return out
end

function monthly_vectors_to_grids(values::Matrix{Float64}, js::Vector{Int}, is::Vector{Int}, grid_shape::Tuple{Int, Int})
    nmonth, nvalid = size(values)
    ny, nx = grid_shape
    out = fill(NaN, nmonth, ny, nx)
    @inbounds for m in 1:nmonth, idx in 1:nvalid
        out[m, js[idx], is[idx]] = values[m, idx]
    end
    return out
end

empty_final_grids() = NamedTuple{FINAL_GRID_KEYS}(ntuple(_ -> Matrix{Float64}(undef, 0, 0), length(FINAL_GRID_KEYS)))
empty_monthly_grids() = NamedTuple{MONTHLY_GRID_KEYS}(ntuple(_ -> Array{Float64}(undef, 0, 0, 0), length(MONTHLY_GRID_KEYS)))

function collect_final_layer_grids(
    domain::SnowpackDomain,
    js::Vector{Int},
    is::Vector{Int},
    grid_shape::Tuple{Int, Int},
    nlayer::Int,
)
    ny, nx = grid_shape
    ncol = length(js)
    n_active = fill(Int32(0), ny, nx)
    layer_density     = fill(NaN, nlayer, ny, nx)
    layer_thickness   = fill(NaN, nlayer, ny, nx)
    layer_snow_mass   = fill(NaN, nlayer, ny, nx)
    layer_liquid_mass = fill(NaN, nlayer, ny, nx)
    layer_temperature_c = fill(NaN, nlayer, ny, nx)
    c = domain.c
    @inbounds for col in 1:ncol
        j, i = js[col], is[col]
        n_active[j, i] = Int32(domain.N[col])
        for k in 1:nlayer
            rho = domain.density[k, col]
            m   = domain.mass[k, col]
            layer_density[k, j, i]       = rho
            layer_snow_mass[k, j, i]     = m
            layer_liquid_mass[k, j, i]   = domain.mass_w[k, col]
            layer_temperature_c[k, j, i] = domain.temperature[k, col] - c.T0
            layer_thickness[k, j, i]     = rho > 0 ? m / rho : 0.0
        end
    end
    return (
        n_active=n_active,
        layer_density=layer_density,
        layer_thickness=layer_thickness,
        layer_snow_mass=layer_snow_mass,
        layer_liquid_mass=layer_liquid_mass,
        layer_temperature_c=layer_temperature_c,
    )
end

function _empty_layer_grids()
    return (
        n_active=Matrix{Int32}(undef, 0, 0),
        layer_density=Array{Float64}(undef, 0, 0, 0),
        layer_thickness=Array{Float64}(undef, 0, 0, 0),
        layer_snow_mass=Array{Float64}(undef, 0, 0, 0),
        layer_liquid_mass=Array{Float64}(undef, 0, 0, 0),
        layer_temperature_c=Array{Float64}(undef, 0, 0, 0),
    )
end

_allocate_step_vectors(active::Bool, ncol::Int) = NamedTuple{OUTPUT_GROUPS.step}(ntuple(_ -> active ? zeros(Float64, ncol) : Float64[], length(OUTPUT_GROUPS.step)))
_allocate_monthly_sums(active::Bool, nmonth_total::Int, ncol::Int) = NamedTuple{MONTHLY_GRID_KEYS}(ntuple(_ -> active ? zeros(Float64, nmonth_total, ncol) : Matrix{Float64}(undef, 0, 0), length(MONTHLY_GRID_KEYS)))

function _step_output_grids(step_vectors, layout)
    return NamedTuple{OUTPUT_GROUPS.step}(ntuple(i -> scatter_to_grid(getfield(step_vectors, OUTPUT_GROUPS.step[i]), layout.js, layout.is, _grid_shape(layout)), length(OUTPUT_GROUPS.step)))
end

function _reset_step_vectors!(step_vectors)
    for key in OUTPUT_GROUPS.step
        isempty(getfield(step_vectors, key)) || fill!(getfield(step_vectors, key), 0.0)
    end
    return
end

function _accumulate_step_diagnostics!(
    summary,
    previous,
    monthly_sums,
    step_vectors,
    month_idx::Int,
    need_monthly_outputs::Bool,
    need_step_outputs::Bool,
)
    current_base = summary.base_mass
    current_smb = summary.smb_ice
    current_runoff = summary.runoff
    delta_base = current_base .- previous.base_mass
    delta_smb = current_smb .- previous.smb_ice
    delta_runoff = current_runoff .- previous.runoff

    if need_monthly_outputs
        monthly_sums.monthly_mean_thickness[month_idx, :] .+= summary.thickness
        monthly_sums.monthly_mean_wet_mass[month_idx, :] .+= summary.wet_mass
        monthly_sums.monthly_mean_bulk_density[month_idx, :] .+= summary.bulk_density
        monthly_sums.monthly_mean_base_mass[month_idx, :] .+= current_base
        monthly_sums.monthly_mean_ice_sheet_smb[month_idx, :] .+= delta_smb
        monthly_sums.monthly_export_to_ice[month_idx, :] .+= delta_base
        monthly_sums.monthly_net_ice_sheet_forcing[month_idx, :] .+= delta_smb
        monthly_sums.monthly_runoff[month_idx, :] .+= delta_runoff
    end
    if need_step_outputs
        step_vectors.step_export_to_ice .+= delta_base
        step_vectors.step_ice_sheet_smb .+= delta_smb
    end
    previous.base_mass .= current_base
    previous.smb_ice .= current_smb
    previous.runoff .= current_runoff
    return
end

function _update_cycle_smb_delta!(last_delta::Vector{Float64}, previous_cycle_smb_ice::Vector{Float64}, domain)
    current = _host_vector(domain.smb_ice; copy_array=true)
    last_delta .= current .- previous_cycle_smb_ice
    previous_cycle_smb_ice .= current
    return
end

function _scatter_final_grids(final_state, domain, deltas, layout)
    final_smb_ice = _host_vector(domain.smb_ice; copy_array=true)
    final_runoff = _host_vector(domain.runoff; copy_array=true)
    return (
        final_thickness=scatter_to_grid(final_state.thickness, layout.js, layout.is, _grid_shape(layout)),
        final_wet_mass=scatter_to_grid(final_state.wet_mass, layout.js, layout.is, _grid_shape(layout)),
        final_bulk_density=scatter_to_grid(final_state.bulk_density, layout.js, layout.is, _grid_shape(layout)),
        final_base_mass=scatter_to_grid(final_state.base_mass, layout.js, layout.is, _grid_shape(layout)),
        final_ice_sheet_smb=scatter_to_grid(final_smb_ice, layout.js, layout.is, _grid_shape(layout)),
        final_runoff=scatter_to_grid(final_runoff, layout.js, layout.is, _grid_shape(layout)),
        last_cycle_delta_thickness=scatter_to_grid(deltas.thickness, layout.js, layout.is, _grid_shape(layout)),
        last_cycle_delta_wet_mass=scatter_to_grid(deltas.wet_mass, layout.js, layout.is, _grid_shape(layout)),
        last_cycle_delta_base_mass=scatter_to_grid(deltas.base_mass, layout.js, layout.is, _grid_shape(layout)),
        last_cycle_delta_ice_sheet_smb=scatter_to_grid(deltas.ice_sheet_smb, layout.js, layout.is, _grid_shape(layout)),
    )
end

const MONTHLY_MEAN_KEYS = (:monthly_mean_thickness, :monthly_mean_wet_mass, :monthly_mean_bulk_density)

function _finalize_monthly_grids(monthly_sums, monthly_count::Vector{Int32}, layout)
    vectors = NamedTuple{MONTHLY_GRID_KEYS}(ntuple(i -> begin
        key = MONTHLY_GRID_KEYS[i]
        data = copy(getfield(monthly_sums, key))
        if key in MONTHLY_MEAN_KEYS
            @inbounds for m in axes(data, 1)
                data[m, :] ./= max(monthly_count[m], 1)
            end
        end
        data
    end, length(MONTHLY_GRID_KEYS)))
    return NamedTuple{MONTHLY_GRID_KEYS}(ntuple(i -> monthly_vectors_to_grids(getfield(vectors, MONTHLY_GRID_KEYS[i]), layout.js, layout.is, _grid_shape(layout)), length(MONTHLY_GRID_KEYS)))
end
