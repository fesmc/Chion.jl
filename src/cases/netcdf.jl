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

@inline _nc_float_attrib(long_name::AbstractString, units::AbstractString) = Dict("long_name" => String(long_name), "units" => String(units))
@inline _nc_text_attrib(long_name::AbstractString) = Dict("long_name" => String(long_name))
_define_nc_output_variable(ds::NCDataset, name::AbstractString, dims; long_name::AbstractString, units::AbstractString) =
    defVar(ds, String(name), Float32, dims; fillvalue=NaN32, attrib=_nc_float_attrib(long_name, units))
_define_nc_int_variable(ds::NCDataset, name::AbstractString, dims; long_name::AbstractString) =
    defVar(ds, String(name), Int32, dims; attrib=_nc_text_attrib(long_name))
_define_nc_double_variable(ds::NCDataset, name::AbstractString, dims; attrib::Dict{String, String}=Dict{String, String}()) =
    defVar(ds, String(name), Float64, dims; attrib=attrib)

const CASE_NC_SPECS = (
    (key=:final_thickness, dims=("y", "x"), name="final_thickness", long_name="Final snow thickness", units="m", integer=false),
    (key=:final_wet_mass, dims=("y", "x"), name="final_wet_mass", long_name="Final snow wet mass", units="mmWE", integer=false),
    (key=:final_bulk_density, dims=("y", "x"), name="final_bulk_density", long_name="Final bulk snow density", units="kg m-3", integer=false),
    (key=:final_base_mass, dims=("y", "x"), name="final_base_mass", long_name="Cumulative firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:final_ice_sheet_smb, dims=("y", "x"), name="final_ice_sheet_smb", long_name="Cumulative net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:final_runoff, dims=("y", "x"), name="final_runoff", long_name="Final cumulative runoff", units="mmWE", integer=false),
    (key=:last_cycle_delta_thickness, dims=("y", "x"), name="last_cycle_delta_thickness", long_name="Last cycle snow-thickness change", units="m", integer=false),
    (key=:last_cycle_delta_wet_mass, dims=("y", "x"), name="last_cycle_delta_wet_mass", long_name="Last cycle wet-mass change", units="mmWE", integer=false),
    (key=:last_cycle_delta_base_mass, dims=("y", "x"), name="last_cycle_delta_base_mass", long_name="Last cycle firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:last_cycle_delta_ice_sheet_smb, dims=("y", "x"), name="last_cycle_delta_ice_sheet_smb", long_name="Last cycle net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:n_active, dims=("y", "x"), name="n_active", long_name="Number of active Chion layers", units="", integer=true),
    (key=:layer_density, dims=("layer", "y", "x"), name="layer_density", long_name="Final Chion layer density", units="kg m-3", integer=false),
    (key=:layer_thickness, dims=("layer", "y", "x"), name="layer_thickness", long_name="Final Chion layer thickness", units="m", integer=false),
    (key=:layer_snow_mass, dims=("layer", "y", "x"), name="layer_snow_mass", long_name="Final Chion layer snow mass", units="kg m-2", integer=false),
    (key=:layer_liquid_mass, dims=("layer", "y", "x"), name="layer_liquid_mass", long_name="Final Chion layer liquid-water mass", units="kg m-2", integer=false),
    (key=:layer_temperature_c, dims=("layer", "y", "x"), name="layer_temperature_c", long_name="Final Chion layer temperature", units="C", integer=false),
    (key=:history_mean_thickness, dims=("cycle",), name="history_mean_thickness", long_name="Cycle-mean snow thickness", units="m", integer=false),
    (key=:history_mean_wet_mass, dims=("cycle",), name="history_mean_wet_mass", long_name="Cycle-mean snow wet mass", units="mmWE", integer=false),
    (key=:history_mean_bulk_density, dims=("cycle",), name="history_mean_bulk_density", long_name="Cycle-mean bulk snow density", units="kg m-3", integer=false),
    (key=:history_mean_base_mass, dims=("cycle",), name="history_mean_base_mass", long_name="Cycle-mean firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:history_mean_abs_delta_thickness, dims=("cycle",), name="history_mean_abs_delta_thickness", long_name="Cycle mean absolute snow-thickness change", units="m", integer=false),
    (key=:history_mean_abs_delta_wet_mass, dims=("cycle",), name="history_mean_abs_delta_wet_mass", long_name="Cycle mean absolute wet-mass change", units="mmWE", integer=false),
    (key=:history_mean_abs_delta_base_mass, dims=("cycle",), name="history_mean_abs_delta_base_mass", long_name="Cycle mean absolute firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:monthly_mean_thickness, dims=("month", "y", "x"), name="monthly_mean_thickness", long_name="Monthly mean snow thickness", units="m", integer=false),
    (key=:monthly_mean_wet_mass, dims=("month", "y", "x"), name="monthly_mean_wet_mass", long_name="Monthly mean snow wet mass", units="mmWE", integer=false),
    (key=:monthly_mean_bulk_density, dims=("month", "y", "x"), name="monthly_mean_bulk_density", long_name="Monthly mean bulk snow density", units="kg m-3", integer=false),
    (key=:monthly_mean_base_mass, dims=("month", "y", "x"), name="monthly_mean_base_mass", long_name="Monthly mean cumulative firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:monthly_mean_ice_sheet_smb, dims=("month", "y", "x"), name="monthly_mean_ice_sheet_smb", long_name="Monthly net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:monthly_export_to_ice, dims=("month", "y", "x"), name="monthly_export_to_ice", long_name="Monthly firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:monthly_net_ice_sheet_forcing, dims=("month", "y", "x"), name="monthly_net_ice_sheet_forcing", long_name="Monthly net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:monthly_runoff, dims=("month", "y", "x"), name="monthly_runoff", long_name="Monthly runoff production", units="mmWE", integer=false),
    (key=:step_export_to_ice, dims=("step", "y", "x"), name="step_export_to_ice", long_name="Annual firn mass exported to the ice model for each written output interval", units="mmWE", integer=false),
    (key=:step_ice_sheet_smb, dims=("step", "y", "x"), name="step_ice_sheet_smb", long_name="Annual net mass forcing to the ice sheet for each written output interval", units="mmWE", integer=false),
)

function _define_selected_nc_variables!(ds::NCDataset, selected::Set{Symbol})
    vars = Dict{Symbol, Any}()
    for spec in CASE_NC_SPECS
        spec.key in selected || continue
        vars[spec.key] = spec.integer ?
            _define_nc_int_variable(ds, spec.name, spec.dims; long_name=spec.long_name) :
            _define_nc_output_variable(ds, spec.name, spec.dims; long_name=spec.long_name, units=spec.units)
    end
    return vars
end

function init_case_netcdf(
    netcdf_path::AbstractString,
    options::RunConfig,
    time_values::Vector{DateTime},
    nlayer::Int,
    layout::GridLayout,
    initial_thickness::Matrix{Float64},
    month_cycle::Vector{Int32},
    month_of_year::Vector{Int32},
    source_month_code::Vector{Int32},
    annual_output_source_indices::Vector{Int32},
    annual_output_source_codes::Vector{Int32},
)
    mkpath(dirname(netcdf_path))
    isdir(netcdf_path) && error("NetCDF output path '$(abspath(netcdf_path))' is a directory; pass a file path ending in `.nc`.")
    ny, nx = _grid_shape(layout)
    max_steps = options.cycles * length(annual_output_source_indices)
    selected = Set(options.netcdf_variables)
    step_cycle = Int32[cyc for cyc in 1:options.cycles for _ in annual_output_source_indices]
    step_source_index = Int32[idx for _ in 1:options.cycles for idx in annual_output_source_indices]
    step_source_code = Int32[code for _ in 1:options.cycles for code in annual_output_source_codes]

    ds = NCDataset(netcdf_path, "c")
    for (name, len) in (("y", ny), ("x", nx), ("layer", max(nlayer, 1)), ("cycle", options.cycles), ("month", length(month_cycle)), ("point", length(layout.js)), ("step", max_steps))
        defDim(ds, name, len)
    end

    _define_nc_double_variable(ds, "x", ("x",); attrib=Dict("units" => "km", "axis" => "X"))[:] = layout.x
    _define_nc_double_variable(ds, "y", ("y",); attrib=Dict("units" => "km", "axis" => "Y"))[:] = layout.y
    _define_nc_int_variable(ds, "layer", ("layer",); long_name="Chion internal layer index from surface downward")[:] = Int32.(collect(1:max(nlayer, 1)))
    _define_nc_int_variable(ds, "cycle", ("cycle",); long_name="Repeated annual forcing cycle index")[:] = Int32.(collect(1:options.cycles))
    _define_nc_int_variable(ds, "month", ("month",); long_name="Sequential monthly output index")[:] = Int32.(collect(1:length(month_cycle)))
    _define_nc_int_variable(ds, "month_cycle", ("month",); long_name="Forcing cycle associated with monthly output")[:] = month_cycle
    _define_nc_int_variable(ds, "month_of_year", ("month",); long_name="Calendar month of the repeated forcing")[:] = month_of_year
    _define_nc_int_variable(ds, "source_month_code", ("month",); long_name="Source forcing month code YYYYMM")[:] = source_month_code
    _define_nc_int_variable(ds, "point", ("point",); long_name="Compact valid cell index")[:] = Int32.(collect(1:length(layout.js)))
    _define_nc_int_variable(ds, "point_j", ("point",); long_name="1-based y-index for each compact valid cell")[:] = Int32.(layout.js)
    _define_nc_int_variable(ds, "point_i", ("point",); long_name="1-based x-index for each compact valid cell")[:] = Int32.(layout.is)
    _define_nc_double_variable(ds, "point_y_km", ("point",); attrib=Dict("long_name" => "Y coordinate for each compact valid cell", "units" => "km"))[:] = layout.y[layout.js]
    _define_nc_double_variable(ds, "point_x_km", ("point",); attrib=Dict("long_name" => "X coordinate for each compact valid cell", "units" => "km"))[:] = layout.x[layout.is]
    _define_nc_int_variable(ds, "step", ("step",); long_name="Sequential yearly output index across repeated annual cycles")[:] = Int32.(collect(1:max_steps))
    _define_nc_int_variable(ds, "step_cycle", ("step",); long_name="Repeated annual forcing cycle index for each yearly output")[:] = step_cycle
    _define_nc_int_variable(ds, "step_source_index", ("step",); long_name="1-based index of the last forcing step included in each yearly output")[:] = step_source_index
    _define_nc_int_variable(ds, "step_source_code", ("step",); long_name="Source forcing timestamp code YYYYMMDDHH for the final step included in each yearly output")[:] = step_source_code
    _define_nc_output_variable(ds, "domain_mask", ("y", "x"); long_name="Domain mask", units="1")[:, :] = Float32.(layout.mask)
    _define_nc_output_variable(ds, "initial_thickness", ("y", "x"); long_name="Initial snow thickness", units="m")[:, :] = Float32.(initial_thickness)

    vars = _define_selected_nc_variables!(ds, selected)
    if any(key -> key in selected, CASE_OUTPUT_GROUPS.step)
        vars[:step_valid] = _define_nc_int_variable(ds, "step_valid", ("step",); long_name="1 where a yearly output record was completed and written, 0 for unused trailing slots")
    end

    ds.attrib["title"] = options.name
    ds.attrib["source_model"] = "Chion"
    ds.attrib["input_label"] = isempty(options.input_label) ? "not provided" : options.input_label
    ds.attrib["forcing_start"] = string(first(time_values))
    ds.attrib["forcing_end"] = string(last(time_values))
    ds.attrib["cycles_completed"] = "pending"
    ds.attrib["status"] = "pending"
    ds.attrib["created"] = string(now())
    return CaseNetCDFWriter(ds, vars, max_steps, options.cycles)
end

function maybe_write_step_output!(writer::CaseNetCDFWriter, step_index::Int, key::Symbol, data::AbstractMatrix{<:Real})
    haskey(writer.vars, key) && (writer.vars[key][step_index, :, :] = Float32.(data))
    return
end

function _write_dataset_var!(var, data::Vector{Float64})
    var[:] = Float32.(data)
    return
end

function _write_dataset_var!(var, data::Vector{Int32})
    var[:] = data
    return
end

function _write_dataset_var!(var, data::AbstractMatrix{Int32})
    var[:, :] = data
    return
end

function _write_dataset_var!(var, data::AbstractMatrix{<:Real})
    var[:, :] = Float32.(data)
    return
end

function _write_dataset_var!(var, data::Array{Float64, 3})
    var[:, :, :] = Float32.(data)
    return
end

function maybe_write_output!(writer::CaseNetCDFWriter, key::Symbol, data)
    haskey(writer.vars, key) && _write_dataset_var!(writer.vars[key], data)
    return
end

function _history_vectors(history::Vector{NamedTuple}, cycles::Int)
    out = Dict(spec.output => fill(NaN, cycles) for spec in HISTORY_OUTPUT_SPECS)
    for rec in history
        idx = rec.cycle
        for spec in HISTORY_OUTPUT_SPECS
            out[spec.output][idx] = getfield(rec, spec.record)
        end
    end
    return out
end

function finalize_case_netcdf!(
    writer::CaseNetCDFWriter,
    final_grids,
    layer_grids,
    history::Vector{NamedTuple},
    monthly_grids,
    status::Symbol,
    cycles_completed::Int,
    steps_written::Int,
)
    for key in FINAL_GRID_KEYS
        maybe_write_output!(writer, key, getfield(final_grids, key))
    end
    for key in LAYER_GRID_KEYS
        maybe_write_output!(writer, key, getfield(layer_grids, key))
    end
    for (key, values) in _history_vectors(history, writer.cycles)
        maybe_write_output!(writer, key, values)
    end
    for key in MONTHLY_GRID_KEYS
        maybe_write_output!(writer, key, getfield(monthly_grids, key))
    end
    if haskey(writer.vars, :step_valid)
        step_valid = zeros(Int32, writer.max_steps)
        step_valid[1:steps_written] .= 1
        writer.vars[:step_valid][:] = step_valid
    end
    writer.dataset.attrib["cycles_completed"] = string(cycles_completed)
    writer.dataset.attrib["status"] = string(status)
    writer.dataset.attrib["steps_written"] = string(steps_written)
    close(writer.dataset)
    return
end
