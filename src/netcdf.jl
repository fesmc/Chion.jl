@inline _nc_attrib(long_name::AbstractString, units::AbstractString="") =
    units == "" ? Dict("long_name" => String(long_name)) : Dict("long_name" => String(long_name), "units" => String(units))

@inline function _def_nc_var(ds::NCDataset, name::AbstractString, dims, spec)
    return spec.integer ?
        defVar(ds, String(name), Int32, dims; attrib=_nc_attrib(spec.long_name)) :
        defVar(ds, String(name), Float32, dims; fillvalue=NaN32, attrib=_nc_attrib(spec.long_name, spec.units))
end

@inline _write_nc_var!(ds::NCDataset, name, dims, spec, value) = (_def_nc_var(ds, name, dims, spec)[:] = value)

const NC_SPECS = (
    (key=:final_thickness,        name="final_thickness",        dims=("x", "y"),       long_name="Final snow thickness", units="m", integer=false),
    (key=:final_wet_mass,         name="final_wet_mass",         dims=("x", "y"),       long_name="Final snow wet mass", units="mmWE", integer=false),
    (key=:final_bulk_density,     name="final_bulk_density",     dims=("x", "y"),       long_name="Final bulk snow density", units="kg m-3", integer=false),
    (key=:final_base_mass,        name="final_base_mass",        dims=("x", "y"),       long_name="Cumulative firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:final_ice_sheet_smb,    name="final_ice_sheet_smb",    dims=("x", "y"),       long_name="Cumulative net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:final_runoff,           name="final_runoff",           dims=("x", "y"),       long_name="Final cumulative runoff", units="mmWE", integer=false),
    (key=:last_cycle_delta_thickness, name="last_cycle_delta_thickness", dims=("x", "y"), long_name="Last cycle snow-thickness change", units="m", integer=false),
    (key=:last_cycle_delta_wet_mass,  name="last_cycle_delta_wet_mass",  dims=("x", "y"), long_name="Last cycle wet-mass change", units="mmWE", integer=false),
    (key=:last_cycle_delta_base_mass, name="last_cycle_delta_base_mass", dims=("x", "y"), long_name="Last cycle firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:last_cycle_delta_ice_sheet_smb, name="last_cycle_delta_ice_sheet_smb", dims=("x", "y"), long_name="Last cycle net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:n_active,               name="n_active",               dims=("x", "y"),       long_name="Number of active Chion layers", units="", integer=true),
    (key=:layer_density,         name="layer_density",         dims=("layer", "x", "y"), long_name="Final Chion layer density", units="kg m-3", integer=false),
    (key=:layer_thickness,       name="layer_thickness",       dims=("layer", "x", "y"), long_name="Final Chion layer thickness", units="m", integer=false),
    (key=:layer_snow_mass,       name="layer_snow_mass",       dims=("layer", "x", "y"), long_name="Final Chion layer snow mass", units="kg m-2", integer=false),
    (key=:layer_liquid_mass,     name="layer_liquid_mass",     dims=("layer", "x", "y"), long_name="Final Chion layer liquid-water mass", units="kg m-2", integer=false),
    (key=:layer_temperature_c,   name="layer_temperature_c",   dims=("layer", "x", "y"), long_name="Final Chion layer temperature", units="C", integer=false),
    (key=:history_mean_thickness,      name="history_mean_thickness",      dims=("cycle",), long_name="Cycle-mean snow thickness", units="m", integer=false),
    (key=:history_mean_wet_mass,       name="history_mean_wet_mass",       dims=("cycle",), long_name="Cycle-mean snow wet mass", units="mmWE", integer=false),
    (key=:history_mean_bulk_density,   name="history_mean_bulk_density",   dims=("cycle",), long_name="Cycle-mean bulk snow density", units="kg m-3", integer=false),
    (key=:history_mean_base_mass,      name="history_mean_base_mass",      dims=("cycle",), long_name="Cycle-mean firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:history_mean_abs_delta_thickness, name="history_mean_abs_delta_thickness", dims=("cycle",), long_name="Cycle mean absolute snow-thickness change", units="m", integer=false),
    (key=:history_mean_abs_delta_wet_mass,  name="history_mean_abs_delta_wet_mass",  dims=("cycle",), long_name="Cycle mean absolute wet-mass change", units="mmWE", integer=false),
    (key=:history_mean_abs_delta_base_mass, name="history_mean_abs_delta_base_mass", dims=("cycle",), long_name="Cycle mean absolute firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:monthly_mean_thickness,      name="monthly_mean_thickness",      dims=("month", "x", "y"), long_name="Monthly mean snow thickness", units="m", integer=false),
    (key=:monthly_mean_wet_mass,       name="monthly_mean_wet_mass",       dims=("month", "x", "y"), long_name="Monthly mean snow wet mass", units="mmWE", integer=false),
    (key=:monthly_mean_bulk_density,   name="monthly_mean_bulk_density",   dims=("month", "x", "y"), long_name="Monthly mean bulk snow density", units="kg m-3", integer=false),
    (key=:monthly_mean_base_mass,      name="monthly_mean_base_mass",      dims=("month", "x", "y"), long_name="Monthly mean cumulative firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:monthly_mean_ice_sheet_smb,  name="monthly_mean_ice_sheet_smb",  dims=("month", "x", "y"), long_name="Monthly net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:monthly_export_to_ice,       name="monthly_export_to_ice",       dims=("month", "x", "y"), long_name="Monthly firn mass exported to the ice model", units="mmWE", integer=false),
    (key=:monthly_net_ice_sheet_forcing, name="monthly_net_ice_sheet_forcing", dims=("month", "x", "y"), long_name="Monthly net mass forcing to the ice sheet", units="mmWE", integer=false),
    (key=:monthly_runoff,              name="monthly_runoff",              dims=("month", "x", "y"), long_name="Monthly runoff production", units="mmWE", integer=false),
    (key=:step_export_to_ice,          name="step_export_to_ice",          dims=("step", "x", "y"), long_name="Annual firn mass exported to the ice model for each written output interval", units="mmWE", integer=false),
    (key=:step_ice_sheet_smb,          name="step_ice_sheet_smb",          dims=("step", "x", "y"), long_name="Annual net mass forcing to the ice sheet for each written output interval", units="mmWE", integer=false),
)

const HISTORY_OUTPUT_SPECS = (
    (output=:history_mean_thickness, record=:mean_thickness),
    (output=:history_mean_wet_mass, record=:mean_wet_mass),
    (output=:history_mean_bulk_density, record=:mean_bulk_density),
    (output=:history_mean_base_mass, record=:mean_base_mass),
    (output=:history_mean_abs_delta_thickness, record=:mean_abs_delta_thickness),
    (output=:history_mean_abs_delta_wet_mass, record=:mean_abs_delta_wet_mass),
    (output=:history_mean_abs_delta_base_mass, record=:mean_abs_delta_base_mass),
)

function _define_selected_nc_variables!(ds::NCDataset, selected::Set{Symbol})
    vars = Dict{Symbol, Any}()
    for spec in NC_SPECS
        spec.key in selected || continue
        vars[spec.key] = _def_nc_var(ds, spec.name, spec.dims, spec)
    end
    return vars
end

function init_netcdf(
    netcdf_path::AbstractString,
    options::RunOptions,
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
    for (name, len) in (("x", nx), ("y", ny), ("layer", max(nlayer, 1)), ("cycle", options.cycles), ("month", length(month_cycle)), ("point", length(layout.js)), ("step", max_steps))
        defDim(ds, name, len)
    end

    for spec in (
        (name="x", dims=("x",), meta=(key=:x, long_name="X coordinate", units="km", integer=false), value=layout.x),
        (name="y", dims=("y",), meta=(key=:y, long_name="Y coordinate", units="km", integer=false), value=layout.y),
        (name="layer", dims=("layer",), meta=(key=:layer, long_name="Chion internal layer index from surface downward", units="", integer=true), value=Int32.(collect(1:max(nlayer, 1)))),
        (name="cycle", dims=("cycle",), meta=(key=:cycle, long_name="Repeated annual forcing cycle index", units="", integer=true), value=Int32.(collect(1:options.cycles))),
        (name="month", dims=("month",), meta=(key=:month, long_name="Sequential monthly output index", units="", integer=true), value=Int32.(collect(1:length(month_cycle)))),
        (name="month_cycle", dims=("month",), meta=(key=:month_cycle, long_name="Forcing cycle associated with monthly output", units="", integer=true), value=month_cycle),
        (name="month_of_year", dims=("month",), meta=(key=:month_of_year, long_name="Calendar month of the repeated forcing", units="", integer=true), value=month_of_year),
        (name="source_month_code", dims=("month",), meta=(key=:source_month_code, long_name="Source forcing month code YYYYMM", units="", integer=true), value=source_month_code),
        (name="point", dims=("point",), meta=(key=:point, long_name="Compact valid cell index", units="", integer=true), value=Int32.(collect(1:length(layout.js)))),
        (name="point_j", dims=("point",), meta=(key=:point_j, long_name="1-based y-index for each compact valid cell", units="", integer=true), value=Int32.(layout.js)),
        (name="point_i", dims=("point",), meta=(key=:point_i, long_name="1-based x-index for each compact valid cell", units="", integer=true), value=Int32.(layout.is)),
        (name="point_y_km", dims=("point",), meta=(key=:point_y_km, long_name="Y coordinate for each compact valid cell", units="km", integer=false), value=layout.y[layout.js]),
        (name="point_x_km", dims=("point",), meta=(key=:point_x_km, long_name="X coordinate for each compact valid cell", units="km", integer=false), value=layout.x[layout.is]),
        (name="step", dims=("step",), meta=(key=:step, long_name="Sequential yearly output index across repeated annual cycles", units="", integer=true), value=Int32.(collect(1:max_steps))),
        (name="step_cycle", dims=("step",), meta=(key=:step_cycle, long_name="Repeated annual forcing cycle index for each yearly output", units="", integer=true), value=step_cycle),
        (name="step_source_index", dims=("step",), meta=(key=:step_source_index, long_name="1-based index of the last forcing step included in each yearly output", units="", integer=true), value=step_source_index),
        (name="step_source_code", dims=("step",), meta=(key=:step_source_code, long_name="Source forcing timestamp code YYYYMMDDHH for the final step included in each yearly output", units="", integer=true), value=step_source_code),
        (name="domain_mask", dims=("x", "y"), meta=(key=:domain_mask, long_name="Domain mask", units="1", integer=false), value=Float32.(permutedims(layout.mask, (2, 1)))),
        (name="initial_thickness", dims=("x", "y"), meta=(key=:initial_thickness, long_name="Initial snow thickness", units="m", integer=false), value=Float32.(permutedims(initial_thickness, (2, 1)))),
    )
        _write_nc_var!(ds, spec.name, spec.dims, spec.meta, spec.value)
    end

    vars = _define_selected_nc_variables!(ds, selected)
    any(key -> key in selected, OUTPUT_GROUPS.step) && (vars[:step_valid] = _def_nc_var(ds, "step_valid", ("step",), (key=:step_valid, long_name="1 where a yearly output record was completed and written, 0 for unused trailing slots", units="", integer=true)))

    ds.attrib["title"] = options.name
    ds.attrib["source_model"] = "Chion"
    ds.attrib["input_label"] = isempty(options.input_label) ? "not provided" : options.input_label
    ds.attrib["forcing_start"] = string(first(time_values))
    ds.attrib["forcing_end"] = string(last(time_values))
    ds.attrib["cycles_completed"] = "pending"
    ds.attrib["status"] = "pending"
    ds.attrib["created"] = string(now())
    return NetCDFWriter(ds, vars, max_steps, options.cycles)
end

@inline _write_dataset_var!(var, data::AbstractVector) = (var[:] = eltype(var) <: Integer ? data : Float32.(data))
@inline _write_dataset_var!(var, data::AbstractMatrix) = (var[:, :] = eltype(var) <: Integer ? permutedims(data, (2, 1)) : Float32.(permutedims(data, (2, 1))))
@inline _write_dataset_var!(var, data::Array{Float64, 3}) = (var[:, :, :] = Float32.(permutedims(data, (1, 3, 2))))

maybe_write_step_output!(writer::NetCDFWriter, step_index::Int, key::Symbol, data::AbstractMatrix{<:Real}) =
    (haskey(writer.vars, key) && (writer.vars[key][step_index, :, :] = Float32.(permutedims(data, (2, 1))); nothing))

maybe_write_output!(writer::NetCDFWriter, key::Symbol, data) = (haskey(writer.vars, key) && (_write_dataset_var!(writer.vars[key], data); nothing))

function _history_vectors(history::Vector{NamedTuple}, cycles::Int)
    out = Dict(spec.output => fill(NaN, cycles) for spec in HISTORY_OUTPUT_SPECS)
    for rec in history, spec in HISTORY_OUTPUT_SPECS
        out[spec.output][rec.cycle] = getfield(rec, spec.record)
    end
    return out
end

function finalize_netcdf!(
    writer::NetCDFWriter,
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
end
