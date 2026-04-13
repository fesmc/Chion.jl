
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

@inline _nc_float_attrib(long_name::AbstractString, units::AbstractString) =
    Dict("long_name" => String(long_name), "units" => String(units))

@inline _nc_text_attrib(long_name::AbstractString) =
    Dict("long_name" => String(long_name))

function _define_nc_output_variable(ds::NCDataset, name::AbstractString, dims; long_name::AbstractString, units::AbstractString)
    return defVar(
        ds,
        String(name),
        Float32,
        dims;
        fillvalue=NaN32,
        attrib=_nc_float_attrib(long_name, units),
    )
end

function _define_nc_int_variable(ds::NCDataset, name::AbstractString, dims; long_name::AbstractString)
    return defVar(ds, String(name), Int32, dims; attrib=_nc_text_attrib(long_name))
end

function _define_nc_double_variable(ds::NCDataset, name::AbstractString, dims; attrib::Dict{String, String}=Dict{String, String}())
    return defVar(ds, String(name), Float64, dims; attrib=attrib)
end

function maybe_define_nc_output_variable!(
    vars::Dict{Symbol, Any},
    selected::Set{Symbol},
    ds::NCDataset,
    dims,
    key::Symbol,
    name::AbstractString,
    long_name::AbstractString,
    units::AbstractString,
)
    if key in selected
        vars[key] = _define_nc_output_variable(ds, name, dims; long_name=long_name, units=units)
    end
    return
end

function maybe_define_nc_int_variable!(
    vars::Dict{Symbol, Any},
    selected::Set{Symbol},
    ds::NCDataset,
    dims,
    key::Symbol,
    name::AbstractString,
    long_name::AbstractString,
)
    if key in selected
        vars[key] = _define_nc_int_variable(ds, name, dims; long_name=long_name)
    end
    return
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
    step_cycle = Vector{Int32}(undef, max_steps)
    step_source_index = Vector{Int32}(undef, max_steps)
    step_source_code = Vector{Int32}(undef, max_steps)
    step_counter = 0
    for cyc in 1:options.cycles
        for annual_idx in eachindex(annual_output_source_indices)
            step_counter += 1
            step_cycle[step_counter] = Int32(cyc)
            step_source_index[step_counter] = annual_output_source_indices[annual_idx]
            step_source_code[step_counter] = annual_output_source_codes[annual_idx]
        end
    end

    ds = NCDataset(netcdf_path, "c")
    defDim(ds, "y", ny)
    defDim(ds, "x", nx)
    defDim(ds, "layer", max(nlayer, 1))
    defDim(ds, "cycle", options.cycles)
    defDim(ds, "month", length(month_cycle))
    defDim(ds, "point", length(layout.js))
    defDim(ds, "step", max_steps)

    var_x = _define_nc_double_variable(ds, "x", ("x",); attrib=Dict("units" => "km", "axis" => "X"))
    var_y = _define_nc_double_variable(ds, "y", ("y",); attrib=Dict("units" => "km", "axis" => "Y"))
    var_layer = _define_nc_int_variable(ds, "layer", ("layer",); long_name="Chion internal layer index from surface downward")
    var_cycle = _define_nc_int_variable(ds, "cycle", ("cycle",); long_name="Repeated annual forcing cycle index")
    var_month = _define_nc_int_variable(ds, "month", ("month",); long_name="Sequential monthly output index")
    var_month_cycle = _define_nc_int_variable(ds, "month_cycle", ("month",); long_name="Forcing cycle associated with monthly output")
    var_month_of_year = _define_nc_int_variable(ds, "month_of_year", ("month",); long_name="Calendar month of the repeated forcing")
    var_source_month_code = _define_nc_int_variable(ds, "source_month_code", ("month",); long_name="Source forcing month code YYYYMM")
    var_point = _define_nc_int_variable(ds, "point", ("point",); long_name="Compact valid cell index")
    var_point_j = _define_nc_int_variable(ds, "point_j", ("point",); long_name="1-based y-index for each compact valid cell")
    var_point_i = _define_nc_int_variable(ds, "point_i", ("point",); long_name="1-based x-index for each compact valid cell")
    var_point_y = _define_nc_double_variable(ds, "point_y_km", ("point",); attrib=Dict("long_name" => "Y coordinate for each compact valid cell", "units" => "km"))
    var_point_x = _define_nc_double_variable(ds, "point_x_km", ("point",); attrib=Dict("long_name" => "X coordinate for each compact valid cell", "units" => "km"))
    var_step = _define_nc_int_variable(ds, "step", ("step",); long_name="Sequential yearly output index across repeated annual cycles")
    var_step_cycle = _define_nc_int_variable(ds, "step_cycle", ("step",); long_name="Repeated annual forcing cycle index for each yearly output")
    var_step_source_index = _define_nc_int_variable(ds, "step_source_index", ("step",); long_name="1-based index of the last forcing step included in each yearly output")
    var_step_source_code = _define_nc_int_variable(ds, "step_source_code", ("step",); long_name="Source forcing timestamp code YYYYMMDDHH for the final step included in each yearly output")

    vars = Dict{Symbol, Any}()
    var_mask = _define_nc_output_variable(ds, "domain_mask", ("y", "x"); long_name="Domain mask", units="1")
    var_init_th = _define_nc_output_variable(ds, "initial_thickness", ("y", "x"); long_name="Initial snow thickness", units="m")
    maybe_define_nc_output_variable!(vars, selected, ds, ("y", "x"), :final_thickness, "final_thickness", "Final snow thickness", "m")
    maybe_define_nc_output_variable!(vars, selected, ds, ("y", "x"), :final_wet_mass, "final_wet_mass", "Final snow wet mass", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ds, ("y", "x"), :final_bulk_density, "final_bulk_density", "Final bulk snow density", "kg m-3")
    maybe_define_nc_output_variable!(vars, selected, ds, ("y", "x"), :final_base_mass, "final_base_mass", "Cumulative firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ds, ("y", "x"), :final_ice_sheet_smb, "final_ice_sheet_smb", "Cumulative net mass forcing to the ice sheet", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ds, ("y", "x"), :final_runoff, "final_runoff", "Final cumulative runoff", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ds, ("y", "x"), :last_cycle_delta_thickness, "last_cycle_delta_thickness", "Last cycle snow-thickness change", "m")
    maybe_define_nc_output_variable!(vars, selected, ds, ("y", "x"), :last_cycle_delta_wet_mass, "last_cycle_delta_wet_mass", "Last cycle wet-mass change", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ds, ("y", "x"), :last_cycle_delta_base_mass, "last_cycle_delta_base_mass", "Last cycle firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ds, ("y", "x"), :last_cycle_delta_ice_sheet_smb, "last_cycle_delta_ice_sheet_smb", "Last cycle net mass forcing to the ice sheet", "mmWE")
    maybe_define_nc_int_variable!(vars, selected, ds, ("y", "x"), :n_active, "n_active", "Number of active Chion layers")
    maybe_define_nc_output_variable!(vars, selected, ds, ("layer", "y", "x"), :layer_density, "layer_density", "Final Chion layer density", "kg m-3")
    maybe_define_nc_output_variable!(vars, selected, ds, ("layer", "y", "x"), :layer_thickness, "layer_thickness", "Final Chion layer thickness", "m")
    maybe_define_nc_output_variable!(vars, selected, ds, ("layer", "y", "x"), :layer_snow_mass, "layer_snow_mass", "Final Chion layer snow mass", "kg m-2")
    maybe_define_nc_output_variable!(vars, selected, ds, ("layer", "y", "x"), :layer_liquid_mass, "layer_liquid_mass", "Final Chion layer liquid-water mass", "kg m-2")
    maybe_define_nc_output_variable!(vars, selected, ds, ("layer", "y", "x"), :layer_temperature_c, "layer_temperature_c", "Final Chion layer temperature", "C")
    maybe_define_nc_output_variable!(vars, selected, ds, ("cycle",), :history_mean_thickness, "history_mean_thickness", "Cycle-mean snow thickness", "m")
    maybe_define_nc_output_variable!(vars, selected, ds, ("cycle",), :history_mean_wet_mass, "history_mean_wet_mass", "Cycle-mean snow wet mass", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ds, ("cycle",), :history_mean_bulk_density, "history_mean_bulk_density", "Cycle-mean bulk snow density", "kg m-3")
    maybe_define_nc_output_variable!(vars, selected, ds, ("cycle",), :history_mean_base_mass, "history_mean_base_mass", "Cycle-mean firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ds, ("cycle",), :history_mean_abs_delta_thickness, "history_mean_abs_delta_thickness", "Cycle mean absolute snow-thickness change", "m")
    maybe_define_nc_output_variable!(vars, selected, ds, ("cycle",), :history_mean_abs_delta_wet_mass, "history_mean_abs_delta_wet_mass", "Cycle mean absolute wet-mass change", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ds, ("cycle",), :history_mean_abs_delta_base_mass, "history_mean_abs_delta_base_mass", "Cycle mean absolute firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ds, ("month", "y", "x"), :monthly_mean_thickness, "monthly_mean_thickness", "Monthly mean snow thickness", "m")
    maybe_define_nc_output_variable!(vars, selected, ds, ("month", "y", "x"), :monthly_mean_wet_mass, "monthly_mean_wet_mass", "Monthly mean snow wet mass", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ds, ("month", "y", "x"), :monthly_mean_bulk_density, "monthly_mean_bulk_density", "Monthly mean bulk snow density", "kg m-3")
    maybe_define_nc_output_variable!(vars, selected, ds, ("month", "y", "x"), :monthly_mean_base_mass, "monthly_mean_base_mass", "Monthly mean cumulative firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ds, ("month", "y", "x"), :monthly_mean_ice_sheet_smb, "monthly_mean_ice_sheet_smb", "Monthly net mass forcing to the ice sheet", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ds, ("month", "y", "x"), :monthly_export_to_ice, "monthly_export_to_ice", "Monthly firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ds, ("month", "y", "x"), :monthly_net_ice_sheet_forcing, "monthly_net_ice_sheet_forcing", "Monthly net mass forcing to the ice sheet", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ds, ("month", "y", "x"), :monthly_runoff, "monthly_runoff", "Monthly runoff production", "mmWE")
    if any(v -> v in selected, (:step_export_to_ice, :step_ice_sheet_smb))
        vars[:step_valid] = _define_nc_int_variable(ds, "step_valid", ("step",); long_name="1 where a yearly output record was completed and written, 0 for unused trailing slots")
    end
    maybe_define_nc_output_variable!(vars, selected, ds, ("step", "y", "x"), :step_export_to_ice, "step_export_to_ice", "Annual firn mass exported to the ice model for each written output interval", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ds, ("step", "y", "x"), :step_ice_sheet_smb, "step_ice_sheet_smb", "Annual net mass forcing to the ice sheet for each written output interval", "mmWE")

    ds.attrib["title"] = options.name
    ds.attrib["source_model"] = "Chion"
    ds.attrib["input_label"] = isempty(options.input_label) ? "not provided" : options.input_label
    ds.attrib["forcing_start"] = string(first(time_values))
    ds.attrib["forcing_end"] = string(last(time_values))
    ds.attrib["cycles_completed"] = "pending"
    ds.attrib["status"] = "pending"
    ds.attrib["created"] = string(now())

    var_x[:] = layout.x
    var_y[:] = layout.y
    var_layer[:] = Int32.(collect(1:max(nlayer, 1)))
    var_cycle[:] = Int32.(collect(1:options.cycles))
    var_month[:] = Int32.(collect(1:length(month_cycle)))
    var_month_cycle[:] = month_cycle
    var_month_of_year[:] = month_of_year
    var_source_month_code[:] = source_month_code
    var_point[:] = Int32.(collect(1:length(layout.js)))
    var_point_j[:] = Int32.(layout.js)
    var_point_i[:] = Int32.(layout.is)
    var_point_y[:] = layout.y[layout.js]
    var_point_x[:] = layout.x[layout.is]
    var_step[:] = Int32.(collect(1:max_steps))
    var_step_cycle[:] = step_cycle
    var_step_source_index[:] = step_source_index
    var_step_source_code[:] = step_source_code
    var_mask[:, :] = Float32.(layout.mask)
    var_init_th[:, :] = Float32.(initial_thickness)
    return CaseNetCDFWriter(ds, vars, max_steps, options.cycles)
end

function maybe_write_step_output!(writer::CaseNetCDFWriter, step_index::Int, key::Symbol, data::AbstractMatrix{<:Real})
    if haskey(writer.vars, key)
        writer.vars[key][step_index, :, :] = Float32.(data)
    end
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
    if haskey(writer.vars, key)
        _write_dataset_var!(writer.vars[key], data)
    end
    return
end

function finalize_case_netcdf!(
    writer::CaseNetCDFWriter,
    final_thickness::Matrix{Float64},
    final_wet_mass::Matrix{Float64},
    final_bulk_density::Matrix{Float64},
    final_base_mass::Matrix{Float64},
    final_ice_sheet_smb::Matrix{Float64},
    last_delta_thickness::Matrix{Float64},
    last_delta_wet_mass::Matrix{Float64},
    last_delta_base_mass::Matrix{Float64},
    last_delta_ice_sheet_smb::Matrix{Float64},
    final_runoff::Matrix{Float64},
    layer_grids,
    history::Vector{NamedTuple},
    monthly_mean_thickness::Array{Float64, 3},
    monthly_mean_wet_mass::Array{Float64, 3},
    monthly_mean_bulk_density::Array{Float64, 3},
    monthly_mean_base_mass::Array{Float64, 3},
    monthly_mean_ice_sheet_smb::Array{Float64, 3},
    monthly_export_to_ice::Array{Float64, 3},
    monthly_net_ice_sheet_forcing::Array{Float64, 3},
    monthly_runoff::Array{Float64, 3},
    status::Symbol,
    cycles_completed::Int,
    steps_written::Int,
)
    hist_th = fill(NaN, writer.cycles)
    hist_wet = fill(NaN, writer.cycles)
    hist_rho = fill(NaN, writer.cycles)
    hist_base = fill(NaN, writer.cycles)
    hist_dth = fill(NaN, writer.cycles)
    hist_dswe = fill(NaN, writer.cycles)
    hist_dbase = fill(NaN, writer.cycles)
    for rec in history
        idx = getfield(rec, :cycle)
        hist_th[idx] = getfield(rec, :mean_thickness)
        hist_wet[idx] = getfield(rec, :mean_wet_mass)
        hist_rho[idx] = getfield(rec, :mean_bulk_density)
        hist_base[idx] = getfield(rec, :mean_base_mass)
        hist_dth[idx] = getfield(rec, :mean_abs_delta_thickness)
        hist_dswe[idx] = getfield(rec, :mean_abs_delta_wet_mass)
        hist_dbase[idx] = getfield(rec, :mean_abs_delta_base_mass)
    end
    maybe_write_output!(writer, :final_thickness, final_thickness)
    maybe_write_output!(writer, :final_wet_mass, final_wet_mass)
    maybe_write_output!(writer, :final_bulk_density, final_bulk_density)
    maybe_write_output!(writer, :final_base_mass, final_base_mass)
    maybe_write_output!(writer, :final_ice_sheet_smb, final_ice_sheet_smb)
    maybe_write_output!(writer, :final_runoff, final_runoff)
    maybe_write_output!(writer, :last_cycle_delta_thickness, last_delta_thickness)
    maybe_write_output!(writer, :last_cycle_delta_wet_mass, last_delta_wet_mass)
    maybe_write_output!(writer, :last_cycle_delta_base_mass, last_delta_base_mass)
    maybe_write_output!(writer, :last_cycle_delta_ice_sheet_smb, last_delta_ice_sheet_smb)
    maybe_write_output!(writer, :n_active, layer_grids.n_active)
    maybe_write_output!(writer, :layer_density, layer_grids.layer_density)
    maybe_write_output!(writer, :layer_thickness, layer_grids.layer_thickness)
    maybe_write_output!(writer, :layer_snow_mass, layer_grids.layer_snow_mass)
    maybe_write_output!(writer, :layer_liquid_mass, layer_grids.layer_liquid_mass)
    maybe_write_output!(writer, :layer_temperature_c, layer_grids.layer_temperature_c)
    maybe_write_output!(writer, :history_mean_thickness, hist_th)
    maybe_write_output!(writer, :history_mean_wet_mass, hist_wet)
    maybe_write_output!(writer, :history_mean_bulk_density, hist_rho)
    maybe_write_output!(writer, :history_mean_base_mass, hist_base)
    maybe_write_output!(writer, :history_mean_abs_delta_thickness, hist_dth)
    maybe_write_output!(writer, :history_mean_abs_delta_wet_mass, hist_dswe)
    maybe_write_output!(writer, :history_mean_abs_delta_base_mass, hist_dbase)
    maybe_write_output!(writer, :monthly_mean_thickness, monthly_mean_thickness)
    maybe_write_output!(writer, :monthly_mean_wet_mass, monthly_mean_wet_mass)
    maybe_write_output!(writer, :monthly_mean_bulk_density, monthly_mean_bulk_density)
    maybe_write_output!(writer, :monthly_mean_base_mass, monthly_mean_base_mass)
    maybe_write_output!(writer, :monthly_mean_ice_sheet_smb, monthly_mean_ice_sheet_smb)
    maybe_write_output!(writer, :monthly_export_to_ice, monthly_export_to_ice)
    maybe_write_output!(writer, :monthly_net_ice_sheet_forcing, monthly_net_ice_sheet_forcing)
    maybe_write_output!(writer, :monthly_runoff, monthly_runoff)
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
