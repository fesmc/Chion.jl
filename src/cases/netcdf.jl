const NC_NOERR = 0
const NC_CLOBBER = 0x0000
const NC_NETCDF4 = 0x1000
const NC_GLOBAL = -1
const NC_FLOAT = 5
const NC_DOUBLE = 6
const NC_INT = 4

const _LIBNETCDF_CACHE = Ref{Union{Nothing, String}}(nothing)

function _netcdf_search_dirs()
    dirs = String[]
    for key in ("NETCDF_DIR", "NETCDF_HOME", "HOMEBREW_PREFIX")
        root = strip(get(ENV, key, ""))
        isempty(root) || push!(dirs, joinpath(expanduser(root), "lib"))
    end
    for key in ("LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "DYLD_FALLBACK_LIBRARY_PATH")
        value = strip(get(ENV, key, ""))
        isempty(value) && continue
        append!(dirs, filter(!isempty, split(value, ':')))
    end
    append!(dirs, Base.DL_LOAD_PATH)
    append!(dirs, ("/opt/homebrew/lib", "/usr/local/lib", "/opt/local/lib", "/usr/lib", "/lib"))
    return unique(filter(isdir, map(expanduser, dirs)))
end

function _netcdf_library_candidates()
    candidates = String[]
    env_lib = strip(get(ENV, "NETCDF_LIB", ""))
    if !isempty(env_lib)
        push!(candidates, expanduser(env_lib))
    end
    append!(candidates, ("libnetcdf", "libnetcdf.so", "libnetcdf.dylib"))
    for dir in _netcdf_search_dirs()
        append!(candidates, (
            joinpath(dir, "libnetcdf"),
            joinpath(dir, "libnetcdf.so"),
            joinpath(dir, "libnetcdf.dylib"),
        ))
    end
    return unique(candidates)
end

function resolve_libnetcdf()
    for candidate in _netcdf_library_candidates()
        try
            handle = Libdl.dlopen(candidate)
            Libdl.dlclose(handle)
            return candidate
        catch
        end
    end
    searched_dirs = join(_netcdf_search_dirs(), ", ")
    suggestion = Sys.isapple() && isfile("/opt/homebrew/lib/libnetcdf.dylib") ?
        " Try `export NETCDF_LIB=/opt/homebrew/lib/libnetcdf.dylib`." : ""
    error(
        "Could not load NetCDF library. Set NETCDF_LIB to the shared library path or load a NetCDF module." *
        (isempty(searched_dirs) ? "" : " Searched: " * searched_dirs * ".") *
        suggestion,
    )
end

function _libnetcdf()
    if isnothing(_LIBNETCDF_CACHE[])
        _LIBNETCDF_CACHE[] = resolve_libnetcdf()
    end
    return _LIBNETCDF_CACHE[]::String
end

@inline function nc_check(code::Integer)
    if code != NC_NOERR
        msg = unsafe_string(ccall((:nc_strerror, _libnetcdf()), Cstring, (Cint,), code))
        error("NetCDF error: $msg")
    end
    return
end

function nc_create(path::AbstractString)
    isdir(path) && error("NetCDF output path '$(abspath(path))' is a directory; pass a file path ending in `.nc`.")
    ncid = Ref{Cint}()
    code = ccall((:nc_create, _libnetcdf()), Cint, (Cstring, Cint, Ref{Cint}), path, NC_CLOBBER | NC_NETCDF4, ncid)
    if code != NC_NOERR
        msg = unsafe_string(ccall((:nc_strerror, _libnetcdf()), Cstring, (Cint,), code))
        error("NetCDF error while creating '$(abspath(path))': $msg")
    end
    return ncid[]
end

function nc_close(ncid::Cint)
    nc_check(ccall((:nc_close, _libnetcdf()), Cint, (Cint,), ncid))
    return
end

function nc_enddef(ncid::Cint)
    nc_check(ccall((:nc_enddef, _libnetcdf()), Cint, (Cint,), ncid))
    return
end

function nc_redef(ncid::Cint)
    nc_check(ccall((:nc_redef, _libnetcdf()), Cint, (Cint,), ncid))
    return
end

function nc_def_dim(ncid::Cint, name::AbstractString, len::Integer)
    dimid = Ref{Cint}()
    nc_check(ccall((:nc_def_dim, _libnetcdf()), Cint, (Cint, Cstring, Csize_t, Ref{Cint}), ncid, name, len, dimid))
    return dimid[]
end

function nc_def_var(ncid::Cint, name::AbstractString, xtype::Integer, dimids::Vector{Cint})
    varid = Ref{Cint}()
    nc_check(ccall((:nc_def_var, _libnetcdf()), Cint, (Cint, Cstring, Cint, Cint, Ptr{Cint}, Ref{Cint}), ncid, name, xtype, length(dimids), dimids, varid))
    return varid[]
end

function nc_put_att_text(ncid::Cint, varid::Integer, name::AbstractString, value::AbstractString)
    nc_check(ccall((:nc_put_att_text, _libnetcdf()), Cint, (Cint, Cint, Cstring, Csize_t, Cstring), ncid, Cint(varid), name, sizeof(value), value))
    return
end

function nc_put_att_float(ncid::Cint, varid::Integer, name::AbstractString, value::Float32)
    buf = Ref{Float32}(value)
    nc_check(ccall((:nc_put_att_float, _libnetcdf()), Cint, (Cint, Cint, Cstring, Cint, Csize_t, Ref{Float32}), ncid, Cint(varid), name, NC_FLOAT, 1, buf))
    return
end

function nc_put_var_double(ncid::Cint, varid::Integer, data::Vector{Float64})
    nc_check(ccall((:nc_put_var_double, _libnetcdf()), Cint, (Cint, Cint, Ptr{Cdouble}), ncid, Cint(varid), data))
    return
end

function nc_put_var_int(ncid::Cint, varid::Integer, data::Vector{Int32})
    nc_check(ccall((:nc_put_var_int, _libnetcdf()), Cint, (Cint, Cint, Ptr{Cint}), ncid, Cint(varid), data))
    return
end

function nc_put_var_int_2d(ncid::Cint, varid::Integer, data::AbstractMatrix{Int32})
    buf = permutedims(data, (2, 1))
    nc_check(ccall((:nc_put_var_int, _libnetcdf()), Cint, (Cint, Cint, Ptr{Cint}), ncid, Cint(varid), buf))
    return
end

function nc_put_var_float_2d(ncid::Cint, varid::Integer, data::AbstractMatrix{<:Real})
    buf = permutedims(Float32.(data), (2, 1))
    nc_check(ccall((:nc_put_var_float, _libnetcdf()), Cint, (Cint, Cint, Ptr{Cfloat}), ncid, Cint(varid), buf))
    return
end

function nc_put_var_float_3d(ncid::Cint, varid::Integer, data::Array{Float64, 3})
    buf = permutedims(Float32.(data), (3, 2, 1))
    nc_check(ccall((:nc_put_var_float, _libnetcdf()), Cint, (Cint, Cint, Ptr{Cfloat}), ncid, Cint(varid), buf))
    return
end

function nc_put_var_float_1d(ncid::Cint, varid::Integer, data::Vector{Float64})
    buf = Float32.(data)
    nc_check(ccall((:nc_put_var_float, _libnetcdf()), Cint, (Cint, Cint, Ptr{Cfloat}), ncid, Cint(varid), buf))
    return
end

function nc_put_vara_float_3d_step_yx(
    ncid::Cint,
    varid::Integer,
    step_index::Integer,
    data::AbstractMatrix{<:Real},
)
    start = Csize_t[Csize_t(step_index - 1), Csize_t(0), Csize_t(0)]
    count = Csize_t[Csize_t(1), Csize_t(size(data, 1)), Csize_t(size(data, 2))]
    buf = permutedims(Float32.(data), (2, 1))
    nc_check(
        ccall(
            (:nc_put_vara_float, _libnetcdf()),
            Cint,
            (Cint, Cint, Ptr{Csize_t}, Ptr{Csize_t}, Ptr{Cfloat}),
            ncid,
            Cint(varid),
            start,
            count,
            buf,
        ),
    )
    return
end

function define_nc_output_variable(ncid::Cint, dimids::Vector{Cint}, name::AbstractString, long_name::AbstractString, units::AbstractString)
    varid = nc_def_var(ncid, name, NC_FLOAT, dimids)
    nc_put_att_text(ncid, varid, "long_name", long_name)
    nc_put_att_text(ncid, varid, "units", units)
    nc_put_att_float(ncid, varid, "_FillValue", Float32(NaN))
    return varid
end

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

function maybe_define_nc_output_variable!(
    vars::Dict{Symbol, Cint},
    selected::Set{Symbol},
    ncid::Cint,
    dims::Vector{Cint},
    key::Symbol,
    name::AbstractString,
    long_name::AbstractString,
    units::AbstractString,
)
    if key in selected
        vars[key] = define_nc_output_variable(ncid, dims, name, long_name, units)
    end
    return
end

function maybe_define_nc_int_variable!(
    vars::Dict{Symbol, Cint},
    selected::Set{Symbol},
    ncid::Cint,
    dims::Vector{Cint},
    key::Symbol,
    name::AbstractString,
    long_name::AbstractString,
)
    if key in selected
        vars[key] = nc_def_var(ncid, name, NC_INT, dims)
        nc_put_att_text(ncid, vars[key], "long_name", long_name)
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
    isfile(netcdf_path) && rm(netcdf_path, force=true)
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

    ncid = nc_create(netcdf_path)
    dim_y = nc_def_dim(ncid, "y", ny)
    dim_x = nc_def_dim(ncid, "x", nx)
    dim_layer = nc_def_dim(ncid, "layer", max(nlayer, 1))
    dim_cycle = nc_def_dim(ncid, "cycle", options.cycles)
    dim_month = nc_def_dim(ncid, "month", length(month_cycle))
    dim_point = nc_def_dim(ncid, "point", length(layout.js))
    dim_step = nc_def_dim(ncid, "step", max_steps)

    dims_yx = Cint[dim_y, dim_x]
    dims_lyx = Cint[dim_layer, dim_y, dim_x]
    dims_c = Cint[dim_cycle]
    dims_myx = Cint[dim_month, dim_y, dim_x]
    dims_p = Cint[dim_point]
    dims_syx = Cint[dim_step, dim_y, dim_x]

    var_x = nc_def_var(ncid, "x", NC_DOUBLE, Cint[dim_x])
    nc_put_att_text(ncid, var_x, "units", "km")
    nc_put_att_text(ncid, var_x, "axis", "X")
    var_y = nc_def_var(ncid, "y", NC_DOUBLE, Cint[dim_y])
    nc_put_att_text(ncid, var_y, "units", "km")
    nc_put_att_text(ncid, var_y, "axis", "Y")
    var_layer = nc_def_var(ncid, "layer", NC_INT, Cint[dim_layer])
    nc_put_att_text(ncid, var_layer, "long_name", "Chion internal layer index from surface downward")
    var_cycle = nc_def_var(ncid, "cycle", NC_INT, dims_c)
    nc_put_att_text(ncid, var_cycle, "long_name", "Repeated annual forcing cycle index")
    var_month = nc_def_var(ncid, "month", NC_INT, Cint[dim_month])
    nc_put_att_text(ncid, var_month, "long_name", "Sequential monthly output index")
    var_month_cycle = nc_def_var(ncid, "month_cycle", NC_INT, Cint[dim_month])
    nc_put_att_text(ncid, var_month_cycle, "long_name", "Forcing cycle associated with monthly output")
    var_month_of_year = nc_def_var(ncid, "month_of_year", NC_INT, Cint[dim_month])
    nc_put_att_text(ncid, var_month_of_year, "long_name", "Calendar month of the repeated forcing")
    var_source_month_code = nc_def_var(ncid, "source_month_code", NC_INT, Cint[dim_month])
    nc_put_att_text(ncid, var_source_month_code, "long_name", "Source forcing month code YYYYMM")
    var_point = nc_def_var(ncid, "point", NC_INT, dims_p)
    nc_put_att_text(ncid, var_point, "long_name", "Compact valid cell index")
    var_point_j = nc_def_var(ncid, "point_j", NC_INT, dims_p)
    nc_put_att_text(ncid, var_point_j, "long_name", "1-based y-index for each compact valid cell")
    var_point_i = nc_def_var(ncid, "point_i", NC_INT, dims_p)
    nc_put_att_text(ncid, var_point_i, "long_name", "1-based x-index for each compact valid cell")
    var_point_y = nc_def_var(ncid, "point_y_km", NC_DOUBLE, dims_p)
    nc_put_att_text(ncid, var_point_y, "long_name", "Y coordinate for each compact valid cell")
    nc_put_att_text(ncid, var_point_y, "units", "km")
    var_point_x = nc_def_var(ncid, "point_x_km", NC_DOUBLE, dims_p)
    nc_put_att_text(ncid, var_point_x, "long_name", "X coordinate for each compact valid cell")
    nc_put_att_text(ncid, var_point_x, "units", "km")
    var_step = nc_def_var(ncid, "step", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step, "long_name", "Sequential yearly output index across repeated annual cycles")
    var_step_cycle = nc_def_var(ncid, "step_cycle", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step_cycle, "long_name", "Repeated annual forcing cycle index for each yearly output")
    var_step_source_index = nc_def_var(ncid, "step_source_index", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step_source_index, "long_name", "1-based index of the last forcing step included in each yearly output")
    var_step_source_code = nc_def_var(ncid, "step_source_code", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step_source_code, "long_name", "Source forcing timestamp code YYYYMMDDHH for the final step included in each yearly output")

    vars = Dict{Symbol, Cint}()
    var_mask = define_nc_output_variable(ncid, dims_yx, "domain_mask", "Domain mask", "1")
    var_init_th = define_nc_output_variable(ncid, dims_yx, "initial_thickness", "Initial snow thickness", "m")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :final_thickness, "final_thickness", "Final snow thickness", "m")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :final_wet_mass, "final_wet_mass", "Final snow wet mass", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :final_bulk_density, "final_bulk_density", "Final bulk snow density", "kg m-3")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :final_base_mass, "final_base_mass", "Cumulative firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :final_ice_sheet_smb, "final_ice_sheet_smb", "Cumulative net mass forcing to the ice sheet", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :final_runoff, "final_runoff", "Final cumulative runoff", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :last_cycle_delta_thickness, "last_cycle_delta_thickness", "Last cycle snow-thickness change", "m")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :last_cycle_delta_wet_mass, "last_cycle_delta_wet_mass", "Last cycle wet-mass change", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :last_cycle_delta_base_mass, "last_cycle_delta_base_mass", "Last cycle firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :last_cycle_delta_ice_sheet_smb, "last_cycle_delta_ice_sheet_smb", "Last cycle net mass forcing to the ice sheet", "mmWE")
    maybe_define_nc_int_variable!(vars, selected, ncid, dims_yx, :n_active, "n_active", "Number of active Chion layers")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_lyx, :layer_density, "layer_density", "Final Chion layer density", "kg m-3")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_lyx, :layer_thickness, "layer_thickness", "Final Chion layer thickness", "m")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_lyx, :layer_snow_mass, "layer_snow_mass", "Final Chion layer snow mass", "kg m-2")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_lyx, :layer_liquid_mass, "layer_liquid_mass", "Final Chion layer liquid-water mass", "kg m-2")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_lyx, :layer_temperature_c, "layer_temperature_c", "Final Chion layer temperature", "C")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_c, :history_mean_thickness, "history_mean_thickness", "Cycle-mean snow thickness", "m")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_c, :history_mean_wet_mass, "history_mean_wet_mass", "Cycle-mean snow wet mass", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_c, :history_mean_bulk_density, "history_mean_bulk_density", "Cycle-mean bulk snow density", "kg m-3")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_c, :history_mean_base_mass, "history_mean_base_mass", "Cycle-mean firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_c, :history_mean_abs_delta_thickness, "history_mean_abs_delta_thickness", "Cycle mean absolute snow-thickness change", "m")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_c, :history_mean_abs_delta_wet_mass, "history_mean_abs_delta_wet_mass", "Cycle mean absolute wet-mass change", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_c, :history_mean_abs_delta_base_mass, "history_mean_abs_delta_base_mass", "Cycle mean absolute firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_myx, :monthly_mean_thickness, "monthly_mean_thickness", "Monthly mean snow thickness", "m")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_myx, :monthly_mean_wet_mass, "monthly_mean_wet_mass", "Monthly mean snow wet mass", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_myx, :monthly_mean_bulk_density, "monthly_mean_bulk_density", "Monthly mean bulk snow density", "kg m-3")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_myx, :monthly_mean_base_mass, "monthly_mean_base_mass", "Monthly mean cumulative firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_myx, :monthly_mean_ice_sheet_smb, "monthly_mean_ice_sheet_smb", "Monthly net mass forcing to the ice sheet", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_myx, :monthly_export_to_ice, "monthly_export_to_ice", "Monthly firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_myx, :monthly_net_ice_sheet_forcing, "monthly_net_ice_sheet_forcing", "Monthly net mass forcing to the ice sheet", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_myx, :monthly_runoff, "monthly_runoff", "Monthly runoff production", "mmWE")
    if any(v -> v in selected, (:step_export_to_ice, :step_ice_sheet_smb))
        vars[:step_valid] = nc_def_var(ncid, "step_valid", NC_INT, Cint[dim_step])
        nc_put_att_text(ncid, vars[:step_valid], "long_name", "1 where a yearly output record was completed and written, 0 for unused trailing slots")
    end
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_syx, :step_export_to_ice, "step_export_to_ice", "Annual firn mass exported to the ice model for each written output interval", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_syx, :step_ice_sheet_smb, "step_ice_sheet_smb", "Annual net mass forcing to the ice sheet for each written output interval", "mmWE")

    nc_put_att_text(ncid, NC_GLOBAL, "title", options.name)
    nc_put_att_text(ncid, NC_GLOBAL, "source_model", "Chion")
    nc_put_att_text(ncid, NC_GLOBAL, "input_label", isempty(options.input_label) ? "not provided" : options.input_label)
    nc_put_att_text(ncid, NC_GLOBAL, "forcing_start", string(first(time_values)))
    nc_put_att_text(ncid, NC_GLOBAL, "forcing_end", string(last(time_values)))
    nc_put_att_text(ncid, NC_GLOBAL, "cycles_completed", "pending")
    nc_put_att_text(ncid, NC_GLOBAL, "status", "pending")
    nc_put_att_text(ncid, NC_GLOBAL, "created", string(now()))

    nc_enddef(ncid)
    nc_put_var_double(ncid, var_x, layout.x)
    nc_put_var_double(ncid, var_y, layout.y)
    nc_put_var_int(ncid, var_layer, Int32.(collect(1:max(nlayer, 1))))
    nc_put_var_int(ncid, var_cycle, Int32.(collect(1:options.cycles)))
    nc_put_var_int(ncid, var_month, Int32.(collect(1:length(month_cycle))))
    nc_put_var_int(ncid, var_month_cycle, month_cycle)
    nc_put_var_int(ncid, var_month_of_year, month_of_year)
    nc_put_var_int(ncid, var_source_month_code, source_month_code)
    nc_put_var_int(ncid, var_point, Int32.(collect(1:length(layout.js))))
    nc_put_var_int(ncid, var_point_j, Int32.(layout.js))
    nc_put_var_int(ncid, var_point_i, Int32.(layout.is))
    nc_put_var_double(ncid, var_point_y, layout.y[layout.js])
    nc_put_var_double(ncid, var_point_x, layout.x[layout.is])
    nc_put_var_int(ncid, var_step, Int32.(collect(1:max_steps)))
    nc_put_var_int(ncid, var_step_cycle, step_cycle)
    nc_put_var_int(ncid, var_step_source_index, step_source_index)
    nc_put_var_int(ncid, var_step_source_code, step_source_code)
    nc_put_var_float_2d(ncid, var_mask, layout.mask)
    nc_put_var_float_2d(ncid, var_init_th, initial_thickness)
    return CaseNetCDFWriter(ncid, vars, max_steps, options.cycles)
end

function maybe_write_step_output!(writer::CaseNetCDFWriter, step_index::Int, key::Symbol, data::AbstractMatrix{<:Real})
    if haskey(writer.vars, key)
        nc_put_vara_float_3d_step_yx(writer.ncid, writer.vars[key], step_index, data)
    end
    return
end

function maybe_write_output!(writer::CaseNetCDFWriter, key::Symbol, data, writer_fn)
    if haskey(writer.vars, key)
        writer_fn(writer.ncid, writer.vars[key], data)
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
    maybe_write_output!(writer, :final_thickness, final_thickness, nc_put_var_float_2d)
    maybe_write_output!(writer, :final_wet_mass, final_wet_mass, nc_put_var_float_2d)
    maybe_write_output!(writer, :final_bulk_density, final_bulk_density, nc_put_var_float_2d)
    maybe_write_output!(writer, :final_base_mass, final_base_mass, nc_put_var_float_2d)
    maybe_write_output!(writer, :final_ice_sheet_smb, final_ice_sheet_smb, nc_put_var_float_2d)
    maybe_write_output!(writer, :final_runoff, final_runoff, nc_put_var_float_2d)
    maybe_write_output!(writer, :last_cycle_delta_thickness, last_delta_thickness, nc_put_var_float_2d)
    maybe_write_output!(writer, :last_cycle_delta_wet_mass, last_delta_wet_mass, nc_put_var_float_2d)
    maybe_write_output!(writer, :last_cycle_delta_base_mass, last_delta_base_mass, nc_put_var_float_2d)
    maybe_write_output!(writer, :last_cycle_delta_ice_sheet_smb, last_delta_ice_sheet_smb, nc_put_var_float_2d)
    maybe_write_output!(writer, :n_active, layer_grids.n_active, nc_put_var_int_2d)
    maybe_write_output!(writer, :layer_density, layer_grids.layer_density, nc_put_var_float_3d)
    maybe_write_output!(writer, :layer_thickness, layer_grids.layer_thickness, nc_put_var_float_3d)
    maybe_write_output!(writer, :layer_snow_mass, layer_grids.layer_snow_mass, nc_put_var_float_3d)
    maybe_write_output!(writer, :layer_liquid_mass, layer_grids.layer_liquid_mass, nc_put_var_float_3d)
    maybe_write_output!(writer, :layer_temperature_c, layer_grids.layer_temperature_c, nc_put_var_float_3d)
    maybe_write_output!(writer, :history_mean_thickness, hist_th, nc_put_var_float_1d)
    maybe_write_output!(writer, :history_mean_wet_mass, hist_wet, nc_put_var_float_1d)
    maybe_write_output!(writer, :history_mean_bulk_density, hist_rho, nc_put_var_float_1d)
    maybe_write_output!(writer, :history_mean_base_mass, hist_base, nc_put_var_float_1d)
    maybe_write_output!(writer, :history_mean_abs_delta_thickness, hist_dth, nc_put_var_float_1d)
    maybe_write_output!(writer, :history_mean_abs_delta_wet_mass, hist_dswe, nc_put_var_float_1d)
    maybe_write_output!(writer, :history_mean_abs_delta_base_mass, hist_dbase, nc_put_var_float_1d)
    maybe_write_output!(writer, :monthly_mean_thickness, monthly_mean_thickness, nc_put_var_float_3d)
    maybe_write_output!(writer, :monthly_mean_wet_mass, monthly_mean_wet_mass, nc_put_var_float_3d)
    maybe_write_output!(writer, :monthly_mean_bulk_density, monthly_mean_bulk_density, nc_put_var_float_3d)
    maybe_write_output!(writer, :monthly_mean_base_mass, monthly_mean_base_mass, nc_put_var_float_3d)
    maybe_write_output!(writer, :monthly_mean_ice_sheet_smb, monthly_mean_ice_sheet_smb, nc_put_var_float_3d)
    maybe_write_output!(writer, :monthly_export_to_ice, monthly_export_to_ice, nc_put_var_float_3d)
    maybe_write_output!(writer, :monthly_net_ice_sheet_forcing, monthly_net_ice_sheet_forcing, nc_put_var_float_3d)
    maybe_write_output!(writer, :monthly_runoff, monthly_runoff, nc_put_var_float_3d)
    if haskey(writer.vars, :step_valid)
        step_valid = zeros(Int32, writer.max_steps)
        step_valid[1:steps_written] .= 1
        nc_put_var_int(writer.ncid, writer.vars[:step_valid], step_valid)
    end
    nc_redef(writer.ncid)
    nc_put_att_text(writer.ncid, NC_GLOBAL, "cycles_completed", string(cycles_completed))
    nc_put_att_text(writer.ncid, NC_GLOBAL, "status", string(status))
    nc_put_att_text(writer.ncid, NC_GLOBAL, "steps_written", string(steps_written))
    nc_enddef(writer.ncid)
    nc_close(writer.ncid)
    return
end
