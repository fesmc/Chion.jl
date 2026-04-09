module ChionMarProblemLoader

using Dates
using HDF5
using Printf
using Base.Threads
using Chion

const SM = Chion.SnowpackModel
const FILL_THRESHOLD = -9.0e18

export default_gris_nc_path
export read_dataset_shapes, read_hdf5_subset, read_hdf5_full
export read_timeslice_2d, read_timeslice_3d
export valid_or, mmwe_day_to_kgm2s
export read_mar_times, choose_time_index, infer_dt_days
export extract_mar_layers, populate_domain_column_from_mar!
export read_full_timeseries_3d, read_first_available_timeseries_3d
export load_gris_mar_problem

function default_gris_nc_path()
    for candidate in (
        "/p/projects/ou/labs/ai/Nils/MARv3.14.3-10km-daily-ERA5-2025.nc",
        "/Users/niboch001/Downloads/MARv3.14.3-10km-daily-ERA5-2026.nc",
    )
        isfile(candidate) && return candidate
    end
    return ""
end

function read_dataset_shapes(nc_path::AbstractString)
    shapes = Dict{String, Vector{Int}}()
    h5open(nc_path, "r") do file
        for name in keys(file)
            obj = file[name]
            if obj isa HDF5.Dataset
                # Preserve the logical dimension order used throughout the examples:
                # time, optional layer, y, x.
                shapes[String(name)] = reverse(collect(size(obj)))
            end
        end
    end
    return shapes
end

function _clean_fill!(A)
    @inbounds for i in eachindex(A)
        if A[i] <= FILL_THRESHOLD
            A[i] = NaN
        end
    end
    return A
end

function read_hdf5_subset(
    nc_path::AbstractString,
    varname::AbstractString,
    full_shape::Vector{Int};
    start::Union{Nothing, Vector{Int}}=nothing,
    count::Union{Nothing, Vector{Int}}=nothing,
)
    start_vec = isnothing(start) ? zeros(Int, length(full_shape)) : start
    count_vec = isnothing(count) ? copy(full_shape) : count
    length(start_vec) == length(full_shape) || error("start rank mismatch for $varname")
    length(count_vec) == length(full_shape) || error("count rank mismatch for $varname")
    data = h5open(nc_path, "r") do file
        dataset = file[varname]
        nd = length(full_shape)
        file_ranges = ntuple(nd) do dim
            logical_dim = nd - dim + 1
            first_index = start_vec[logical_dim] + 1
            last_index = start_vec[logical_dim] + count_vec[logical_dim]
            first_index:last_index
        end
        raw = dataset[file_ranges...]
        logical = Float64.(raw)
        return nd > 1 ? permutedims(logical, nd:-1:1) : logical
    end
    _clean_fill!(data)
    return data
end

function read_hdf5_full(nc_path::AbstractString, varname::AbstractString, shapes::Dict{String, Vector{Int}})
    haskey(shapes, varname) || error("Variable '$varname' not found in $nc_path.")
    return read_hdf5_subset(nc_path, varname, shapes[varname])
end

function read_timeslice_2d(
    nc_path::AbstractString,
    varname::AbstractString,
    time_index::Int,
    shapes::Dict{String, Vector{Int}},
)
    shape = shapes[varname]
    if length(shape) == 3
        start = [time_index - 1, 0, 0]
        count = [1, shape[2], shape[3]]
    elseif length(shape) == 4
        start = [time_index - 1, 0, 0, 0]
        count = [1, 1, shape[3], shape[4]]
    else
        error("Variable '$varname' does not have a supported rank for 2D slicing.")
    end
    data = read_hdf5_subset(nc_path, varname, shape; start=start, count=count)
    return dropdims(data; dims=Tuple(findall(==(1), size(data))))
end

function read_timeslice_3d(
    nc_path::AbstractString,
    varname::AbstractString,
    time_index::Int,
    shapes::Dict{String, Vector{Int}},
)
    shape = shapes[varname]
    length(shape) == 4 || error("Variable '$varname' does not have the expected 4D layout.")
    start = [time_index - 1, 0, 0, 0]
    count = [1, shape[2], shape[3], shape[4]]
    data = read_hdf5_subset(nc_path, varname, shape; start=start, count=count)
    return dropdims(data; dims=(1,))
end

@inline valid_or(default::Float64, x::Float64) = isfinite(x) ? x : default
@inline mmwe_day_to_kgm2s(x::Float64) = isfinite(x) ? max(x, 0.0) / 86_400.0 : 0.0

function read_mar_times(nc_path::AbstractString, shapes::Dict{String, Vector{Int}})
    yyyy = round.(Int, vec(read_hdf5_full(nc_path, "YYYY", shapes)))
    mm = round.(Int, vec(read_hdf5_full(nc_path, "MM", shapes)))
    dd = round.(Int, vec(read_hdf5_full(nc_path, "DD", shapes)))
    hh = round.(Int, vec(read_hdf5_full(nc_path, "HH", shapes)))
    ntime = length(yyyy)
    codes = Vector{String}(undef, ntime)
    times = Vector{DateTime}(undef, ntime)
    for i in 1:ntime
        codes[i] = @sprintf("%04d%02d%02d%02d", yyyy[i], mm[i], dd[i], hh[i])
        times[i] = DateTime(yyyy[i], mm[i], dd[i], hh[i])
    end
    return codes, times
end

function choose_time_index(config, date_codes::Vector{String})
    ntime = length(date_codes)
    if !isempty(strip(config.date_code))
        wanted = strip(config.date_code)
        for i in eachindex(date_codes)
            if date_codes[i] == wanted
                return i
            end
        end
        error("DATE=$(wanted) is not present in the file.")
    end
    1 <= config.time_index <= ntime || error("--time-index must be between 1 and $ntime.")
    return config.time_index
end

function infer_dt_days(time_values::Vector{DateTime}, time_index::Int)
    if length(time_values) == 1
        return 1.0
    elseif time_index < length(time_values)
        return Dates.value(time_values[time_index + 1] - time_values[time_index]) / (1000 * 60 * 60 * 24)
    else
        return Dates.value(time_values[time_index] - time_values[time_index - 1]) / (1000 * 60 * 60 * 24)
    end
end

function extract_mar_layers(
    total_height::Float64,
    density_profile::AbstractVector{<:Real},
    temperature_profile_c::AbstractVector{<:Real},
    liquid_water_profile::AbstractVector{<:Real},
    outlay_bounds::AbstractMatrix{<:Real};
    ntot::Int=80,
    c::SM.SnowpackPhysicalConstants=SM.SnowpackPhysicalConstants(),
    mass_split::Float64=SM.DEFAULT_MASS_SPLIT,
)
    if !isfinite(total_height) || total_height <= 0.0
        return (
            N=0,
            mass=Float64[],
            mass_w=Float64[],
            density=Float64[],
            temperature=Float64[],
        )
    end

    layer_mass = Float64[]
    layer_mass_w = Float64[]
    layer_density = Float64[]
    layer_temperature = Float64[]

    function append_restart_layer!(
        snow_mass::Float64,
        liquid_mass::Float64,
        density::Float64,
        temperature::Float64,
    )
        snow_mass <= 0.0 && return
        # Preserve the native MAR layering on restart. Pre-splitting into
        # Chion-sized chunks shifts the day-zero vertical grid.
        push!(layer_mass, snow_mass)
        push!(layer_mass_w, liquid_mass)
        push!(layer_density, density)
        push!(layer_temperature, temperature)
        return
    end

    function merge_restart_bottom_pair!()
        n = length(layer_mass)
        n >= 2 || return

        lower_mass = layer_mass[n - 1]
        bottom_mass = layer_mass[n]
        combined_mass = lower_mass + bottom_mass
        if combined_mass <= 0.0
            layer_mass[n - 1] = 0.0
            layer_mass_w[n - 1] += layer_mass_w[n]
            layer_density[n - 1] = max(layer_density[n - 1], layer_density[n])
            layer_temperature[n - 1] = min(layer_temperature[n - 1], layer_temperature[n])
        else
            layer_mass[n - 1] = combined_mass
            layer_mass_w[n - 1] += layer_mass_w[n]
            layer_density[n - 1] =
                (lower_mass * layer_density[n - 1] + bottom_mass * layer_density[n]) / combined_mass
            layer_temperature[n - 1] =
                (lower_mass * layer_temperature[n - 1] + bottom_mass * layer_temperature[n]) / combined_mass
        end

        pop!(layer_mass)
        pop!(layer_mass_w)
        pop!(layer_density)
        pop!(layer_temperature)
        return
    end

    nlayers = size(outlay_bounds, 1)
    for k in 1:nlayers
        lower = max(outlay_bounds[k, 1], 0.0)
        upper = outlay_bounds[k, 2]
        depth_cap = if k == nlayers && total_height > upper
            total_height
        else
            min(total_height, upper)
        end
        layer_thickness = max(depth_cap - lower, 0.0)
        layer_thickness <= 0.0 && continue

        rho = clamp(valid_or(300.0, Float64(density_profile[k])), 50.0, c.rho_i)
        temp_k = clamp(valid_or(-10.0, Float64(temperature_profile_c[k])) + c.T0, 200.0, c.T0)
        water_fraction = max(valid_or(0.0, Float64(liquid_water_profile[k])), 0.0)
        snow_mass = rho * layer_thickness
        liquid_mass = water_fraction * snow_mass
        append_restart_layer!(snow_mass, liquid_mass, rho, temp_k)
    end

    while length(layer_mass) > ntot
        merge_restart_bottom_pair!()
    end

    return (
        N=length(layer_mass),
        mass=layer_mass,
        mass_w=layer_mass_w,
        density=layer_density,
        temperature=layer_temperature,
    )
end

function populate_domain_column_from_mar!(
    domain::SM.SnowpackDomain,
    idx::Int,
    total_height::Float64,
    density_profile::AbstractVector{<:Real},
    temperature_profile_c::AbstractVector{<:Real},
    liquid_water_profile::AbstractVector{<:Real},
    outlay_bounds::AbstractMatrix{<:Real},
)
    layers = extract_mar_layers(
        total_height,
        density_profile,
        temperature_profile_c,
        liquid_water_profile,
        outlay_bounds;
        ntot=domain.Ntot,
        c=domain.c,
        mass_split=domain.mass_split,
    )

    domain.N[idx] = layers.N
    @views domain.mass[:, idx] .= 0.0
    @views domain.mass_w[:, idx] .= 0.0
    @views domain.density[:, idx] .= 0.0
    @views domain.temperature[:, idx] .= domain.c.T0
    domain.mass_base[idx] = 0.0
    domain.smb_ice[idx] = 0.0
    domain.runoff[idx] = 0.0
    domain.snow_cover[idx] = 0.0
    domain.albedo_dynamic[idx] = layers.N > 0 ? domain.c.alpha_dry : domain.c.alpha_ice
    if layers.N > 0
        @inbounds for k in 1:layers.N
            domain.mass[k, idx] = layers.mass[k]
            domain.mass_w[k, idx] = layers.mass_w[k]
            domain.density[k, idx] = layers.density[k]
            domain.temperature[k, idx] = layers.temperature[k]
        end
        domain.Tsrf[idx] = domain.temperature[1, idx]
    else
        domain.Tsrf[idx] = domain.c.T0
    end
    SM.compute_auxiliary!(domain, idx)
    return nothing
end

function read_full_timeseries_3d(nc_path::AbstractString, varname::AbstractString, shapes::Dict{String, Vector{Int}})
    data = read_hdf5_full(nc_path, varname, shapes)
    if ndims(data) == 4
        size(data, 2) == 1 || error("Variable '$varname' has an unexpected non-singleton vertical dimension.")
        return dropdims(data; dims=(2,))
    elseif ndims(data) == 3
        return data
    end
    error("Variable '$varname' does not have a supported timeseries layout.")
end

function read_first_available_timeseries_3d(
    nc_path::AbstractString,
    candidate_names::Vector{String},
    shapes::Dict{String, Vector{Int}},
)
    for name in candidate_names
        if haskey(shapes, name)
            return (name=name, data=read_full_timeseries_3d(nc_path, name, shapes))
        end
    end
    return nothing
end

function load_gris_mar_problem(
    nc_path::AbstractString;
    mask_threshold::Float64=50.0,
    turbulent_flux_sign::Float64=1.0,
    ntot::Integer=20,
    physics::SM.SnowpackPhysicalConstants{Float64}=Chion.physics(),
)
    definition = Chion.mar_case(
        nc_path;
        mask_threshold=mask_threshold,
        turbulent_flux_sign=turbulent_flux_sign,
        physics=physics,
        ntot=ntot,
    )
    return (
        domain=definition.domain,
        forcing=definition.forcing,
        layout=definition.layout,
        wind_forcing_message=isempty(definition.notes) ? "" : definition.notes[1],
        nvalid=hasproperty(definition.metadata, :ncol) ? definition.metadata.ncol : SM.column_count(definition.domain),
        mask_threshold=mask_threshold,
        nc_path=String(nc_path),
    )
end

end
