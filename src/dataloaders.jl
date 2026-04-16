using HDF5

const FILL_THRESHOLD = -9.0e18

struct LoadedProblem{D,F,L,M}
    domain::D
    forcing::F
    layout::L
    notes::Vector{String}
    metadata::M
end

export LoadedProblem
export read_dataset_shapes, read_hdf5_subset, read_hdf5_full
export read_timeslice_2d, read_timeslice_3d
export valid_or, mmwe_day_to_kgm2s
export read_forcing_times, choose_time_index, infer_dt_days
export extract_forcing_file_layers, populate_domain_column_from_forcing_file!
export read_full_timeseries_3d, read_first_available_timeseries_3d
export load_gris_forcing_file_problem

function read_dataset_shapes(nc_path::AbstractString)
    shapes = Dict{String, Vector{Int}}()
    h5open(nc_path, "r") do file
        for name in keys(file)
            obj = file[name]
            if obj isa HDF5.Dataset
                # Preserve the logical dimension order used by Chion:
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

function read_forcing_times(nc_path::AbstractString, shapes::Dict{String, Vector{Int}})
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

function extract_forcing_file_layers(
    total_height::Float64,
    density_profile::AbstractVector{<:Real},
    temperature_profile_c::AbstractVector{<:Real},
    liquid_water_profile::AbstractVector{<:Real},
    outlay_bounds::AbstractMatrix{<:Real};
    ntot::Int=80,
    c::SnowpackPhysicalConstants=SnowpackPhysicalConstants(),
    mass_split::Float64=DEFAULT_MASS_SPLIT,
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
        # Preserve the native restart layering on restart.
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

function populate_domain_column_from_forcing_file!(
    domain::SnowpackDomain,
    idx::Int,
    total_height::Float64,
    density_profile::AbstractVector{<:Real},
    temperature_profile_c::AbstractVector{<:Real},
    liquid_water_profile::AbstractVector{<:Real},
    outlay_bounds::AbstractMatrix{<:Real},
)
    layers = extract_forcing_file_layers(
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
    compute_auxiliary!(domain, idx)
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

function load_gris_forcing_file_problem(
    nc_path::AbstractString;
    mask_threshold::Union{Nothing, Float64}=50.0,
    ntot::Integer=20,
    physics::SnowpackPhysicalConstants{Float64}=physics(),
)
    source_path = abspath(String(nc_path))
    isfile(source_path) || error("Forcing file was not found: $(source_path)")

    shapes = read_dataset_shapes(source_path)
    required_variables = ("x", "y", "MSK", "OUTLAY_bnds", "TT", "SF", "RF", "SWD", "LWD", "SHF", "LHF", "ZN3", "RO1", "TI1", "WA1", "YYYY", "MM", "DD", "HH")
    missing = String[var for var in required_variables if !haskey(shapes, var)]
    isempty(missing) || error(
        "Forcing file $(abspath(source_path)) is missing required variables: $(join(missing, ", "))."
    )

    time_values = read_forcing_times(source_path, shapes)[2]
    dt_days = [infer_dt_days(time_values, t) for t in eachindex(time_values)]

    x = vec(read_hdf5_full(source_path, "x", shapes))
    y = vec(read_hdf5_full(source_path, "y", shapes))
    mask = read_hdf5_full(source_path, "MSK", shapes)
    outlay_bounds = read_hdf5_full(source_path, "OUTLAY_bnds", shapes)

    tt_full = read_full_timeseries_3d(source_path, "TT", shapes)
    sf_full = read_full_timeseries_3d(source_path, "SF", shapes)
    rf_full = read_full_timeseries_3d(source_path, "RF", shapes)
    swd_full = read_full_timeseries_3d(source_path, "SWD", shapes)
    lwd_full = read_full_timeseries_3d(source_path, "LWD", shapes)
    shf_full = read_full_timeseries_3d(source_path, "SHF", shapes)
    lhf_full = read_full_timeseries_3d(source_path, "LHF", shapes)

    u_wind_info = read_first_available_timeseries_3d(source_path, ["UU", "U10"], shapes)
    v_wind_info = read_first_available_timeseries_3d(source_path, ["VV", "V10"], shapes)
    wind_full = if !isnothing(u_wind_info) && !isnothing(v_wind_info)
        hypot.(u_wind_info.data, v_wind_info.data)
    else
        nothing
    end
    wind_note = if isnothing(wind_full)
        "Wind forcing: file wind components not found; using default 5.0 m s^-1."
    else
        @sprintf(
            "Wind forcing: |V| from components %s and %s.",
            u_wind_info.name,
            v_wind_info.name,
        )
    end

    zn3_init = read_timeslice_2d(source_path, "ZN3", 1, shapes)
    ro1_init = read_timeslice_3d(source_path, "RO1", 1, shapes)
    ti1_init = read_timeslice_3d(source_path, "TI1", 1, shapes)
    wa1_init = read_timeslice_3d(source_path, "WA1", 1, shapes)

    ny, nx = size(mask)
    threshold = isnothing(mask_threshold) ? -Inf : Float64(mask_threshold)
    valid_mask = falses(ny, nx)
    @inbounds for j in 1:ny, i in 1:nx
        valid_mask[j, i] =
            isfinite(mask[j, i]) &&
            mask[j, i] >= threshold &&
            isfinite(tt_full[1, j, i])
    end
    valid_indices = findall(valid_mask)
    nvalid = length(valid_indices)
    nvalid > 0 || error(
        "No valid forcing-file grid cells remain after applying `mask_threshold=$(threshold)`."
    )
    ntime = length(time_values)

    js = Vector{Int}(undef, nvalid)
    is = Vector{Int}(undef, nvalid)
    domain = SnowpackDomain(ncol=nvalid, Ntot=Int(ntot), c=physics)
    tair_k = Matrix{Float64}(undef, nvalid, ntime)
    snow_rate = Matrix{Float64}(undef, nvalid, ntime)
    rain_rate = Matrix{Float64}(undef, nvalid, ntime)
    s_boa = Matrix{Float64}(undef, nvalid, ntime)
    q_lw = Matrix{Float64}(undef, nvalid, ntime)
    has_q_lw = fill(false, nvalid, ntime)
    q_sh = Matrix{Float64}(undef, nvalid, ntime)
    has_q_sh = fill(false, nvalid, ntime)
    q_lh = Matrix{Float64}(undef, nvalid, ntime)
    has_q_lh = fill(false, nvalid, ntime)
    wind_speed = Matrix{Float64}(undef, nvalid, ntime)

    @threads :static for idx in eachindex(valid_indices)
        j, i = Tuple(valid_indices[idx])
        js[idx] = j
        is[idx] = i
        populate_domain_column_from_forcing_file!(
            domain,
            idx,
            Float64(zn3_init[j, i]),
            @view(ro1_init[:, j, i]),
            @view(ti1_init[:, j, i]),
            @view(wa1_init[:, j, i]),
            outlay_bounds,
        )
        for t in 1:ntime
            tair_k[idx, t] = valid_or(-15.0, Float64(tt_full[t, j, i])) + domain.c.T0
            snow_rate[idx, t] = mmwe_day_to_kgm2s(Float64(sf_full[t, j, i]))
            rain_rate[idx, t] = mmwe_day_to_kgm2s(Float64(rf_full[t, j, i]))
            s_boa[idx, t] = valid_or(0.0, Float64(swd_full[t, j, i]))
            q_lw_ij = Float64(lwd_full[t, j, i])
            q_sh_ij = Float64(shf_full[t, j, i])
            q_lh_ij = Float64(lhf_full[t, j, i])
            has_q_lw[idx, t] = isfinite(q_lw_ij)
            has_q_sh[idx, t] = isfinite(q_sh_ij)
            has_q_lh[idx, t] = isfinite(q_lh_ij)
            q_lw[idx, t] = has_q_lw[idx, t] ? q_lw_ij : 0.0
            q_sh[idx, t] = has_q_sh[idx, t] ? q_sh_ij : 0.0
            q_lh[idx, t] = has_q_lh[idx, t] ? q_lh_ij : 0.0
            wind_speed[idx, t] = isnothing(wind_full) ? 5.0 : valid_or(5.0, Float64(wind_full[t, j, i]))
        end
    end

    forcing = ForcingData(
        time_values=time_values,
        dt_days=dt_days,
        air_temperature=tair_k,
        snowfall_rate=snow_rate,
        rainfall_rate=rain_rate,
        shortwave_down=s_boa,
        wind_speed=wind_speed,
        q_lw_down=q_lw,
        has_q_lw_down=has_q_lw,
        q_sh=q_sh,
        has_q_sh=has_q_sh,
        q_lh=q_lh,
        has_q_lh=has_q_lh,
    )
    layout = GridLayout(x, y, js, is, mask)
    metadata = (
        format=:prescribed,
        source=:file,
        path=source_path,
        ncol=nvalid,
        ntime=ntime,
        ntot=Int(ntot),
        masked_by_script=true,
        mask_threshold=isnothing(mask_threshold) ? nothing : Float64(mask_threshold),
    )
    notes = [wind_note]
    return LoadedProblem(
        domain,
        forcing,
        layout,
        notes,
        metadata,
    )
end
