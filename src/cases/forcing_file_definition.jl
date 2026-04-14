const _FORCING_FILE_FILL_THRESHOLD = -9.0e18

function _forcing_file_read_dataset_shapes(path::AbstractString)
    shapes = Dict{String, Vector{Int}}()
    h5open(path, "r") do file
        for name in keys(file)
            obj = file[name]
            if obj isa HDF5.Dataset
                shapes[String(name)] = reverse(collect(size(obj)))
            end
        end
    end
    return shapes
end

function _forcing_file_clean_fill!(A)
    @inbounds for i in eachindex(A)
        if A[i] <= _FORCING_FILE_FILL_THRESHOLD
            A[i] = NaN
        end
    end
    return A
end

function _forcing_file_read_hdf5_subset(
    path::AbstractString,
    varname::AbstractString,
    full_shape::Vector{Int};
    start::Union{Nothing, Vector{Int}}=nothing,
    count::Union{Nothing, Vector{Int}}=nothing,
)
    start_vec = isnothing(start) ? zeros(Int, length(full_shape)) : start
    count_vec = isnothing(count) ? copy(full_shape) : count
    length(start_vec) == length(full_shape) || error("start rank mismatch for variable '$varname'.")
    length(count_vec) == length(full_shape) || error("count rank mismatch for variable '$varname'.")
    data = h5open(path, "r") do file
        haskey(file, varname) || error("Variable '$varname' was not found in $(abspath(path)).")
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
    _forcing_file_clean_fill!(data)
    return data
end

function _forcing_file_read_hdf5_full(path::AbstractString, varname::AbstractString, shapes::Dict{String, Vector{Int}})
    haskey(shapes, varname) || error("Variable '$varname' was not found in $(abspath(path)).")
    return _forcing_file_read_hdf5_subset(path, varname, shapes[varname])
end

function _forcing_file_read_timeslice_2d(
    path::AbstractString,
    varname::AbstractString,
    time_index::Int,
    shapes::Dict{String, Vector{Int}},
)
    haskey(shapes, varname) || error("Variable '$varname' was not found in $(abspath(path)).")
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
    data = _forcing_file_read_hdf5_subset(path, varname, shape; start=start, count=count)
    return dropdims(data; dims=Tuple(findall(==(1), size(data))))
end

function _forcing_file_read_timeslice_3d(
    path::AbstractString,
    varname::AbstractString,
    time_index::Int,
    shapes::Dict{String, Vector{Int}},
)
    haskey(shapes, varname) || error("Variable '$varname' was not found in $(abspath(path)).")
    shape = shapes[varname]
    length(shape) == 4 || error("Variable '$varname' does not have the expected 4D layout.")
    start = [time_index - 1, 0, 0, 0]
    count = [1, shape[2], shape[3], shape[4]]
    data = _forcing_file_read_hdf5_subset(path, varname, shape; start=start, count=count)
    return dropdims(data; dims=(1,))
end

@inline _forcing_file_valid_or(default::Float64, x::Float64) = isfinite(x) ? x : default
@inline _forcing_file_mmwe_day_to_kgm2s(x::Float64) = isfinite(x) ? max(x, 0.0) / 86_400.0 : 0.0

function _forcing_file_read_times(path::AbstractString, shapes::Dict{String, Vector{Int}})
    yyyy = round.(Int, vec(_forcing_file_read_hdf5_full(path, "YYYY", shapes)))
    mm = round.(Int, vec(_forcing_file_read_hdf5_full(path, "MM", shapes)))
    dd = round.(Int, vec(_forcing_file_read_hdf5_full(path, "DD", shapes)))
    hh = round.(Int, vec(_forcing_file_read_hdf5_full(path, "HH", shapes)))
    ntime = length(yyyy)
    times = Vector{DateTime}(undef, ntime)
    for i in 1:ntime
        times[i] = DateTime(yyyy[i], mm[i], dd[i], hh[i])
    end
    return times
end

function _forcing_file_infer_dt_days(time_values::Vector{DateTime}, time_index::Int)
    if length(time_values) == 1
        return 1.0
    elseif time_index < length(time_values)
        return Dates.value(time_values[time_index + 1] - time_values[time_index]) / (1000 * 60 * 60 * 24)
    else
        return Dates.value(time_values[time_index] - time_values[time_index - 1]) / (1000 * 60 * 60 * 24)
    end
end

function _forcing_file_extract_layers(
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

        rho = clamp(_forcing_file_valid_or(300.0, Float64(density_profile[k])), 50.0, c.rho_i)
        temp_k = clamp(_forcing_file_valid_or(-10.0, Float64(temperature_profile_c[k])) + c.T0, 200.0, c.T0)
        water_fraction = max(_forcing_file_valid_or(0.0, Float64(liquid_water_profile[k])), 0.0)
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

function _forcing_file_populate_domain_column_from_restart!(
    domain::SnowpackDomain,
    idx::Int,
    total_height::Float64,
    density_profile::AbstractVector{<:Real},
    temperature_profile_c::AbstractVector{<:Real},
    liquid_water_profile::AbstractVector{<:Real},
    outlay_bounds::AbstractMatrix{<:Real},
)
    extracted = _forcing_file_extract_layers(
        total_height,
        density_profile,
        temperature_profile_c,
        liquid_water_profile,
        outlay_bounds;
        ntot=domain.Ntot,
        c=domain.c,
        mass_split=domain.mass_split,
    )
    domain.N[idx] = extracted.N
    if extracted.N > 0
        domain.mass[1:extracted.N, idx] .= extracted.mass
        domain.mass_w[1:extracted.N, idx] .= extracted.mass_w
        domain.density[1:extracted.N, idx] .= extracted.density
        domain.temperature[1:extracted.N, idx] .= extracted.temperature
        if extracted.N < domain.Ntot
            domain.mass[(extracted.N + 1):end, idx] .= 0.0
            domain.mass_w[(extracted.N + 1):end, idx] .= 0.0
            domain.density[(extracted.N + 1):end, idx] .= DEFAULT_DENSITY_INIT
            domain.temperature[(extracted.N + 1):end, idx] .= domain.c.T0 - 10.0
        end
        domain.Tsrf[idx] = extracted.temperature[1]
    else
        domain.mass[:, idx] .= 0.0
        domain.mass_w[:, idx] .= 0.0
        domain.density[:, idx] .= DEFAULT_DENSITY_INIT
        domain.temperature[:, idx] .= domain.c.T0 - 10.0
        domain.Tsrf[idx] = domain.c.T0 - 10.0
    end
    domain.mass_base[idx] = 0.0
    domain.smb_ice[idx] = 0.0
    domain.runoff[idx] = 0.0
    domain.snow_cover[idx] = 0.0
    domain.albedo_dynamic[idx] = domain.c.alpha_dry
    return domain
end

function _forcing_file_read_full_timeseries_3d(
    path::AbstractString,
    varname::AbstractString,
    shapes::Dict{String, Vector{Int}},
)
    haskey(shapes, varname) || error("Variable '$varname' was not found in $(abspath(path)).")
    shape = shapes[varname]
    if length(shape) == 4
        data = _forcing_file_read_hdf5_full(path, varname, shapes)
        size(data, 2) == 1 || error("Variable '$varname' has an unexpected non-singleton vertical dimension.")
        return dropdims(data; dims=(2,))
    elseif length(shape) == 3
        return _forcing_file_read_hdf5_full(path, varname, shapes)
    end
    error("Variable '$varname' does not have a supported timeseries layout.")
end

function _forcing_file_read_first_available_timeseries_3d(
    path::AbstractString,
    candidate_names::Vector{String},
    shapes::Dict{String, Vector{Int}},
)
    for name in candidate_names
        if haskey(shapes, name)
            return (name=name, data=_forcing_file_read_full_timeseries_3d(path, name, shapes))
        end
    end
    return nothing
end

"""
    _prescribed_definition_from_forcing_file(forcing_file; physics=physics(), ntot=20)

Load a reusable [`CaseDefinition`](@ref) from a prepared external forcing file.
"""
function _prescribed_definition_from_forcing_file(
    forcing_file::AbstractString;
    physics::SnowpackPhysicalConstants{Float64}=physics(),
    ntot::Integer=20,
)
    Int(ntot) > 0 || error("`ntot` must be positive.")
    isempty(strip(forcing_file)) && error("`forcing_file` must point to a prepared forcing file.")
    source_path = abspath(String(forcing_file))
    isfile(source_path) || error("Forcing file was not found: $(source_path)")

    shapes = _forcing_file_read_dataset_shapes(source_path)
    required_variables = ("x", "y", "MSK", "OUTLAY_bnds", "TT", "SF", "RF", "SWD", "LWD", "SHF", "LHF", "ZN3", "RO1", "TI1", "WA1", "YYYY", "MM", "DD", "HH")
    missing = String[var for var in required_variables if !haskey(shapes, var)]
    isempty(missing) || error(
        "Forcing file $(abspath(source_path)) is missing required variables: $(join(missing, ", "))."
    )

    time_values = _forcing_file_read_times(source_path, shapes)
    dt_days = [_forcing_file_infer_dt_days(time_values, t) for t in eachindex(time_values)]

    x = vec(_forcing_file_read_hdf5_full(source_path, "x", shapes))
    y = vec(_forcing_file_read_hdf5_full(source_path, "y", shapes))
    mask = _forcing_file_read_hdf5_full(source_path, "MSK", shapes)
    outlay_bounds = _forcing_file_read_hdf5_full(source_path, "OUTLAY_bnds", shapes)

    tt_full = _forcing_file_read_full_timeseries_3d(source_path, "TT", shapes)
    sf_full = _forcing_file_read_full_timeseries_3d(source_path, "SF", shapes)
    rf_full = _forcing_file_read_full_timeseries_3d(source_path, "RF", shapes)
    swd_full = _forcing_file_read_full_timeseries_3d(source_path, "SWD", shapes)
    lwd_full = _forcing_file_read_full_timeseries_3d(source_path, "LWD", shapes)
    shf_full = _forcing_file_read_full_timeseries_3d(source_path, "SHF", shapes)
    lhf_full = _forcing_file_read_full_timeseries_3d(source_path, "LHF", shapes)

    u_wind_info = _forcing_file_read_first_available_timeseries_3d(source_path, ["UU", "U10"], shapes)
    v_wind_info = _forcing_file_read_first_available_timeseries_3d(source_path, ["VV", "V10"], shapes)
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

    zn3_init = _forcing_file_read_timeslice_2d(source_path, "ZN3", 1, shapes)
    ro1_init = _forcing_file_read_timeslice_3d(source_path, "RO1", 1, shapes)
    ti1_init = _forcing_file_read_timeslice_3d(source_path, "TI1", 1, shapes)
    wa1_init = _forcing_file_read_timeslice_3d(source_path, "WA1", 1, shapes)

    ny, nx = size(mask)
    valid_mask = falses(ny, nx)
    @inbounds for j in 1:ny, i in 1:nx
        valid_mask[j, i] =
            all(isfinite, @view(tt_full[:, j, i])) &&
            all(isfinite, @view(sf_full[:, j, i])) &&
            all(isfinite, @view(rf_full[:, j, i])) &&
            all(isfinite, @view(swd_full[:, j, i]))
    end
    valid_indices = findall(valid_mask)
    nvalid = length(valid_indices)
    nvalid > 0 || error("No valid forcing-file grid cells remain after filtering for finite required prescribed forcing.")
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
        _forcing_file_populate_domain_column_from_restart!(
            domain,
            idx,
            Float64(zn3_init[j, i]),
            @view(ro1_init[:, j, i]),
            @view(ti1_init[:, j, i]),
            @view(wa1_init[:, j, i]),
            outlay_bounds,
        )
        for t in 1:ntime
            tair_k[idx, t] = _forcing_file_valid_or(-15.0, Float64(tt_full[t, j, i])) + domain.c.T0
            snow_rate[idx, t] = _forcing_file_mmwe_day_to_kgm2s(Float64(sf_full[t, j, i]))
            rain_rate[idx, t] = _forcing_file_mmwe_day_to_kgm2s(Float64(rf_full[t, j, i]))
            s_boa[idx, t] = _forcing_file_valid_or(0.0, Float64(swd_full[t, j, i]))
            q_lw_ij = Float64(lwd_full[t, j, i])
            q_sh_ij = Float64(shf_full[t, j, i])
            q_lh_ij = Float64(lhf_full[t, j, i])
            has_q_lw[idx, t] = isfinite(q_lw_ij)
            has_q_sh[idx, t] = isfinite(q_sh_ij)
            has_q_lh[idx, t] = isfinite(q_lh_ij)
            q_lw[idx, t] = has_q_lw[idx, t] ? q_lw_ij : 0.0
            q_sh[idx, t] = has_q_sh[idx, t] ? q_sh_ij : 0.0
            q_lh[idx, t] = has_q_lh[idx, t] ? q_lh_ij : 0.0
            wind_speed[idx, t] = isnothing(wind_full) ? 5.0 : _forcing_file_valid_or(5.0, Float64(wind_full[t, j, i]))
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
    )
    return CaseDefinition(
        domain,
        forcing;
        layout=layout,
        input_label=source_path,
        notes=[wind_note],
        metadata=metadata,
    )
end
