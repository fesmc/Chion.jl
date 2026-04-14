#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Printf
using Base.Threads
using Chion

const SM = Chion.SnowpackModel
const DEFAULT_GRIS_FORCING_FILE_PATH = let path = ""
    for candidate in (
        "/p/projects/ou/labs/ai/Nils/MARv3.14.3-10km-daily-ERA5-2025.nc",
        "/Users/niboch001/Downloads/MARv3.14.3-10km-daily-ERA5-2025.nc",
    )
        if isfile(candidate)
            path = candidate
            break
        end
    end
    path
end
const DEFAULT_GRIS_FORCING_FILE_OUTPUT_DIR = joinpath(@__DIR__, "..", "plots", "gris_forcing_file_case")

function _script_definition_from_forcing_file(
    forcing_file::AbstractString;
    mask_threshold::Union{Nothing, Float64}=50.0,
    ntot::Integer=20,
    physics::SM.SnowpackPhysicalConstants{Float64}=Chion.physics(),
)
    if isnothing(mask_threshold)
        return Chion.prescribed_case(
            forcing_file=forcing_file,
            physics=physics,
            ntot=ntot,
            run=Chion.RunConfig(
                write_outputs=false,
                write_netcdf=false,
                cycles=1,
            ),
        ).definition
    end

    source_path = abspath(String(forcing_file))
    isfile(source_path) || error("Forcing file was not found: $(source_path)")

    shapes = Chion._forcing_file_read_dataset_shapes(source_path)
    required_variables = ("x", "y", "MSK", "OUTLAY_bnds", "TT", "SF", "RF", "SWD", "LWD", "SHF", "LHF", "ZN3", "RO1", "TI1", "WA1", "YYYY", "MM", "DD", "HH")
    missing = String[var for var in required_variables if !haskey(shapes, var)]
    isempty(missing) || error(
        "Forcing file $(abspath(source_path)) is missing required variables: $(join(missing, ", "))."
    )

    time_values = Chion._forcing_file_read_times(source_path, shapes)
    dt_days = [Chion._forcing_file_infer_dt_days(time_values, t) for t in eachindex(time_values)]

    x = vec(Chion._forcing_file_read_hdf5_full(source_path, "x", shapes))
    y = vec(Chion._forcing_file_read_hdf5_full(source_path, "y", shapes))
    mask = Chion._forcing_file_read_hdf5_full(source_path, "MSK", shapes)
    outlay_bounds = Chion._forcing_file_read_hdf5_full(source_path, "OUTLAY_bnds", shapes)

    tt_full = Chion._forcing_file_read_full_timeseries_3d(source_path, "TT", shapes)
    sf_full = Chion._forcing_file_read_full_timeseries_3d(source_path, "SF", shapes)
    rf_full = Chion._forcing_file_read_full_timeseries_3d(source_path, "RF", shapes)
    swd_full = Chion._forcing_file_read_full_timeseries_3d(source_path, "SWD", shapes)
    lwd_full = Chion._forcing_file_read_full_timeseries_3d(source_path, "LWD", shapes)
    shf_full = Chion._forcing_file_read_full_timeseries_3d(source_path, "SHF", shapes)
    lhf_full = Chion._forcing_file_read_full_timeseries_3d(source_path, "LHF", shapes)

    u_wind_info = Chion._forcing_file_read_first_available_timeseries_3d(source_path, ["UU", "U10"], shapes)
    v_wind_info = Chion._forcing_file_read_first_available_timeseries_3d(source_path, ["VV", "V10"], shapes)
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

    zn3_init = Chion._forcing_file_read_timeslice_2d(source_path, "ZN3", 1, shapes)
    ro1_init = Chion._forcing_file_read_timeslice_3d(source_path, "RO1", 1, shapes)
    ti1_init = Chion._forcing_file_read_timeslice_3d(source_path, "TI1", 1, shapes)
    wa1_init = Chion._forcing_file_read_timeslice_3d(source_path, "WA1", 1, shapes)

    ny, nx = size(mask)
    valid_mask = falses(ny, nx)
    @inbounds for j in 1:ny, i in 1:nx
        valid_mask[j, i] =
            isfinite(mask[j, i]) &&
            mask[j, i] >= Float64(mask_threshold) &&
            isfinite(tt_full[1, j, i])
    end
    valid_indices = findall(valid_mask)
    nvalid = length(valid_indices)
    nvalid > 0 || error(
        "No valid forcing-file grid cells remain after applying `mask_threshold=$(Float64(mask_threshold))`."
    )
    ntime = length(time_values)

    js = Vector{Int}(undef, nvalid)
    is = Vector{Int}(undef, nvalid)
    domain = SM.SnowpackDomain(ncol=nvalid, Ntot=Int(ntot), c=physics)
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
        Chion._forcing_file_populate_domain_column_from_restart!(
            domain,
            idx,
            Float64(zn3_init[j, i]),
            @view(ro1_init[:, j, i]),
            @view(ti1_init[:, j, i]),
            @view(wa1_init[:, j, i]),
            outlay_bounds,
        )
        for t in 1:ntime
            tair_k[idx, t] = Chion._forcing_file_valid_or(-15.0, Float64(tt_full[t, j, i])) + domain.c.T0
            snow_rate[idx, t] = Chion._forcing_file_mmwe_day_to_kgm2s(Float64(sf_full[t, j, i]))
            rain_rate[idx, t] = Chion._forcing_file_mmwe_day_to_kgm2s(Float64(rf_full[t, j, i]))
            s_boa[idx, t] = Chion._forcing_file_valid_or(0.0, Float64(swd_full[t, j, i]))
            q_lw_ij = Float64(lwd_full[t, j, i])
            q_sh_ij = Float64(shf_full[t, j, i])
            q_lh_ij = Float64(lhf_full[t, j, i])
            has_q_lw[idx, t] = isfinite(q_lw_ij)
            has_q_sh[idx, t] = isfinite(q_sh_ij)
            has_q_lh[idx, t] = isfinite(q_lh_ij)
            q_lw[idx, t] = has_q_lw[idx, t] ? q_lw_ij : 0.0
            q_sh[idx, t] = has_q_sh[idx, t] ? q_sh_ij : 0.0
            q_lh[idx, t] = has_q_lh[idx, t] ? q_lh_ij : 0.0
            wind_speed[idx, t] = isnothing(wind_full) ? 5.0 : Chion._forcing_file_valid_or(5.0, Float64(wind_full[t, j, i]))
        end
    end

    forcing = Chion.ForcingData(
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
    layout = Chion.GridLayout(x, y, js, is, mask)
    metadata = (
        format=:prescribed,
        source=:file,
        path=source_path,
        ncol=nvalid,
        ntime=ntime,
        ntot=Int(ntot),
        masked_by_script=true,
        mask_threshold=Float64(mask_threshold),
    )
    return Chion.CaseDefinition(
        domain,
        forcing;
        layout=layout,
        input_label=source_path,
        notes=[wind_note],
        metadata=metadata,
    )
end

function load_gris_forcing_file_problem(
    forcing_file::AbstractString;
    mask_threshold::Union{Nothing, Float64}=50.0,
    ntot::Integer=15,
    physics::SM.SnowpackPhysicalConstants{Float64}=Chion.physics(),
)
    definition = _script_definition_from_forcing_file(
        forcing_file;
        mask_threshold=mask_threshold,
        ntot=ntot,
        physics=physics,
    )
    return (
        domain=definition.domain,
        forcing=definition.forcing,
        layout=definition.layout,
        wind_forcing_message=isempty(definition.notes) ? "" : definition.notes[1],
        nvalid=hasproperty(definition.metadata, :ncol) ? definition.metadata.ncol : SM.column_count(definition.domain),
        forcing_file=String(forcing_file),
        mask_threshold=mask_threshold,
    )
end

function run_gris_forcing_file_case(
    forcing_file::AbstractString;
    io::IO=stdout,
    mask_threshold::Union{Nothing, Float64}=50.0,
    ntot::Integer=20,
    run::Chion.RunConfig=Chion.RunConfig(),
    physics::SM.SnowpackPhysicalConstants{Float64}=Chion.physics(),
)
    definition = _script_definition_from_forcing_file(
        forcing_file;
        mask_threshold=mask_threshold,
        ntot=ntot,
        physics=physics,
    )
    if !isnothing(mask_threshold)
        println(io, "Applied MAR mask threshold $(Float64(mask_threshold)) during script-side loading.")
    end

    d0 = definition.domain
    """empty_domain = SM.SnowpackDomain(
        c=d0.c,
        ncol=d0.ncol,
        Ntot=15,
        mass_max=d0.mass_max,
        mass_split=d0.mass_split,
        mass_min=d0.mass_min,
        rho_max=d0.rho_max,
        density_init=d0.c.rho_s,
        temperature_init=d0.c.T0,
    )
    empty_domain.Tsrf .= d0.c.T0
    empty_domain.albedo_dynamic .= d0.c.alpha_dry

    definition = Chion.CaseDefinition(
        empty_domain,
        definition.forcing;
        layout=definition.layout,
        input_label=definition.input_label,
        notes=definition.notes,
        metadata=definition.metadata,
    )"""
    case = Chion.SnowpackCase(
        definition;
        run=Chion.RunConfig(
            name=isempty(run.name) ? "GrIS forcing file case" : run.name,
            input_label=isempty(run.input_label) ? abspath(forcing_file) : run.input_label,
            output_dir=run.output_dir,
            netcdf_path=run.netcdf_path,
            write_outputs=run.write_outputs,
            write_netcdf=run.write_netcdf,
            netcdf_variables=run.netcdf_variables,
            cycles=run.cycles,
            history_stride=run.history_stride,
            backend=run.backend,
        ),
    )
    return Chion.run_case(case; io=io)
end
