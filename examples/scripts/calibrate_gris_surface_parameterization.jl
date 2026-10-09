#!/usr/bin/env julia
# Calibrate one parameterized surface process against MAR with real daily forcing,
# keeping every other surface flux prescribed so the process is isolated:
#   CHION_CALIBRATION=albedo    all fluxes prescribed (MAR SW net), dynamic albedo
#                               diagnosed alongside; scored on the absorbed
#                               shortwave implied by Chion's albedo, SWD·(1−α),
#                               against SWD·(1−α_MAR).
#   CHION_CALIBRATION=states    only build and cache the shared starting states
#   CHION_CALIBRATION=sensible  MAR SW net, LWD and LHF prescribed, SEMIX sensible
#                               heat parameterized (coupled to Ts); scored on the
#                               daily sensible heat flux against MAR SHF.
#   CHION_CALIBRATION=diagnose  no calibration: free runs (prescribed, and parameterized
#                               SEB with MAR albedo) from the same starting state, with
#                               monthly SW, LW, SH, LH and melt by elevation band vs MAR.
# Snowpack starting states come from one prescribed spin-up + transient, cached
# and shared by both calibrations. Calibration years are scored first; the best
# sets are then re-scored on independent validation years.
include(joinpath(@__DIR__, "run_gris_fully_prescribed_transient.jl"))
using Serialization, Printf, Dates

parse_years(s) = (r = parse.(Int, split(s, ':')); r[1]:r[end])
const CAL = (
    mode=Symbol(get(ENV, "CHION_CALIBRATION", "albedo")),
    out=get(ENV, "CHION_OUTPUT_DIR", joinpath(@__DIR__, "..", "..", "output", "gris_surface_calibration")),
    state_cache=get(ENV, "CHION_CALIBRATION_STATE_CACHE",
        joinpath(@__DIR__, "..", "..", "output", "gris_surface_calibration", "prescribed_states.jls")),
    forcing_dir=get(ENV, "CHION_FORCING_DIR", "/p/projects/ou/labs/ai/Nils/MAR_daily/mar_daily"),
    cal_years=parse_years(get(ENV, "CHION_CALIBRATION_YEARS", "1980:1989")),
    val_years=parse_years(get(ENV, "CHION_VALIDATION_YEARS", "2010:2019")),
    n_validate=parse(Int, get(ENV, "CHION_VALIDATE_BEST", "5")),
    months=4:9,
)
const SECONDS_PER_DAY = 86_400.0
daily_file(year) = joinpath(CAL.forcing_dir, "MARv3.14.3-10km-daily-ERA5-$(year).nc")

"""Same model as `build_fully_prescribed_model`, with process parameters overridden."""
calibration_model(grid; overrides...) = BESSIModel(grid; merge(
    (Ntot=CONFIG.snow_layers, albedo=:prescribed, diurnal_shortwave_max_substeps=CONFIG.diurnal_substeps),
    model_env_overrides(), (; overrides...))...)

"""Rebuild `forcing` with selected fields replaced (SnowpackForcing is immutable)."""
function rebuild(f; kw...)
    fields = (time_values=f.time_values, dt_days=f.dt_days, air_temperature=f.air_temperature,
        snowfall_rate=f.snowfall_rate, rainfall_rate=f.rainfall_rate, shortwave_down=f.shortwave_down,
        wind_speed=f.wind_speed, q_sw_net=f.q_sw_net, has_q_sw_net=f.has_q_sw_net, q_lw_down=f.q_lw_down,
        has_q_lw_down=f.has_q_lw_down, q_sh=f.q_sh, has_q_sh=f.has_q_sh, q_lh=f.q_lh, has_q_lh=f.has_q_lh,
        relative_humidity=f.relative_humidity, has_relative_humidity=f.has_relative_humidity,
        air_pressure=f.air_pressure, surface_height=f.surface_height, prescribed_albedo=f.prescribed_albedo,
        has_prescribed_albedo=f.has_prescribed_albedo, latitude_deg=f.latitude_deg)
    return SnowpackForcing(; merge(fields, (; kw...))...)
end

"""The step kernel uses the parameters stored in the state, so bind each trial's parameters."""
rebind(state::Chion.BESSIState, model) = Chion.BESSIState(model.parameters, state.ncol,
    (deepcopy(getfield(state, name)) for name in Chion._BESSI_ARRAY_FIELD_NAMES)...)
rebind(snapshot::NamedTuple, model) = Chion.BESSIState(model.parameters, snapshot.ncol, deepcopy.(snapshot.arrays)...)

"""Parameter-free copy of a state's arrays, so the cache survives changes to the model constants."""
snapshot(state) = (ncol=state.ncol, arrays=Tuple(deepcopy(getfield(state, name)) for name in Chion._BESSI_ARRAY_FIELD_NAMES))

# ---------- shared prescribed starting states ----------
function prescribed_states(grid)
    isfile(CAL.state_cache) && return deserialize(CAL.state_cache)
    model = calibration_model(grid)
    clim = load_mar(CONFIG.spinup_file)
    spin = Simulation(model; forcing=fully_prescribed_forcing(CONFIG.spinup_file, clim), years=CONFIG.spinup_years,
        backend=CONFIG.backend, write_netcdf=false, compute_year_metrics=false, name="calibration_spinup")
    run!(spin).status == :complete || error("Calibration spin-up failed")
    state = spin.now
    starts = (first(CAL.cal_years), first(CAL.val_years))
    states = Dict{Int,Any}()
    for year in 1940:(maximum(starts) - 1)
        year in starts && (states[year] = snapshot(state))
        loaded = load_mar(daily_file(year))
        sim = Simulation(model; forcing=fully_prescribed_forcing(daily_file(year), loaded), state, years=1,
            backend=CONFIG.backend, write_netcdf=false, compute_year_metrics=false, name="calibration_prescribed_$(year)")
        run!(sim).status == :complete || error("Prescribed transient failed in $year")
        state = sim.now
        @info "Prescribed starting-state transient" year
        GC.gc()
    end
    states[maximum(starts)] = snapshot(state)
    mkpath(dirname(CAL.state_cache))
    serialize(CAL.state_cache, states)
    return states
end

# ---------- trial forcings and diagnostics ----------
function trial_forcing(year, setting)
    path = daily_file(year)
    loaded = load_mar(path)
    base = fully_prescribed_forcing(path, loaded)
    if CAL.mode == :albedo
        # Albedo evolves with the dynamic scheme but SW stays MAR's net shortwave.
        return rebuild(base; has_prescribed_albedo=false), base
    end
    wind = setting.wind == :mar ? read_lowest_wind_level(path, loaded.grid, length(base.time_values)) : base.wind_speed
    return rebuild(base; has_q_sh=false, wind_speed=wind), base
end

"""Run one year day by day, accumulating monthly sums per column."""
function run_year!(model, forcing, reference, state, acc)
    sim = Simulation(model; forcing, state=rebind(state, model), years=1, backend=CONFIG.backend,
        write_netcdf=false, compute_year_metrics=false, name="calibration_trial")
    integ = init_integrator(sim; io=devnull)
    runtime = integ.model_runtime.data
    previous_sh = copy(Array(runtime.workspace.sensible_heat_energy))
    for k in eachindex(forcing.dt_days)
        step!(integ)
        m = month(forcing.time_values[k])
        if !(m in CAL.months)
            CAL.mode == :sensible && copyto!(previous_sh, Array(runtime.workspace.sensible_heat_energy))
            continue
        end
        if CAL.mode == :albedo
            swd = @view reference.shortwave_down[:, k]
            alpha = Array(runtime.state.albedo)
            acc.chion[m, :] .+= swd .* (1 .- alpha)
            acc.mar[m, :] .+= @view reference.q_sw_net[:, k]
            acc.albedo_chion[m, :] .+= swd .* alpha
            acc.albedo_mar[m, :] .+= swd .* clamp.(@view(reference.prescribed_albedo[:, k]), 0, 1)
            acc.weight[m, :] .+= swd
        else
            sh = Array(runtime.workspace.sensible_heat_energy)
            acc.chion[m, :] .+= (sh .- previous_sh) ./ (forcing.dt_days[k] * SECONDS_PER_DAY)
            copyto!(previous_sh, sh)
            acc.mar[m, :] .+= @view reference.q_sh[:, k]
        end
        acc.days[m, :] .+= 1
    end
    finalize!(integ)
    return sim.now
end

new_acc(ncol) = (chion=zeros(12, ncol), mar=zeros(12, ncol), albedo_chion=zeros(12, ncol),
    albedo_mar=zeros(12, ncol), weight=zeros(12, ncol), days=zeros(12, ncol))

"""Scores in W m⁻²: area-weighted mean bias and RMSE of monthly means over the season, total and by elevation band."""
function score(acc, area, height)
    months = collect(CAL.months)
    days = acc.days[months, :]
    diff = (acc.chion[months, :] .- acc.mar[months, :]) ./ max.(days, 1)
    w = repeat(area', length(months))
    band(sel) = (s = repeat(sel', length(months)); (sum(diff[s] .* w[s]) / sum(w[s]), sqrt(sum(diff[s] .^ 2 .* w[s]) / sum(w[s]))))
    all_bias, all_rmse = band(trues(length(area)))
    low_bias, low_rmse = band(height .< 1500)
    high_bias, high_rmse = band(height .>= 1500)
    albedo_bias = CAL.mode == :albedo ?
        sum((acc.albedo_chion[months, :] .- acc.albedo_mar[months, :]) .* w ./ max.(acc.weight[months, :], 1)) / sum(w) : NaN
    return (; bias=all_bias, rmse=all_rmse, low_bias, low_rmse, high_bias, high_rmse, albedo_bias)
end

function evaluate(setting, years, start_state, grid, area, height, forcings)
    model = calibration_model(grid; setting.model...)
    acc = new_acc(length(area))
    state = start_state
    for year in years
        forcing, reference = forcings[year]
        state = run_year!(model, forcing, reference, state, acc)
    end
    return score(acc, area, height)
end

function settings()
    if CAL.mode == :albedo
        grid_values = Iterators.product((0.81, 0.84), (0.60, 0.65, 0.70, 0.75), (0.40, 0.45, 0.50, 0.55))
        return [(label=@sprintf("dry=%.2f wet=%.2f ice=%.2f", d, w, i),
                 model=(albedo=:dynamic, alpha_dry=d, alpha_wet=w, alpha_ice=i), wind=:constant)
                for (d, w, i) in grid_values if w < d]
    end
    grid_values = Iterators.product((1.5, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0), (2.0, 5.0, 10.0, 20.0, 40.0), (:constant, :mar))
    return [(label=@sprintf("factor=%.2f stable=%.0f wind=%s", f, s, w),
             model=(semix_sensible_exchange_factor=f, semix_stable_coefficient=s), wind=w) for (f, s, w) in grid_values]
end

function load_forcings(years, settings_list)
    winds = unique(getfield.(settings_list, :wind))
    Dict((y, w) => trial_forcing(y, (wind=w,)) for y in years for w in winds)
end

function write_scores(path, rows)
    open(path, "w") do io
        println(io, "label,bias_wm2,rmse_wm2,low_bias_wm2,low_rmse_wm2,high_bias_wm2,high_rmse_wm2,albedo_bias")
        for (s, r) in rows
            println(io, join((s.label, r.bias, r.rmse, r.low_bias, r.low_rmse, r.high_bias, r.high_rmse, r.albedo_bias), ','))
        end
    end
end

function calibration_main()
    out = joinpath(CAL.out, String(CAL.mode))
    mkpath(out)
    base_loaded = load_mar(daily_file(first(CAL.cal_years)))
    grid = base_loaded.grid
    area = read_mar_static_columns(daily_file(first(CAL.cal_years)), "AREA", grid)
    height = read_mar_static_columns(daily_file(first(CAL.cal_years)), "SH", grid)
    states = prescribed_states(grid)
    CAL.mode == :states && return println("Starting states cached in ", CAL.state_cache)
    CAL.mode == :diagnose && return diagnose_main(grid, area, height, states)
    candidates = settings()
    # Optional cap, for quick tests of the harness.
    max_trials = parse(Int, get(ENV, "CHION_CALIBRATION_MAX_TRIALS", "0"))
    max_trials > 0 && (candidates = candidates[1:min(end, max_trials)])

    # Calibration years: cache forcing once per year and wind choice.
    cache = load_forcings(CAL.cal_years, candidates)
    rows = Tuple{Any,Any}[]
    for (i, s) in enumerate(candidates)
        forcings = Dict(y => cache[(y, s.wind)] for y in CAL.cal_years)
        r = evaluate(s, CAL.cal_years, states[first(CAL.cal_years)], grid, area, height, forcings)
        push!(rows, (s, r))
        @info "Calibration trial" i n=length(candidates) s.label r.bias r.rmse
        write_scores(joinpath(out, "calibration_scores_$(first(CAL.cal_years))_$(last(CAL.cal_years)).csv"), sort(rows; by=x -> x[2].rmse))
    end
    cache = nothing; GC.gc()

    # Validation years for the best sets and the current default.
    best = first.(sort(rows; by=x -> x[2].rmse)[1:min(CAL.n_validate, length(rows))])
    default = (label="current default", model=CAL.mode == :albedo ? (albedo=:dynamic,) : (;), wind=:constant)
    validation = [default; best]
    cache = load_forcings(CAL.val_years, validation)
    vrows = Tuple{Any,Any}[]
    for s in validation
        forcings = Dict(y => cache[(y, s.wind)] for y in CAL.val_years)
        push!(vrows, (s, evaluate(s, CAL.val_years, states[first(CAL.val_years)], grid, area, height, forcings)))
        @info "Validation trial" s.label vrows[end][2].rmse
    end
    write_scores(joinpath(out, "validation_scores_$(first(CAL.val_years))_$(last(CAL.val_years)).csv"), vrows)
    open(joinpath(out, "configuration.csv"), "w") do io
        println(io, "setting,value\ncalibration,$(CAL.mode)\ncalibration_years,$(CAL.cal_years)\nvalidation_years,$(CAL.val_years)")
        println(io, "scored_months,$(CAL.months)")
        write_model_configuration(io, calibration_model(grid))
        println(io, "metric,$(CAL.mode == :albedo ? "absorbed SW SWD*(1-alpha) vs MAR SWD-SWU" : "sensible heat flux vs MAR SHF"), monthly means, area weighted, W/m2")
    end
end

# ---------- flux diagnosis and coupled calibration of the parameterized SEB ----------
# CHION_DIAGNOSE_SET selects the free-running configurations (all from the same
# prescribed starting state):
#   reference  prescribed, and parameterized SEB with MAR albedo (graybody and cloud-proxy LW)
#   sensible   cloud-proxy LW, MAR albedo, grid of SEMIX sensible-heat settings
#   sensible_wide  as sensible, stronger stable damping and constant wind only
#   albedo     cloud-proxy LW, sensible heat from the model defaults (or CHION_SEMIX_*,
#              CHION_WIND_SOURCE), grid of dynamic-albedo settings (fully parameterized)
const DIAG_TERMS = (:sw_abs, :lw_net, :sh, :lh, :melt, :runoff, :refreezing, :smb)
const DIAG_MASS_TERMS = (:melt, :runoff, :refreezing, :smb)
const DIAG_BANDS = ((0, 1000), (1000, 1500), (1500, 2000), (2000, 4000))
const DIAG_SET = Symbol(get(ENV, "CHION_DIAGNOSE_SET", "reference"))

env_wind() = Symbol(get(ENV, "CHION_WIND_SOURCE", "constant"))

function diagnose_configs()
    cloud = (longwave_scheme=:cloud_proxy,)
    if DIAG_SET == :reference
        return [
            (label="prescribed", model=(;), forcing=:prescribed, wind=:constant),
            (label="mar_albedo_graybody", model=(longwave_scheme=:graybody,), forcing=:mar_albedo, wind=:constant),
            (label="mar_albedo_cloud_lw", model=cloud, forcing=:mar_albedo, wind=:constant),
            (label="mar_albedo_cloud_lw_mar_wind", model=cloud, forcing=:mar_albedo, wind=:mar),
        ]
    elseif DIAG_SET == :sensible
        return [(label=@sprintf("factor=%.2f stable=%.0f wind=%s", f, b, w),
                 model=merge(cloud, (semix_sensible_exchange_factor=f, semix_stable_coefficient=b)), forcing=:mar_albedo, wind=w)
                for (f, b, w) in Iterators.product((1.5, 2.0, 2.5, 3.0), (5.0, 10.0, 20.0, 40.0), (:constant, :mar))]
    elseif DIAG_SET == :sensible_wide
        # Stronger stable damping with larger neutral exchange, constant wind only.
        return [(label=@sprintf("factor=%.2f stable=%.0f wind=constant", f, b),
                 model=merge(cloud, (semix_sensible_exchange_factor=f, semix_stable_coefficient=b)), forcing=:mar_albedo, wind=:constant)
                for (f, b) in Iterators.product((2.5, 3.0, 3.5, 4.0, 5.0), (40.0, 80.0, 160.0))]
    elseif DIAG_SET == :albedo
        return [(label=@sprintf("dry=%.2f wet=%.2f ice=%.2f", d, w, i),
                 model=merge(cloud, (albedo=:dynamic, alpha_dry=d, alpha_wet=w, alpha_ice=i)),
                 forcing=:parameterized, wind=env_wind())
                for (d, w, i) in Iterators.product((0.81, 0.84, 0.87), (0.65, 0.70, 0.75, 0.80), (0.40, 0.50)) if w < d]
    end
    error("CHION_DIAGNOSE_SET must be reference, sensible, sensible_wide or albedo.")
end

function diagnose_forcing(year, config)
    path = daily_file(year)
    loaded = load_mar(path)
    base = fully_prescribed_forcing(path, loaded)
    config.forcing == :prescribed && return base, base, loaded.grid
    wind = config.wind == :mar ? read_lowest_wind_level(path, loaded.grid, length(base.time_values)) : base.wind_speed
    # Every surface flux is parameterized; MAR albedo stays prescribed unless the albedo is calibrated.
    return rebuild(base; has_q_sw_net=false, has_q_lw_down=false, has_q_sh=false, has_q_lh=false, wind_speed=wind,
        has_prescribed_albedo=config.forcing == :mar_albedo), base, loaded.grid
end

"""Monthly sums (per column) of Chion's daily-mean fluxes and mass terms next to MAR's."""
function diagnose_year!(model, forcing, reference, path, grid, state, acc)
    ntime = length(forcing.time_values)
    mar = Dict(:melt => read_mar_columns(path, "ME", grid, ntime), :runoff => read_mar_columns(path, "RU", grid, ntime),
        :refreezing => read_mar_columns(path, "RZ", grid, ntime), :smb => read_mar_columns(path, "SMB", grid, ntime))
    sim = Simulation(model; forcing, state=rebind(state, model), years=1, backend=CONFIG.backend,
        write_netcdf=false, compute_year_metrics=false, name="diagnose")
    integ = init_integrator(sim; io=devnull)
    runtime = integ.model_runtime.data
    read_now() = (sw_abs=Array(runtime.workspace.absorbed_shortwave_energy), lw_net=Array(runtime.workspace.net_longwave_energy),
        sh=Array(runtime.workspace.sensible_heat_energy), vapor=Array(runtime.state.vapor_mass), melt=Array(runtime.state.melt),
        runoff=Array(runtime.state.runoff), refreezing=Array(runtime.state.refreezing))
    previous = read_now()
    c = model.parameters.c
    for k in 1:ntime
        step!(integ)
        now = read_now()
        m = month(forcing.time_values[k])
        dt = forcing.dt_days[k] * SECONDS_PER_DAY
        for term in (:sw_abs, :lw_net, :sh)
            acc.chion[term][m, :] .+= (getfield(now, term) .- getfield(previous, term)) ./ dt
        end
        # Latent heat from the booked vapour mass (same conversion as the step).
        acc.chion[:lh][m, :] .+= (now.vapor .- previous.vapor) .* (c.Lv + c.Lm) ./ dt
        for term in (:melt, :runoff, :refreezing)
            acc.chion[term][m, :] .+= getfield(now, term) .- getfield(previous, term)
        end
        # Surface SMB as in the monthly output: precipitation − runoff + vapour.
        precipitation = (@view(forcing.snowfall_rate[:, k]) .+ @view(forcing.rainfall_rate[:, k])) .* dt
        acc.chion[:smb][m, :] .+= precipitation .- (now.runoff .- previous.runoff) .+ (now.vapor .- previous.vapor)
        for term in DIAG_MASS_TERMS
            acc.mar[term][m, :] .+= @view mar[term][:, k]
        end
        acc.mar[:sw_abs][m, :] .+= @view reference.q_sw_net[:, k]
        # MAR has no LWU, so its net longwave stays NaN; compare configs instead.
        acc.mar[:lw_net][m, :] .= NaN
        acc.mar[:sh][m, :] .+= @view reference.q_sh[:, k]
        acc.mar[:lh][m, :] .+= @view reference.q_lh[:, k]
        acc.days[m, :] .+= 1
        previous = now
    end
    finalize!(integ)
    return sim.now
end

function diagnose_main(grid, area, height, states)
    out = joinpath(CAL.out, "diagnose")
    mkpath(out)
    nyears = length(CAL.cal_years)
    tag = "$(DIAG_SET)_$(first(CAL.cal_years))_$(last(CAL.cal_years))"
    summary = []
    configs = diagnose_configs()
    max_trials = parse(Int, get(ENV, "CHION_CALIBRATION_MAX_TRIALS", "0"))
    max_trials > 0 && (configs = configs[1:min(end, max_trials)])
    open(joinpath(out, "flux_bands_$(tag).csv"), "w") do io
        println(io, "config,season,band,term,chion,mar,chion_minus_mar,units")
        for (i, config) in enumerate(configs)
            model = calibration_model(grid; config.model...)
            ncol = length(area)
            acc = (chion=Dict(t => zeros(12, ncol) for t in DIAG_TERMS), mar=Dict(t => zeros(12, ncol) for t in DIAG_TERMS),
                days=zeros(12, ncol))
            state = states[first(CAL.cal_years)]
            for year in CAL.cal_years
                forcing, reference, loaded_grid = diagnose_forcing(year, config)
                state = diagnose_year!(model, forcing, reference, daily_file(year), loaded_grid, state, acc)
                GC.gc()
            end
            row = Dict{String,Float64}()
            for (season, months) in (("JJA", 6:8), ("year", 1:12)), band in DIAG_BANDS, term in DIAG_TERMS
                sel = (height .>= band[1]) .& (height .< band[2])
                w = area[sel]
                if term in DIAG_MASS_TERMS
                    # Gt/yr summed over the season and band.
                    val = a -> sum(sum(a[months, sel]; dims=1)[:] .* w) * 1e-6 / nyears
                    units = "Gt/yr"
                else
                    # Area-weighted mean W/m2 over the season (sums of daily means / days).
                    val = a -> sum(sum(a[months, sel]; dims=1)[:] ./ sum(acc.days[months, sel]; dims=1)[:] .* w) / sum(w)
                    units = "W/m2"
                end
                a, b = val(acc.chion[term]), val(acc.mar[term])
                row["$(season)_$(term)_$(band[1])"] = a - b
                println(io, join((config.label, season, "$(band[1])-$(band[2])", term, a, b, a - b, units), ','))
            end
            flush(io)
            push!(summary, (config.label, row))
            @info "Diagnosis" i n=length(configs) config.label melt=sum(row["year_melt_$(b[1])"] for b in DIAG_BANDS)
            write_diagnose_summary(joinpath(out, "summary_$(tag).csv"), summary)
        end
    end
    println("Diagnosis written to ", out)
end

"""One row per configuration, ranked by the RMS of the annual melt and SMB errors over the elevation bands."""
function write_diagnose_summary(path, summary)
    bands = first.(DIAG_BANDS)
    score(row) = sqrt(sum(row["year_melt_$b"]^2 + row["year_smb_$b"]^2 for b in bands) / (2 * length(bands)))
    open(path, "w") do io
        cols = vcat(["year_$(t)_total" for t in DIAG_MASS_TERMS], ["year_melt_$b" for b in bands], ["year_smb_$b" for b in bands],
            ["JJA_sh_$b" for b in bands], ["JJA_sw_abs_$b" for b in bands])
        println(io, join(vcat(["config", "score_gt"], cols), ','))
        for (label, row) in sort(summary; by=x -> score(x[2]))
            values = [startswith(c, "year_") && endswith(c, "_total") ?
                      sum(row["year_$(split(c, '_')[2])_$b"] for b in bands) : row[c] for c in cols]
            println(io, join(vcat([label, string(score(row))], string.(values)), ','))
        end
    end
end

calibration_main()
