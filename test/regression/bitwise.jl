# Bitwise regression harness for refactors that must not change results.
#
# Runs a matrix of model configurations (every scheme, prescribed and parameterized surface
# fluxes, Float64 and Float32, PDD and ITM) on synthetic forcing, and optionally one or two
# Greenland years on MAR forcing with outputs, then saves (`save`) or compares (`compare`) every
# state array bit for bit.
#
#   julia --project=. --threads 16 test/regression/bitwise.jl save    output/regression/base [cpu|gpu] [synthetic,gris]
#   julia --project=. --threads 16 test/regression/bitwise.jl compare output/regression/base [cpu|gpu] [synthetic,gris]
#
# Not part of `Pkg.test()`: the Greenland part reads MAR files and takes a few minutes.
using Chion
using Dates
using NCDatasets
using Printf
using Serialization

const MODE = get(ARGS, 1, "compare")
const DIR = get(ARGS, 2, "output/regression/base")
const BACKEND = Symbol(get(ARGS, 3, "cpu"))
const PARTS = split(get(ARGS, 4, "synthetic,gris"), ",")
const MAR = "/p/projects/ou/labs/ai/Nils/MAR_daily/mar_daily/MARv3.14.3-10km-daily-ERA5-{year}.nc"
BACKEND === :gpu && @eval using CUDA
backend() = Chion._ka_backend_from(BACKEND === :gpu ? :gpu : :threads)

# ---------------------------------------------------------------------------
# Synthetic forcing: 48 columns from cold accumulation to warm ablation, seasonal cycle
# ---------------------------------------------------------------------------
const NCOL = 48

function synthetic_forcing(; years = 3, dt_days = 1.0, prescribed = false, extra = (;))
    nstep = round(Int, years * 365 / dt_days)
    t = [(k - 0.5) * dt_days for k in 1:nstep]
    col = 1:NCOL
    latitude = [60.0 + 20.0 * (c - 1) / (NCOL - 1) for c in col]
    offset = [-22.0 + 30.0 * (c - 1) / (NCOL - 1) for c in col]          # °C
    season(d) = -cos(2π * d / 365)                                          # −1 winter, +1 summer
    weather(c, k) = 3.0 * sin(0.37 * k + 1.3 * c) + 1.5 * sin(2.1 * k + 0.4 * c)
    T = [offset[c] + 12.0 * season(t[k]) + weather(c, k) for c in col, k in 1:nstep]
    precip = [max(0.0, 2.5 + 2.0 * sin(0.21 * k + c) + 0.03 * c) for c in col, k in 1:nstep]
    snow = [T[c, k] < 1.0 ? precip[c, k] : 0.0 for c in col, k in 1:nstep]
    rain = [T[c, k] < 1.0 ? 0.0 : precip[c, k] for c in col, k in 1:nstep]
    sw = [max(0.0, 180.0 + 160.0 * season(t[k]) - 2.0 * (latitude[c] - 60.0) + 25.0 * sin(0.9 * k + c))
          for c in col, k in 1:nstep]
    wind = [3.0 + 2.0 * abs(sin(0.13 * k + 0.7 * c)) for c in col, k in 1:nstep]
    time_values = [DateTime(1990, 1, 1) + Millisecond(round(Int, (k - 0.5) * dt_days * 86_400_000)) for k in 1:nstep]
    kw = (; dt_days = fill(dt_days, nstep), ncol = NCOL, air_temperature_c = T, snowfall_mm_day = snow,
        rainfall_mm_day = rain, shortwave_down = sw, wind_speed = wind, latitude_deg = latitude,
        surface_height = [200.0 + 60.0 * c for c in col], time_values)
    if prescribed
        lw = [250.0 + 2.5 * T[c, k] + 20.0 * sin(0.5 * k + c) for c in col, k in 1:nstep]
        sh = [5.0 + 10.0 * sin(0.3 * k + 0.2c) for c in col, k in 1:nstep]
        lh = [-3.0 + 6.0 * sin(0.17 * k + 0.5c) for c in col, k in 1:nstep]
        kw = merge(kw, (; q_lw_down = lw, q_sh = sh, q_lh = lh, q_sw_net = 0.25 .* sw))
    end
    return SnowpackForcing(; merge(kw, extra)...)
end

function prescribed_albedo_forcing()
    f = synthetic_forcing()
    nstep = length(f.dt_days)
    alb = [0.55 + 0.3 * abs(sin(0.05 * k + c)) for c in 1:NCOL, k in 1:nstep]
    return synthetic_forcing(; extra = (; prescribed_albedo = alb))
end

function semix_forcing()
    nstep = 3 * 365
    extra = (; cloud = [0.5 + 0.4 * sin(0.2k + c) for c in 1:NCOL, k in 1:nstep],
        dust_deposition = fill(1e-9, NCOL, nstep), z_sur_std = fill(20.0, NCOL, nstep))
    return synthetic_forcing(; extra)
end

# Each case: (params keywords, layering keywords, number type, forcing).
const BESSI_CASES = [
    "default" => (() -> (;), (;), Float64, synthetic_forcing),
    "float32" => (() -> (;), (;), Float32, synthetic_forcing),
    "albedo_constant" => (() -> (; albedo_scheme = ConstantAlbedo()), (;), Float64, synthetic_forcing),
    "albedo_aging" => (() -> (; albedo_scheme = AgingAlbedo()), (;), Float64, synthetic_forcing),
    "albedo_semix_dang" => (() -> (; albedo_scheme = SemixAlbedo()), (;), Float64, semix_forcing),
    "albedo_semix_ww" => (() -> (; albedo_scheme = SemixAlbedo(snow_bands = WarrenWiscombeSnowAlbedo())), (;), Float64, semix_forcing),
    "albedo_prescribed" => (() -> (; albedo_scheme = PrescribedAlbedo()), (;), Float64, prescribed_albedo_forcing),
    "seb_bessi" => (() -> (; surface_energy_balance_scheme = BESSISurfaceEnergyBalance()), (;), Float64, synthetic_forcing),
    "turbulence_bessi" => (() -> (; turbulent_flux_scheme = BESSITurbulentFluxes()), (;), Float64, synthetic_forcing),
    "seb_turbulence_bessi" => (() -> (; surface_energy_balance_scheme = BESSISurfaceEnergyBalance(),
        turbulent_flux_scheme = BESSITurbulentFluxes()), (;), Float64, synthetic_forcing),
    "longwave_graybody" => (() -> (; longwave_scheme = GraybodyLongwave()), (;), Float64, synthetic_forcing),
    "densification_htessel" => (() -> (; densification_scheme = HTESSELDensification()), (;), Float64, synthetic_forcing),
    "fresh_snow_parameterized" => (() -> (; fresh_snow_density_scheme = ParameterizedFreshSnowDensity()), (;), Float64, synthetic_forcing),
    "no_diurnal" => (() -> (; diurnal_cycle = NoDiurnalCycle()), (;), Float64, synthetic_forcing),
    "diurnal_sw_only" => (() -> (; diurnal_cycle = DiurnalSubsteps(temperature_cycle = NoDiurnalTemperatureCycle())),
        (;), Float64, synthetic_forcing),
    "no_substrate_coarse" => (() -> (;), (; ice_substrate_layers = 0, near_surface_layer_max_thicknesses_m = (Inf, Inf, Inf, Inf)),
        Float64, synthetic_forcing),
    "prescribed_fluxes" => (() -> (;), (;), Float64, () -> synthetic_forcing(; prescribed = true)),
    "prescribed_fluxes_bessi" => (() -> (; surface_energy_balance_scheme = BESSISurfaceEnergyBalance(),
        turbulent_flux_scheme = BESSITurbulentFluxes()), (;), Float64, () -> synthetic_forcing(; prescribed = true)),
    "three_hourly" => (() -> (;), (;), Float64, () -> synthetic_forcing(; years = 1, dt_days = 0.125)),
    "float32_semix_htessel" => (() -> (; albedo_scheme = SemixAlbedo(), densification_scheme = HTESSELDensification()),
        (;), Float32, semix_forcing),
]

const BULK_CASES = [
    "pdd_simple" => (PDDModel, () -> PDDParams(), Float64),
    "pdd_pism" => (PDDModel, () -> PDDParams(method = PISMPDD()), Float64),
    "pdd_float32" => (PDDModel, () -> PDDParams(NF = Float32), Float32),
    "itm" => (ITMModel, () -> ITMParams(), Float64),
]

# Every array of the host state and of the cache diagnostics, by name.
function state_arrays(model)
    update_diagnostics!(model)
    state = cpu_state(model)
    out = Dict{Symbol, Any}(name => Array(getfield(state, name)) for name in fieldnames(typeof(state)))
    if model isa BESSIModel
        for name in (:thickness, :wet_mass, :bulk_density, :liquid_water)
            out[name] = Array(getproperty(model.cache, name))
        end
    end
    return out
end

# The whole synthetic forcing twice, year by year (the forcing of the baseline ran twice in one piece).
twice(model) = repeat(forcing_years(model.forcing.source), 2)

function run_bessi(params_kwargs, layering_kwargs, NF, forcing)
    grid = Grid(NCOL; NF, layering = Layering(; layering_kwargs...))
    model = BESSIModel(grid; forcing, params = Params(; NF, params_kwargs...), backend = backend())
    run!(Simulation(model; years = twice(model), verbose = false))
    return state_arrays(model)
end

function run_bulk(M, params, NF)
    forcing = synthetic_forcing(; extra = (; ice_thickness = [c <= NCOL ÷ 2 ? 0.0 : 1500.0 for c in 1:NCOL],
        annual_pdd = [max(0.0, 40.0 * (c - 20)) for c in 1:NCOL]))
    model = M(Grid(NCOL; NF); forcing, params, backend = backend())
    run!(Simulation(model; years = twice(model), verbose = false))
    return state_arrays(model)
end

# ---------------------------------------------------------------------------
# Greenland on MAR: two years with outputs, parameterized and prescribed surface fluxes
# ---------------------------------------------------------------------------
function run_gris(label, surface_fluxes, NF, outdir)
    forcing = Forcing(MAR; surface_fluxes, years = 1980:1981)
    model = BESSIModel(forcing; NF, backend = backend())
    mkpath(outdir)
    outputs = [
        Output(joinpath(outdir, "$(label)_monthly.nc"); variables = :monthly, frequency = :monthly),
        Output(joinpath(outdir, "$(label)_budget.csv"); variables = :budget, frequency = :yearly, aggregate = :area_sum),
        Output(joinpath(outdir, "$(label)_firn.nc"); variables = :firn, frequency = Month(12), kind = :snapshot),
    ]
    run!(Simulation(model; years = 1980:1981, output = outputs, name = label, verbose = false))
    return state_arrays(model)
end

# ---------------------------------------------------------------------------
# Compare
# ---------------------------------------------------------------------------
bits(x::AbstractArray{Float64}) = reinterpret(UInt64, x)
bits(x::AbstractArray{Float32}) = reinterpret(UInt32, x)
bits(x) = x

function compare_arrays(case, base, new)
    ok = true
    for name in sort!(collect(union(keys(base), keys(new))))
        if !haskey(base, name) || !haskey(new, name)
            println("  $case.$name: present only in ", haskey(base, name) ? "baseline" : "new run")
            ok = false
            continue
        end
        a, b = base[name], new[name]
        if size(a) != size(b) || eltype(a) != eltype(b)
            println("  $case.$name: size/type $(size(a)) $(eltype(a)) vs $(size(b)) $(eltype(b))")
            ok = false
        elseif bits(a) != bits(b)
            n = count(bits(a) .!= bits(b))
            d = eltype(a) <: AbstractFloat ? maximum(abs.(filter(isfinite, Float64.(a) .- Float64.(b))); init = 0.0) : NaN
            @printf("  %s.%s: %d of %d values differ (max abs diff %.3e)\n", case, name, n, length(a), d)
            ok = false
        end
    end
    return ok
end

function compare_files(base_dir, new_dir)
    ok = true
    for f in sort(readdir(base_dir))
        endswith(f, ".jls") && continue
        a, b = joinpath(base_dir, f), joinpath(new_dir, f)
        if !isfile(b)
            println("  file $f missing"); ok = false; continue
        end
        if endswith(f, ".csv")
            # The wall-time-free CSV must match exactly.
            read(a, String) == read(b, String) || (println("  file $f differs"); ok = false)
        elseif endswith(f, ".nc")
            NCDataset(a) do da
                NCDataset(b) do db
                    for v in keys(da)
                        haskey(db, v) || (println("  $f:$v missing"); ok = false; continue)
                        x, y = Array(da[v].var[:]), Array(db[v].var[:])
                        isequal(x, y) || (println("  $f:$v differs"); ok = false)
                    end
                end
            end
        end
    end
    return ok
end

function main()
    root = joinpath(DIR, String(BACKEND))
    newroot = MODE == "save" ? root : joinpath(DIR * "_new", String(BACKEND))
    mkpath(root); mkpath(newroot)
    results = Dict{String, Any}()
    t0 = time()
    if "synthetic" in PARTS
        for (name, (kw, layering, NF, f)) in BESSI_CASES
            results[name] = run_bessi(kw(), layering, NF, f())
            @printf("ran %-26s %6.1f s\n", name, time() - t0); flush(stdout)
        end
        for (name, (M, params, NF)) in BULK_CASES
            results[name] = run_bulk(M, params(), NF)
            @printf("ran %-26s %6.1f s\n", name, time() - t0); flush(stdout)
        end
    end
    if "gris" in PARTS
        for (label, mode, nf) in (("gris_param", :parameterized, Float64), ("gris_prescribed", :prescribed, Float64),
                                  ("gris_param_f32", :parameterized, Float32))
            results[label] = run_gris(label, mode, nf, newroot)
            @printf("ran %-26s %6.1f s\n", label, time() - t0); flush(stdout)
        end
    end
    if MODE == "save"
        for (name, arrays) in results
            serialize(joinpath(root, name * ".jls"), arrays)
        end
        println("saved $(length(results)) cases to $root")
        return
    end
    failed = String[]
    for (name, arrays) in sort!(collect(results); by = first)
        path = joinpath(root, name * ".jls")
        isfile(path) || (println("  no baseline for $name"); push!(failed, name); continue)
        compare_arrays(name, deserialize(path), arrays) || push!(failed, name)
    end
    "gris" in PARTS && (compare_files(root, newroot) || push!(failed, "gris output files"))
    if isempty(failed)
        println("BITWISE IDENTICAL: $(length(results)) cases ($BACKEND)")
    else
        println("DIFFERENCES in: ", join(failed, ", "))
        exit(1)
    end
end

main()
