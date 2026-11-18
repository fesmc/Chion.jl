using Test
using Dates
using NCDatasets
using Chion

"""Daily forcing for four columns on a 2×2 grid, with warm summers."""
function _experiment_forcing(years = 2001:2002)
    times = [DateTime(y, 1, 1, 12) + Day(d) for y in years for d in 0:(Dates.daysinyear(y) - 1)]
    n = length(times)
    doy = dayofyear.(times)
    temperature = [-12.0 + 14.0 * sin(2π * (d - 100) / 365) + c for c in 0:3, d in doy]
    forcing = SnowpackForcing(dt_days = ones(n), ncol = 4, time_values = times,
        air_temperature_c = temperature, snowfall_mm_day = 2.0, rainfall_mm_day = 0.2,
        shortwave_down = [max(0.0, 250.0 * sin(2π * (d - 80) / 365)) for _ in 1:4, d in doy],
        wind_speed = 4.0, latitude_deg = 70.0)
    grid = Grid(4; x = [0.0, 1.0], y = [0.0, 1.0], js = [1, 1, 2, 2], is = [1, 2, 1, 2], mask = ones(2, 2))
    return forcing, grid
end

_csv(path) = (lines = readlines(path); (split(lines[1], ','), [parse.(Float64, split(l, ',')[3:end]) for l in lines[2:end]]))
_layering8() = Layering(Ntot = 8)

@testset "Simulation outputs" begin
    snowpack_forcing, grid = _experiment_forcing()
    forcing = Forcing(snowpack_forcing; grid, area = [1.0, 2.0, 3.0, 4.0])
    @test forcing_years(forcing) == [2001, 2002]
    model = BESSIModel(forcing; layering = _layering8())
    dir = mktempdir()
    out(name; kw...) = Output(joinpath(dir, name); kw...)
    sim = Simulation(model; spinup = CycledSpinup(2001:2001; cycles = 1), verbose = false, output = [
        out("daily.csv"; variables = :budget, frequency = :daily, aggregate = :area_sum),
        out("monthly.csv"; variables = :budget, frequency = :monthly, aggregate = :area_sum),
        out("yearly.csv"; variables = [:budget; :swe; :albedo], frequency = :yearly, aggregate = :area_sum),
        out("monthly.nc"; variables = [:melt, :swe, :albedo], frequency = :monthly),
        out("snap.nc"; variables = [:swe, :density], frequency = Month(6), kind = :snapshot),
        out("regions.csv"; variables = [:melt, :swe], frequency = :yearly, aggregate = :area_sum,
            regions = (west = [true, false, true, false], east = [false, true, false, true]))])
    @test run!(sim) === sim
    @test sim.years == [2001, 2002] && finished(sim)

    # Period totals of fluxes add up across frequencies.
    hd, daily = _csv(joinpath(dir, "daily.csv"))
    _, monthly = _csv(joinpath(dir, "monthly.csv"))
    hy, yearly = _csv(joinpath(dir, "yearly.csv"))
    @test length(daily) == 365 + 365 && length(monthly) == 24 && length(yearly) == 2
    @test hd[3] == "surface_smb (Gt)"
    for j in 1:6
        @test sum(r[j] for r in daily[1:365]) ≈ yearly[1][j] rtol = 1e-12
        @test sum(r[j] for r in monthly[13:24]) ≈ yearly[2][j] rtol = 1e-12
    end
    # SMB = precipitation − runoff + vapour exchange.
    for r in yearly
        @test r[1] ≈ r[2] - r[4] + r[6] rtol = 1e-12
    end
    @test yearly[1][3] > 0                       # some melt in the warm columns

    # Regions add up to the total.
    _, reg = _csv(joinpath(dir, "regions.csv"))
    @test reg[2][1] + reg[2][2] ≈ yearly[2][3] rtol = 1e-12

    # Gridded output: monthly melt sums to the area-summed yearly melt.
    ds = NCDataset(joinpath(dir, "monthly.nc"))
    @test size(ds["melt"]) == (2, 2, 24)
    melt = Float64.(ds["melt"][:, :, 13:24])
    area = [1.0 3.0; 2.0 4.0]                     # (x, y) of columns 1–4
    @test sum(melt .* area) * 1e-6 ≈ yearly[2][3] rtol = 1e-6
    @test ds["time"][1] == DateTime(2001, 1, 1, 12)
    close(ds)
    ds = NCDataset(joinpath(dir, "snap.nc"))
    @test size(ds["density"]) == (2, 2, 8, 4)
    close(ds)

    # The final state equals a run without any output.
    plain = BESSIModel(forcing; layering = _layering8())
    run!(Simulation(plain; spinup = CycledSpinup(2001:2001; cycles = 1), verbose = false))
    @test plain.state.mass == model.state.mass
    @test plain.state.runoff == model.state.runoff
    # A finished simulation cannot run on.
    @test_throws ArgumentError time_step!(sim)
end

@testset "run! resumes; time_step! and integrate! step the same" begin
    snowpack_forcing, grid = _experiment_forcing()
    forcing = Forcing(snowpack_forcing; grid)
    whole = BESSIModel(forcing; layering = _layering8())
    run!(Simulation(whole; verbose = false))

    parts = BESSIModel(forcing; layering = _layering8())
    sim = Simulation(parts; verbose = false)
    run!(sim; years = 1)
    @test !finished(sim) && sim.clock.year_index == 2 && sim.clock.step == 1
    @test sim.clock.iteration == 365
    run!(sim; years = 1)
    @test finished(sim)
    @test parts.state.mass == whole.state.mass && parts.state.Tsrf == whole.state.Tsrf

    steps = BESSIModel(forcing; layering = _layering8())
    sim = Simulation(steps; verbose = false)
    while !finished(sim)
        time_step!(sim)
    end
    @test sim.clock.iteration == 730
    @test steps.state.mass == whole.state.mass && steps.state.runoff == whole.state.runoff

    # integrate! steps the loaded year without a simulation.
    direct = BESSIModel(forcing; layering = _layering8())
    integrate!(direct, 1:365)
    update_forcing!(direct, 2002)
    integrate!(direct, 1:365)
    @test direct.state.mass == whole.state.mass
    @test_throws ArgumentError integrate!(direct, 1:400)
end

@testset "Simulation and model construction" begin
    snowpack_forcing, grid = _experiment_forcing(2001:2001)
    forcing = Forcing(snowpack_forcing; grid)
    model = BESSIModel(forcing)
    @test_throws ArgumentError Simulation(model; years = [2001, 2003])           # no forcing for 2003
    @test_throws ArgumentError Simulation(model; spinup = CycledSpinup([1999]))
    @test_throws ArgumentError Simulation(model; output = [Output("a.nc"), Output("a.nc")])
    @test_throws ArgumentError BESSIModel(Grid(3); forcing)                     # column count
    @test_throws ArgumentError BESSIModel(Grid(4; NF = Float32); forcing, params = Params())   # number type
    @test_throws ArgumentError Params(no_such_field = 1.0)
    @test_throws ArgumentError Params(albedo_scheme = AgingAlbedo(), alpha_wet = 0.9, alpha_dry = 0.8)
    # Flat forwarding: state, cache, params, grid.
    @test model.mass === model.state.mass
    @test model.thickness === model.cache.thickness
    @test model.alpha_dry === model.params.alpha_dry
    @test model.Ntot == model.grid.layering.Ntot
    @test model.year == 2001
    @test_throws ErrorException model.no_such_property
    # to_backend on the CPU copies the model.
    copied = to_backend(model, CPU())
    @test copied.state.mass == model.state.mass && copied.state.mass !== model.state.mass
    # A restart: the state of one model continues in another.
    run!(Simulation(model; verbose = false))
    continued = BESSIModel(forcing; state = cpu_state(model))
    @test continued.state.mass == model.state.mass
end

@testset "Run directory, metadata and restarts" begin
    snowpack_forcing, grid = _experiment_forcing()
    forcing = Forcing(snowpack_forcing; grid, area = fill(10.0, 4))
    params = Params(albedo_scheme = AgingAlbedo(cold_timescale_days = 12), alpha_dry = 0.83)
    dir = joinpath(mktempdir(), "run")
    model = BESSIModel(forcing; layering = _layering8(), params)
    sim = Simulation(model; name = "dirtest", verbose = false, output = OutputConfig(; dir, restart_every = 1,
        outputs = [Output("budget.csv"; variables = :budget, frequency = :yearly, aggregate = :area_sum)]))
    @test sim.output.outputs[1].path == joinpath(dir, "budget.csv")
    run!(sim)
    @test isfile(joinpath(dir, "budget.csv"))
    @test isfile(joinpath(dir, "restart_2001.jls")) && isfile(joinpath(dir, "restart_2002.jls"))
    log = read(joinpath(dir, "log.txt"), String)
    @test occursin("year 2001 done", log) && occursin("year 2002 done", log) && occursin("restart written", log)
    # The metadata is a run_config file: the same grid and parameters back.
    meta = Chion.TOML.parsefile(joinpath(dir, "chion_run.toml"))
    @test meta["name"] == "dirtest" && meta["years"] == [2001, 2002] && meta["columns"] == 4
    M, back = Chion._config_params(meta["model"])
    @test M === BESSIModel && back == params
    @test Chion._config_layering(meta["grid"]) == model.grid.layering

    # A restart after the first year continues to the state of the whole run.
    continued = BESSIModel(forcing; layering = _layering8(), params)
    sim2 = Simulation(continued; years = [2002], restart = joinpath(dir, "restart_2001.jls"), verbose = false)
    @test sim2.io.spun_up
    run!(sim2)
    @test continued.state.mass == model.state.mass && continued.state.runoff == model.state.runoff
    @test_throws ArgumentError Simulation(PDDModel(forcing); restart = joinpath(dir, "restart_2001.jls"))
end

@testset "Forcing edited in place and backend names" begin
    snowpack_forcing, grid = _experiment_forcing(2001:2001)
    edited = BESSIModel(Forcing(snowpack_forcing; grid); backend = :threads)
    replaced = BESSIModel(Forcing(deepcopy(snowpack_forcing); grid))
    data = edited.forcing.data
    data.air_temperature .+= 1.5                           # e.g. a coupled model's warming
    sync_forcing!(edited)
    @test edited.forcing.fields.air_temperature == data.air_temperature
    set_input!(replaced, :air_temperature, replaced.forcing.data.air_temperature .+ 1.5)
    integrate!(edited, 1:365); integrate!(replaced, 1:365)
    @test edited.state.mass == replaced.state.mass
    @test_throws ArgumentError BESSIModel(Grid(4); forcing = snowpack_forcing, backend = :tpu)
end

@testset "DebugConfig catches non-finite states" begin
    snowpack_forcing, grid = _experiment_forcing(2001:2001)
    forcing = Forcing(snowpack_forcing; grid)
    sim = Simulation(BESSIModel(forcing); verbose = false, debug = DebugConfig(check_nans = true))
    time_step!(sim)                                     # finite: no error
    sim.model.state.Tsrf[3] = NaN
    err = try time_step!(sim); nothing; catch e; e; end
    @test err isa ErrorException && occursin("Non-finite", err.msg) && occursin("in column 3 at 2001-01-02", err.msg)
    # Inactive columns are not checked; checks off by default.
    model = BESSIModel(forcing)
    set_active_mask!(model, [true, true, false, true])
    model.state.Tsrf[3] = NaN
    time_step!(Simulation(model; verbose = false, debug = DebugConfig(check_nans = true)))
    @test isnan(model.state.Tsrf[3])
    model.state.Tsrf[1] = NaN
    @test time_step!(Simulation(model; verbose = false)) isa Simulation
end

@testset "Active columns" begin
    snowpack_forcing, grid = _experiment_forcing(2001:2001)
    forcing = Forcing(snowpack_forcing; grid)
    model = BESSIModel(forcing; layering = _layering8())
    sim = Simulation(model; verbose = false)
    set_active_mask!(sim, [true, false, true, true])
    run!(sim)
    @test model.state.smb_ice[2] == 0 && model.state.N[2] == 0
    @test all(model.state.smb_ice[[1, 3, 4]] .!= 0)
    # A 2-D mask on the grid; the newly inactive column is reset.
    set_active_mask!(model, [1 0; 1 1])                 # (y, x): column 2 (x 2, y 1) inactive
    @test model.cache.active == [true, false, true, true]
    set_active_mask!(model, [true, true, false, true])
    @test model.state.N[3] == 0 && model.state.smb_ice[3] == 0
    @test_throws ArgumentError set_active_mask!(model, falses(4))
end

@testset "ClimatologySpinup and yearly_step!" begin
    snowpack_forcing, grid = _experiment_forcing(2001:2002)
    forcing = Forcing(snowpack_forcing; grid)
    climatology = Forcing(slice_forcing(snowpack_forcing, 1:365); grid)
    model = BESSIModel(forcing; layering = _layering8())
    sim = Simulation(model; spinup = ClimatologySpinup(climatology; years = 2), verbose = false)
    coupling = yearly_step!(sim)
    @test coupling.year == 2001
    reference = BESSIModel(forcing; layering = _layering8())
    integrate!(reference, 1:365)                    # spin-up: the climatology (= 2001) twice
    integrate!(reference, 1:365)
    smb_before = copy(reference.state.smb_ice)
    integrate!(reference, 1:365)
    @test coupling.smb_ice_yearly ≈ reference.state.smb_ice .- smb_before
    @test size(coupling.smb_ice_yearly_grid) == (2, 2)
    @test all(200 .< coupling.mean_T_srf_K .< 274)
end

@testset "Simulation configuration (TOML)" begin
    snowpack_forcing, grid = _experiment_forcing(2001:2001)
    dir = mktempdir()
    path = joinpath(dir, "forcing.nc")
    NCDataset(path, "c") do ds
        defDim(ds, "x", 2); defDim(ds, "y", 2); defDim(ds, "time", 365)
        defVar(ds, "x", [0.0, 1.0], ("x",)); defVar(ds, "y", [0.0, 1.0], ("y",))
        defVar(ds, "time", snowpack_forcing.time_values, ("time",))
        field(f) = permutedims(reshape(f, 2, 2, :), (1, 2, 3))
        defVar(ds, "TT", field(snowpack_forcing.air_temperature .- 273.15), ("x", "y", "time"))
        defVar(ds, "SF", field(snowpack_forcing.snowfall_rate .* 86400), ("x", "y", "time"))
        defVar(ds, "RF", field(snowpack_forcing.rainfall_rate .* 86400), ("x", "y", "time"))
        defVar(ds, "SWD", field(snowpack_forcing.shortwave_down), ("x", "y", "time"))
        defVar(ds, "LAT", fill(70.0, 2, 2), ("x", "y"))
        defVar(ds, "CELL_AREA", fill(100.0, 2, 2), ("x", "y"))
    end
    write(joinpath(dir, "run.toml"), """
        name = "toml_test"
        years = [2001, 2001]
        output_dir = "out"
        [forcing]
        path = "forcing.nc"
        format = "generic"
        area = "CELL_AREA"
        [grid]
        Ntot = 8
        [model]
        alpha_dry = 0.83
        albedo_scheme = {type = "aging", cold_timescale_days = 15}
        [model.diurnal_cycle]
        type = "substeps"
        max_substeps = 4
        temperature_cycle = "none"
        [spinup]
        type = "none"
        [[output]]
        file = "budget.csv"
        variables = "budget"
        frequency = "1 month"
        aggregate = "area_sum"
        """)
    sim = run_config(joinpath(dir, "run.toml"); dry_run = true)
    p = sim.model.params
    @test sim.model.Ntot == 8 && p.alpha_dry == 0.83
    @test p.albedo_scheme == AgingAlbedo(cold_timescale_days = 15.0)
    @test p.diurnal_cycle == DiurnalSubsteps(max_substeps = 4, temperature_cycle = NoDiurnalTemperatureCycle())
    sim = run_config(joinpath(dir, "run.toml"); set = ["model.surface_energy_balance_scheme" => "bessi"])
    @test finished(sim) && sim.model.params.surface_energy_balance_scheme === BESSISurfaceEnergyBalance()
    _, rows = _csv(joinpath(dir, "out", "budget.csv"))
    @test length(rows) == 12
    @test all(r -> r[2] > 0, rows)               # precipitation in every month
    @test isfile(joinpath(dir, "out", "run.toml"))  # the configuration is kept with the results
    # The Antarctic setup and the bulk models by name.
    cfg = Dict{String, Any}("setup" => "AISv1", "number_type" => "Float32")
    M, params = Chion._config_params(cfg)
    @test M === BESSIModel && params == ParamsAISv1(NF = Float32)
    M, params = Chion._config_params(Dict{String, Any}("type" => "pdd", "method" => "pism", "ddf_snow" => 4))
    @test M === PDDModel && params.method === PISMPDD() && params.ddf_snow === 4.0
end

@testset "TOML scheme names" begin
    S = Chion
    @test S._scheme_from_config(S.AbstractAlbedoScheme, "semix") === SemixAlbedo()
    @test S._scheme_from_config(S.AbstractLongwaveScheme, "graybody") === GraybodyLongwave()
    @test S._scheme_from_config(S.AbstractAlbedoScheme,
        Dict{String, Any}("type" => "semix", "snow_bands" => "warren_wiscombe", "frac_vu" => 0.5)) ==
        SemixAlbedo(snow_bands = WarrenWiscombeSnowAlbedo(), frac_vu = 0.5)
    err = try S._scheme_from_config(S.AbstractAlbedoScheme, "nonsense"); catch e; e; end
    @test err isa ArgumentError
    @test occursin("Unknown albedo scheme \"nonsense\"", err.msg)
    @test occursin("\"constant\", \"dynamic\", \"prescribed\", \"aging\", \"semix\"", err.msg)
    @test_throws ArgumentError S._scheme_from_config(S.AbstractAlbedoScheme, Dict{String, Any}("frac_vu" => 0.5))
    @test_throws ArgumentError S._scheme_from_config(S.AbstractAlbedoScheme,
        Dict{String, Any}("type" => "dynamic", "frac_vu" => 0.5))          # no coefficients
    # Every name round-trips.
    for (kind, _, schemes) in S._SCHEME_TABLE, (name, T) in pairs(schemes)
        scheme = S._scheme_from_config(kind, String(name))
        @test scheme isa T
        @test S._scheme_name(scheme) === name
    end
end

@testset "relabel_year" begin
    f, grid = _experiment_forcing(2004:2004)            # a leap year
    g = relabel_year(f, 2005)
    @test year(first(g.time_values)) == 2005 && year(last(g.time_values)) == 2005
    @test count(t -> month(t) == 2 && day(t) == 28, g.time_values) == 2   # 29 February → 28 February
    @test g.day_of_year == f.day_of_year && g.dt_days == f.dt_days
    function run(forcing)
        model = BESSIModel(Grid(grid; layering = _layering8()); forcing = Forcing(forcing; grid))
        run!(Simulation(model; verbose = false))
        return model.state
    end
    a, b = run(f), run(g)
    @test Array(a.mass) == Array(b.mass) && Array(a.runoff) == Array(b.runoff)
end

@testset "PDD and ITM budgets" begin
    f, grid = _experiment_forcing(2001:2001)
    for M in (PDDModel, ITMModel)
        forcing = M === ITMModel ? SnowpackForcing(; time_values = f.time_values, dt_days = f.dt_days,
            air_temperature = f.air_temperature, snowfall_rate = f.snowfall_rate, rainfall_rate = f.rainfall_rate,
            shortwave_down = f.shortwave_down, latitude_deg = f.latitude_deg, surface_height = fill(1500.0, 4),
            ice_thickness = fill(1000.0, 4), annual_pdd = fill(50.0, 4)) : f
        dir = mktempdir()
        variables = [:surface_smb, :precipitation, :runoff]
        model = M(grid; forcing = Forcing(forcing; grid, area = fill(100.0, 4)))
        run!(Simulation(model; verbose = false,
            output = Output(joinpath(dir, "budget.csv"); variables, frequency = :yearly, aggregate = :area_sum)))
        _, rows = _csv(joinpath(dir, "budget.csv"))
        smb, precip, runoff = rows[1]
        @test precip > 0 && runoff >= 0
        @test smb ≈ precip - runoff                      # no vapour exchange in PDD and ITM
        @test haskey(smbstats(model), :total_Gt)
    end
end

@testset "Output choices as types" begin
    o = Output("a.nc"; frequency = :monthly, aggregate = :columns, kind = :snapshot)
    @test o.frequency === Monthly() && o.aggregate === Columns() && o.kind === Snapshot() && o.format === Chion.NetCDFFormat()
    @test Output("b.csv"; frequency = Month(6), aggregate = AreaSum()).frequency === Every(Month(6))
    @test Output("c.nc"; frequency = Yearly(), kind = Mean()).kind === Mean()
    @test Chion._label(Every(Month(6))) == "6 months" && Chion._label(AreaSum()) == "area_sum"
    @test_throws ArgumentError Output("d.nc"; frequency = :hourly)
    @test_throws ArgumentError Output("d.nc"; aggregate = :sum)
    @test_throws ArgumentError Output("d.txt")
    @test_throws ArgumentError Output("d.nc"; frequency = Day(0))
    @test_throws ArgumentError Output("d.nc"; regions = (a = [true],))           # regions need area aggregation
end

@testset "configuration overrides" begin
    @test Chion._set_value("0.84") == 0.84 && Chion._set_value("prescribed") == "prescribed"
    @test Chion._set_value("[1940, 1969]") == [1940, 1969]
    cfg = Chion._apply_overrides!(Dict{String, Any}("model" => Dict{String, Any}("type" => "bessi")),
        ["model.alpha_dry" => "0.84", "spinup.type" => "none", "model.albedo_scheme.type" => "aging"])
    @test cfg["model"]["alpha_dry"] == 0.84 && cfg["spinup"]["type"] == "none"
    @test cfg["model"]["albedo_scheme"]["type"] == "aging"
end
