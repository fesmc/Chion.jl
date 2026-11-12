using Test
using Dates
using NCDatasets
using Chion

"""Two years of daily forcing for four columns on a 2×2 grid, with warm summers."""
function _experiment_forcing(years=2001:2002)
    times = [DateTime(y, 1, 1, 12) + Day(d) for y in years for d in 0:(Dates.daysinyear(y) - 1)]
    n = length(times)
    doy = dayofyear.(times)
    temperature = [-12.0 + 14.0 * sin(2π * (d - 100) / 365) + c for c in 0:3, d in doy]
    forcing = SnowpackForcing(dt_days=ones(n), ncol=4, time_values=times,
        air_temperature_c=temperature, snowfall_mm_day=2.0, rainfall_mm_day=0.2,
        shortwave_down=[max(0.0, 250.0 * sin(2π * (d - 80) / 365)) for _ in 1:4, d in doy],
        wind_speed=4.0, latitude_deg=70.0)
    grid = SnowpackGrid(4; x=[0.0, 1.0], y=[0.0, 1.0], js=[1, 1, 2, 2], is=[1, 2, 1, 2], mask=ones(2, 2))
    return forcing, grid
end

_csv(path) = (lines = readlines(path); (split(lines[1], ','), [parse.(Float64, split(l, ',')[3:end]) for l in lines[2:end]]))

@testset "Experiment outputs" begin
    snowpack_forcing, grid = _experiment_forcing()
    forcing = Forcing(snowpack_forcing; grid, area=[1.0, 2.0, 3.0, 4.0])
    @test forcing_years(forcing) == [2001, 2002]
    model = BESSIModel(forcing; Ntot=8)
    dir = mktempdir()
    out(name; kw...) = Output(joinpath(dir, name); kw...)
    exp = Experiment(model, forcing; spinup=CycledSpinup(2001:2001; cycles=1), verbose=false, outputs=[
        out("daily.csv"; variables=:budget, frequency=:daily, aggregate=:area_sum),
        out("monthly.csv"; variables=:budget, frequency=:monthly, aggregate=:area_sum),
        out("yearly.csv"; variables=[:budget; :swe; :albedo], frequency=:yearly, aggregate=:area_sum),
        out("monthly.nc"; variables=[:melt, :swe, :albedo], frequency=:monthly),
        out("snap.nc"; variables=[:swe, :density], frequency=Month(6), kind=:snapshot),
        out("regions.csv"; variables=[:melt, :swe], frequency=:yearly, aggregate=:area_sum,
            regions=(west=[true, false, true, false], east=[false, true, false, true]))])
    result = run!(exp)
    @test result.years == [2001, 2002]

    # Period totals of fluxes add up across frequencies.
    hd, daily = _csv(joinpath(dir, "daily.csv"))
    _, monthly = _csv(joinpath(dir, "monthly.csv"))
    hy, yearly = _csv(joinpath(dir, "yearly.csv"))
    @test length(daily) == 365 + 365 && length(monthly) == 24 && length(yearly) == 2
    @test hd[3] == "surface_smb (Gt)"
    for j in 1:6
        @test sum(r[j] for r in daily[1:365]) ≈ yearly[1][j] rtol=1e-12
        @test sum(r[j] for r in monthly[13:24]) ≈ yearly[2][j] rtol=1e-12
    end
    # SMB = precipitation − runoff + vapour exchange.
    for r in yearly
        @test r[1] ≈ r[2] - r[4] + r[6] rtol=1e-12
    end
    @test yearly[1][3] > 0                       # some melt in the warm columns

    # Regions add up to the total.
    _, reg = _csv(joinpath(dir, "regions.csv"))
    @test reg[2][1] + reg[2][2] ≈ yearly[2][3] rtol=1e-12

    # Gridded output: monthly melt sums to the area-summed yearly melt.
    ds = NCDataset(joinpath(dir, "monthly.nc"))
    @test size(ds["melt"]) == (2, 2, 24)
    melt = Float64.(ds["melt"][:, :, 13:24])
    area = [1.0 3.0; 2.0 4.0]                     # (x, y) of columns 1–4
    @test sum(melt .* area) * 1e-6 ≈ yearly[2][3] rtol=1e-6
    @test ds["time"][1] == DateTime(2001, 1, 1, 12)
    close(ds)
    ds = NCDataset(joinpath(dir, "snap.nc"))
    @test size(ds["density"]) == (2, 2, 8, 4)
    close(ds)

    # The final state equals a run without any output.
    plain = run!(Experiment(model, forcing; spinup=CycledSpinup(2001:2001; cycles=1), verbose=false))
    @test plain.state.mass == result.state.mass
    @test plain.state.runoff == result.state.runoff
end

@testset "Experiment configuration (TOML)" begin
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
        [model]
        Ntot = 8
        [spinup]
        type = "none"
        [[output]]
        file = "budget.csv"
        variables = "budget"
        frequency = "1 month"
        aggregate = "area_sum"
        """)
    result = run_config(joinpath(dir, "run.toml"))
    _, rows = _csv(joinpath(dir, "out", "budget.csv"))
    @test length(rows) == 12
    @test all(r -> r[2] > 0, rows)               # precipitation in every month
    @test isfile(joinpath(dir, "out", "run.toml"))  # the configuration is kept with the results
end

@testset "relabel_year" begin
    f, grid = _experiment_forcing(2004:2004)            # a leap year
    g = relabel_year(f, 2005)
    @test year(first(g.time_values)) == 2005 && year(last(g.time_values)) == 2005
    @test count(t -> month(t) == 2 && day(t) == 28, g.time_values) == 2   # 29 February → 28 February
    @test g.day_of_year == f.day_of_year && g.dt_days == f.dt_days
    model = BESSIModel(grid; Ntot=8)
    run(forcing) = run!(Experiment(model, Forcing(forcing; grid); verbose=false)).state
    a, b = run(f), run(g)
    @test Array(a.mass) == Array(b.mass) && Array(a.runoff) == Array(b.runoff)
end

@testset "PDD and ITM budgets" begin
    f, grid = _experiment_forcing(2001:2001)
    for model in (PDDModel(grid), ITMModel(grid))
        forcing = model isa ITMModel ? SnowpackForcing(; time_values=f.time_values, dt_days=f.dt_days,
            air_temperature=f.air_temperature, snowfall_rate=f.snowfall_rate, rainfall_rate=f.rainfall_rate,
            shortwave_down=f.shortwave_down, latitude_deg=f.latitude_deg, surface_height=fill(1500.0, 4),
            ice_thickness=fill(1000.0, 4), annual_pdd=fill(50.0, 4)) : f
        dir = mktempdir()
        variables = [:surface_smb, :precipitation, :runoff]
        run!(Experiment(model, Forcing(forcing; grid, area=fill(100.0, 4)); verbose=false,
            outputs=Output(joinpath(dir, "budget.csv"); variables, frequency=:yearly, aggregate=:area_sum)))
        _, rows = _csv(joinpath(dir, "budget.csv"))
        smb, precip, runoff = rows[1]
        @test precip > 0 && runoff >= 0
        @test smb ≈ precip - runoff                      # no vapour exchange in PDD and ITM
    end
end

@testset "configuration overrides" begin
    @test Chion._set_value("0.84") == 0.84 && Chion._set_value("prescribed") == "prescribed"
    @test Chion._set_value("[1940, 1969]") == [1940, 1969]
    cfg = Chion._apply_overrides!(Dict{String, Any}("model" => Dict{String, Any}("type" => "bessi")),
        ["model.alpha_dry" => "0.84", "spinup.type" => "none"])
    @test cfg["model"]["alpha_dry"] == 0.84 && cfg["spinup"]["type"] == "none"
end
