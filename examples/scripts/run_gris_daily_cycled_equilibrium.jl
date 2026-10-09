#!/usr/bin/env julia

# Daily-forcing "equilibrium" without the climatology: the 1940-1980 daily
# climatology has snowfall on almost every day, which keeps a dynamic albedo
# artificially fresh. Instead the real daily MAR years
# CHION_CYCLE_FIRST_YEAR..CHION_CYCLE_LAST_YEAR are repeated CHION_CYCLES times
# from the initial state. The final cycle is compared with MAR's mean over the
# same years; earlier cycles document the drift towards equilibrium. The forcing
# mode and surface setup follow run_gris_fully_prescribed_monthly_transient.jl
# (CHION_FORCING_MODE and optional model overrides such as CHION_ALPHA_*).
include(joinpath(@__DIR__, "run_gris_fully_prescribed_monthly_transient.jl"))

const CYCLE_YEARS = parse(Int, get(ENV, "CHION_CYCLE_FIRST_YEAR", "1940")):parse(Int, get(ENV, "CHION_CYCLE_LAST_YEAR", "1979"))
const CYCLES = parse(Int, get(ENV, "CHION_CYCLES", "5"))
const CYCLE_VARS = (:melt, :runoff, :refreezing, :surface_smb, :vapor_mass)
const CYCLE_OUT = MONTHLY_CONFIG.output_dir

cycle_year_path(year) = joinpath(MONTHLY_CONFIG.forcing_dir, "MARv3.14.3-$(DOMAIN_SETUP.file_tag)-daily-ERA5-$(year).nc")

function write_cycle_plots(rows)
    cycles = 1:CYCLES
    panels = Any[]
    for (name, label) in ((:melt, "Melt"), (:runoff, "Runoff"), (:refreezing, "Refreezing"), (:surface_smb, "Surface SMB"))
        bias = [mean(getproperty(r.chion, name) - getproperty(r.mar, name) for r in rows if r.cycle == c) for c in cycles]
        p = plot(cycles, bias; c="#2a78d6", lw=2, marker=:circle, ms=6, msw=0, label="", xticks=cycles,
            xlabel="Cycle of $(first(CYCLE_YEARS))–$(last(CYCLE_YEARS))", ylabel="Chion − MAR [Gt yr⁻¹]",
            title=label, framestyle=:box, grid=:y, gridalpha=0.15)
        plot!(p, [first(cycles), last(cycles)], [0, 0]; c=:gray, lw=1, label="", primary=false)
        push!(panels, p)
    end
    savefig(plot(panels...; layout=(2, 2), size=(1200, 800), left_margin=5Plots.mm,
            plot_title="Daily cycled equilibrium ($(FORCING_MODE)): drift of the cycle-mean bias", plot_titlefontsize=12),
        joinpath(CYCLE_OUT, "cycle_drift.pdf"))

    final = filter(r -> r.cycle == CYCLES, rows)
    years = getproperty.(final, :year)
    panels = Any[]
    for (name, label) in ((:melt, "Melt"), (:runoff, "Runoff"), (:refreezing, "Refreezing"), (:surface_smb, "Surface SMB"))
        p = plot(years, getproperty.(getproperty.(final, :mar), name); c=:black, lw=2.5, label="MAR", xlabel="Year",
            ylabel="Gt yr⁻¹", title=label, framestyle=:box, grid=:y, gridalpha=0.15, legend=name == :melt ? :topleft : false)
        plot!(p, years, getproperty.(getproperty.(final, :chion), name); c="#eb6834", lw=2, label="Chion (final cycle)")
        push!(panels, p)
    end
    savefig(plot(panels...; layout=(2, 2), size=(1200, 800), left_margin=5Plots.mm), joinpath(CYCLE_OUT, "final_cycle_budgets.pdf"))
end

function daily_cycled_main()
    mkpath(CYCLE_OUT)
    first_loaded = load_mar(cycle_year_path(first(CYCLE_YEARS)))
    model = build_mode_model(first_loaded.grid)
    area_km2 = read_mar_static_columns(cycle_year_path(first(CYCLE_YEARS)), "AREA", model.grid)
    state = initial_state(model)
    mar = Dict{Int,Any}()
    rows = NamedTuple[]
    open(joinpath(CYCLE_OUT, "cycled_budgets.csv"), "w") do io
        println(io, "cycle,year,", join(("mar_$(v)_gt,chion_$(v)_gt" for v in CYCLE_VARS), ','))
        for cycle in 1:CYCLES, year in CYCLE_YEARS
            path = cycle_year_path(year)
            loaded = load_mar(path)
            same_columns(model.grid, loaded.grid) || error("Forcing grid differs in $(basename(path)).")
            forcing = mode_forcing(path, loaded)
            before = state_budget(state)
            sim = Simulation(model; forcing, state, years=1, backend=CONFIG.backend, write_netcdf=false,
                compute_year_metrics=false, name="daily_cycle$(cycle)_$(year)")
            run!(sim).status == :complete || error("Cycle $cycle year $year failed")
            state = sim.now
            reference = get!(() -> mar_annual_budget(path, model.grid, forcing, area_km2), mar, year)
            chion = chion_annual_budget(before, state_budget(state), forcing, area_km2)
            push!(rows, (; cycle, year, mar=reference, chion))
            println(io, cycle, ',', year, ',', join((string(getproperty(reference, v), ',', getproperty(chion, v)) for v in CYCLE_VARS), ','))
            flush(io)
            @info "Completed daily cycle year" cycle year
            GC.gc()
        end
    end

    final = filter(r -> r.cycle == CYCLES, rows)
    open(joinpath(CYCLE_OUT, "equilibrium_annual_integrated_budgets.csv"), "w") do io
        println(io, "variable,MAR_Gt,Chion_Gt,Chion_minus_MAR_Gt")
        for v in CYCLE_VARS
            mar_mean = mean(getproperty(r.mar, v) for r in final)
            chion_mean = mean(getproperty(r.chion, v) for r in final)
            println(io, join((v, mar_mean, chion_mean, chion_mean - mar_mean), ','))
        end
    end
    write_cycle_plots(rows)
    open(joinpath(CYCLE_OUT, "configuration.csv"), "w") do io
        println(io, "setting,value\nforcing,$(FORCING_MODE)_daily_MAR_cycled\ncycle_years,$(first(CYCLE_YEARS))-$(last(CYCLE_YEARS))\ncycles,$CYCLES")
        println(io, "spinup,none (cycles start from the initial state)")
        write_model_configuration(io, model)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    daily_cycled_main()
end
