#!/usr/bin/env julia

# Native 3-hourly "equilibrium": no 3-hourly climatology exists, so after the
# shared 200-y daily-climatology spin-up (8 reconstructed 3-hour steps) the
# native MAR years CHION_NATIVE_CYCLE_FIRST_YEAR..LAST_YEAR are repeated
# CHION_NATIVE_CYCLES times. The final cycle is compared with MAR's mean over
# the same years; earlier cycles document the drift towards equilibrium.
include(joinpath(@__DIR__, "run_gris_native_3hourly_transient.jl"))

const CYCLE_YEARS = parse(Int, get(ENV, "CHION_NATIVE_CYCLE_FIRST_YEAR", "1980")):parse(Int, get(ENV, "CHION_NATIVE_CYCLE_LAST_YEAR", "1989"))
const CYCLES = parse(Int, get(ENV, "CHION_NATIVE_CYCLES", "6"))
const CYCLE_VARS = (:melt, :runoff, :refreezing, :surface_smb, :vapor_mass)

native_year_path(year) = joinpath(TRANSIENT_NATIVE_DIR, "MARv3.14.3-$(year).nc")
daily_year_path(year) = joinpath(TRANSIENT_DAILY_DIR, "MARv3.14.3-10km-daily-ERA5-$(year).nc")

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
            plot_title="Native 3-hourly cycled equilibrium: drift of the cycle-mean bias", plot_titlefontsize=12),
        joinpath(CFG.out, "cycle_drift.pdf"))

    final = filter(r -> r.cycle == CYCLES, rows)
    years = getproperty.(final, :year)
    panels = Any[]
    for (name, label) in ((:melt, "Melt"), (:runoff, "Runoff"), (:refreezing, "Refreezing"), (:surface_smb, "Surface SMB"))
        p = plot(years, getproperty.(getproperty.(final, :mar), name); c=:black, lw=2.5, label="MAR", xlabel="Year",
            ylabel="Gt yr⁻¹", title=label, framestyle=:box, grid=:y, gridalpha=0.15, legend=name == :melt ? :topleft : false)
        plot!(p, years, getproperty.(getproperty.(final, :chion), name); c="#eb6834", lw=2, label="Chion, native 3-hourly (final cycle)")
        push!(panels, p)
    end
    savefig(plot(panels...; layout=(2, 2), size=(1200, 800), left_margin=5Plots.mm), joinpath(CFG.out, "final_cycle_budgets.pdf"))
end

function cycled_main()
    mkpath(CFG.out)
    clim, grid = load_daily(CFG.climatology)
    m = model(grid)
    area = static_field(daily_year_path(first(CYCLE_YEARS)), "AREA", grid)
    spin = Simulation(m; forcing=expand_daily(clim), years=CFG.spinup_years, backend=CFG.backend,
        write_netcdf=false, compute_year_metrics=false, name="native_cycled_equilibrium_spinup")
    run!(spin).status == :complete || error("200-year climatology spin-up failed")
    state = spin.now

    mar = Dict{Int,Any}()
    rows = NamedTuple[]
    monthly_dir = joinpath(CFG.out, "final_cycle_monthly")
    open(joinpath(CFG.out, "cycled_budgets.csv"), "w") do io
        println(io, "cycle,year,", join(("mar_$(v)_gt,chion_$(v)_gt" for v in CYCLE_VARS), ','))
        for cycle in 1:CYCLES, year in CYCLE_YEARS
            forcing, this_grid = load_native(native_year_path(year))
            same_native_grid(grid, this_grid) || error("Grid mismatch in $year")
            final = cycle == CYCLES
            final && mkpath(monthly_dir)
            before = budget(state)
            sim = Simulation(m; forcing, state=deepcopy(state), years=1, backend=CFG.backend,
                write_netcdf=final, netcdf_variables=:monthly,
                netcdf_path=final ? joinpath(monthly_dir, "native_$(year)_monthly.nc") : nothing,
                compute_year_metrics=false, name="native_cycle$(cycle)_$(year)")
            run!(sim).status == :complete || error("Cycle $cycle year $year failed")
            state = sim.now
            reference = get!(mar, year) do
                merge(mar_native_budget(native_year_path(year), grid, area, length(forcing.dt_days)),
                    mar_daily_hydro_budget(daily_year_path(year), grid, area))
            end
            chion = chion_budget_delta(before, budget(state), forcing, area)
            push!(rows, (; cycle, year, mar=reference, chion))
            println(io, cycle, ',', year, ',', join((string(getproperty(reference, v), ',', getproperty(chion, v)) for v in CYCLE_VARS), ','))
            flush(io)
            @info "Completed native cycle year" cycle year
            GC.gc()
        end
    end

    final = filter(r -> r.cycle == CYCLES, rows)
    open(joinpath(CFG.out, "equilibrium_annual_integrated_budgets.csv"), "w") do io
        println(io, "variable,MAR_Gt,Chion_Gt,Chion_minus_MAR_Gt")
        for v in CYCLE_VARS
            mar_mean = mean(getproperty(r.mar, v) for r in final)
            chion_mean = mean(getproperty(r.chion, v) for r in final)
            println(io, join((v, mar_mean, chion_mean, chion_mean - mar_mean), ','))
        end
    end
    write_cycle_plots(rows)
    open(joinpath(CFG.out, "configuration.csv"), "w") do io
        println(io, "setting,value\nforcing,native_3hourly_MAR_cycled\ncycle_years,$(first(CYCLE_YEARS))-$(last(CYCLE_YEARS))\ncycles,$CYCLES")
        println(io, "spinup,$(CFG.spinup_years)y daily climatology, 8 reconstructed 3-hour steps\nseb_scheme,semix")
        println(io, "near_surface_layer_max_thicknesses_m,$(join(m.parameters.near_surface_layer_max_thicknesses_m, '|'))")
        println(io, "ice_substrate_layers,$(m.parameters.ice_substrate_layers)\nice_substrate_top_thickness_m,$(m.parameters.ice_substrate_top_thickness_m)")
        println(io, "refreezing_reference,daily_MAR_RZ\nvapor_reference,daily_MAR_SU")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    cycled_main()
end
