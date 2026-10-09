#!/usr/bin/env julia

# Run a continuous native-3-hourly MAR transient after shared 200-y spin-up.
include(joinpath(@__DIR__, "run_gris_1981_three_hourly_comparison.jl"))
using Printf
using Plots

const TRANSIENT_NATIVE_DIR = get(ENV, "CHION_3HOURLY_DIR", "/p/projects/ou/labs/ai/Nils/MAR_daily/mar_3hourly")
const TRANSIENT_DAILY_DIR = get(ENV, "CHION_DAILY_DIR", "/p/projects/ou/labs/ai/Nils/MAR_daily/mar_daily")
same_native_grid(a, b) = a.x == b.x && a.y == b.y && a.js == b.js && a.is == b.is

function gt_sum(field_values, area)
    return sum(field_values .* area) * 1e-6
end

function mar_native_budget(path, grid, area, ntime)
    total(name) = gt_sum(vec(sum(field(path, name, grid, ntime); dims=2)), area)
    return (surface_smb=total("SMB"), melt=total("ME"), runoff=total("RU"))
end

function mar_daily_hydro_budget(path, grid, area)
    # Reference budgets need only the calendar, RZ and SU, not every forcing field.
    times = NCDataset(path) do ds
        Chion._read_time_values(ds, "TIME", length(ds["TIME"]))
    end
    dt_days = Chion.infer_dt_days(times)
    n = length(dt_days)
    annual(name) = vec(sum(field(path, name, grid, n) .* reshape(dt_days, 1, :); dims=2))
    return (refreezing=gt_sum(annual("RZ"), area), vapor_mass=gt_sum(-annual("SU"), area))
end

function chion_budget_delta(before, after, forcing, area)
    change(name) = getfield(after, name) .- getfield(before, name)
    precipitation = vec(sum((forcing.snowfall_rate .+ forcing.rainfall_rate) .*
        reshape(forcing.dt_days .* 86_400, 1, :); dims=2))
    runoff = change(:runoff)
    vapor = change(:vapor_mass)
    return (surface_smb=gt_sum(precipitation .- runoff .+ vapor, area),
        melt=gt_sum(change(:melt), area), runoff=gt_sum(runoff, area),
        refreezing=gt_sum(change(:refreezing), area), vapor_mass=gt_sum(vapor, area))
end

function write_native_transient_plot(path, rows)
    years = getproperty.(rows, :year)
    panels = Any[]
    for (name, label) in ((:melt, "Melt [Gt/yr]"), (:runoff, "Runoff [Gt/yr]"),
                          (:refreezing, "Refreezing [Gt/yr]"), (:surface_smb, "Surface SMB [Gt/yr]"))
        mar = getproperty.(getproperty.(rows, :mar), name)
        chion = getproperty.(getproperty.(rows, :chion), name)
        push!(panels, plot(years, mar; label="MAR", lw=2.5, c=:black, xlabel="Year", ylabel=label,
            title=label, legend=:best, framestyle=:box))
        plot!(panels[end], years, chion; label="Chion, native 3-hour", lw=2, c=:royalblue)
    end
    savefig(plot(panels...; layout=(2,2), size=(1200,800)), path)
end

function native_transient_main()
    mkpath(CFG.out)
    files = filter(p -> occursin(r"^MARv3\.14\.3-(\d{4})\.nc$", basename(p)), readdir(TRANSIENT_NATIVE_DIR; join=true))
    sort!(files; by=p -> parse(Int, only(match(r"-(\d{4})\.nc$", basename(p)).captures)))
    years = [parse(Int, only(match(r"-(\d{4})\.nc$", basename(p)).captures)) for p in files]
    first(years) == 1980 && last(years) == 2017 || error("Expected complete native years 1980--2017, found $(first(years))--$(last(years))")

    clim, grid = load_daily(CFG.climatology)
    m = model(grid)
    area = static_field(joinpath(TRANSIENT_DAILY_DIR, "MARv3.14.3-10km-daily-ERA5-1980.nc"), "AREA", grid)
    spin = Simulation(m; forcing=expand_daily(clim), years=CFG.spinup_years, backend=CFG.backend,
        write_netcdf=false, compute_year_metrics=false, name="shared_200y_climatology_spinup_native_3h")
    run!(spin).status == :complete || error("200-year climatology spin-up failed")
    state = spin.now

    csv = joinpath(CFG.out, "native_3hourly_transient_1980_2017_budgets.csv")
    rows = NamedTuple[]
    open(csv, "w") do io
        println(io, "year,mar_melt_gt,chion_melt_gt,mar_runoff_gt,chion_runoff_gt,mar_refreezing_gt,chion_refreezing_gt,mar_surface_smb_gt,chion_surface_smb_gt,mar_vapor_mass_gt,chion_vapor_mass_gt")
        for (path, year) in zip(files, years)
            forcing, this_grid = load_native(path)
            same_native_grid(grid, this_grid) || error("Grid mismatch in $(basename(path))")
            dpath = joinpath(TRANSIENT_DAILY_DIR, "MARv3.14.3-10km-daily-ERA5-$(year).nc")
            isfile(dpath) || error("Missing daily reference $dpath")
            ntime = length(forcing.dt_days)
            before = budget(state)
            monthly_dir = get(ENV, "CHION_NATIVE_MONTHLY_DIR", "")
            isempty(monthly_dir) || mkpath(monthly_dir)
            sim = Simulation(m; forcing, state=deepcopy(state), years=1, backend=CFG.backend,
                write_netcdf=!isempty(monthly_dir), netcdf_variables=:monthly,
                netcdf_path=isempty(monthly_dir) ? nothing : joinpath(monthly_dir, "native_$(year)_monthly.nc"),
                compute_year_metrics=false, name="native_3hourly_transient_$(year)")
            result = run!(sim)
            result.status == :complete || error("Native transient year $year failed: $(result.status)")
            state = sim.now
            after = budget(state)
            mar = merge(mar_native_budget(path, grid, area, ntime), mar_daily_hydro_budget(dpath, grid, area))
            chion = chion_budget_delta(before, after, forcing, area)
            push!(rows, (;year, mar, chion))
            println(io, join((year, mar.melt, chion.melt, mar.runoff, chion.runoff,
                mar.refreezing, chion.refreezing, mar.surface_smb, chion.surface_smb,
                mar.vapor_mass, chion.vapor_mass), ','))
            flush(io)
            @info "Completed continuous native 3-hour transient year" year
            GC.gc()
        end
    end
    write_native_transient_plot(joinpath(CFG.out, "native_3hourly_transient_1980_2017_budgets.pdf"), rows)
    open(joinpath(CFG.out, "native_3hourly_transient_configuration.csv"), "w") do io
        println(io, "setting,value\nforcing,native_3hourly_MAR_1980_2017\nspinup_years,$(CFG.spinup_years)\nbackend,$(CFG.backend)\nseb_scheme,semix\nlayers,15\nnear_surface_profile_m,0.02|0.05|0.10|0.30\nice_substrate_layers,$(m.parameters.ice_substrate_layers)\nice_substrate_top_thickness_m,$(m.parameters.ice_substrate_top_thickness_m)\ntransient_steps,actual_3_hourly\nrefreezing_reference,daily_MAR_RZ\nvapor_reference,daily_MAR_SU")
    end
    @info "Wrote native 3-hourly transient budgets and plot" csv
end

if abspath(PROGRAM_FILE) == @__FILE__
    native_transient_main()
end
