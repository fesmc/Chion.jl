#!/usr/bin/env julia

# Shared-climatology-spinup comparison of native and reconstructed 3-hourly MAR forcing.
const COMPARISON_SCRIPT_START_NS = time_ns()
import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Chion
using Dates
using NCDatasets
using Statistics

# Include JIT compilation and setup in phase timings, outside Simulation timers.
function timed_phase(f, label)
    Base.@nospecialize f
    println("PHASE_START ", label); flush(stdout)
    measured = @timed Base.invokelatest(f)
    println("PHASE_END ", label, " wall_s=", measured.time,
        " compile_s=", measured.compile_time, " recompile_s=", measured.recompile_time,
        " gc_s=", measured.gctime, " allocated_bytes=", measured.bytes)
    flush(stdout)
    return measured.value
end
println("SCRIPT_IMPORTS wall_s=", (time_ns() - COMPARISON_SCRIPT_START_NS) / 1e9); flush(stdout)

const CFG = (
    climatology=get(ENV, "CHION_CLIMATOLOGY_FILE", "/p/projects/ou/labs/ai/Nils/MAR3.14/MARv3.14.3-10km-daily-ERA5-1940-1980_daily_climatology.nc"),
    daily=get(ENV, "CHION_DAILY_FILE", "/p/projects/ou/labs/ai/Nils/MAR_daily/mar_daily/MARv3.14.3-10km-daily-ERA5-1981.nc"),
    native=get(ENV, "CHION_3HOURLY_FILE", "/p/projects/ou/labs/ai/Nils/MAR_daily/mar_3hourly/MARv3.14.3-1981.nc"),
    out=get(ENV, "CHION_OUTPUT_DIR", joinpath(@__DIR__, "..", "..", "output", "gris_1981_3hourly_comparison_v2")),
    spinup_years=parse(Int, get(ENV, "CHION_SPINUP_YEARS", "200")),
    backend=Symbol(get(ENV, "CHION_BACKEND", "gpu")), mask=50.0,
)

columns(A, grid, ntime, name, dims) = begin
    B = Chion._column_matrix(Chion._as_time_y_x(Float64.(A), dims, ntime, length(grid.y), length(grid.x), name))
    B[grid.is .+ (grid.js .- 1) .* length(grid.x), :]
end

function field(path, name, grid, ntime)
    NCDataset(path) do ds
        columns(ds[name].var[ntuple(_ -> (:), ndims(ds[name].var))...], grid, ntime, name, dimnames(ds[name].var))
    end
end

function static_field(path, name, grid)
    NCDataset(path) do ds
        A = Chion._as_y_x(Float64.(ds[name].var[:,:]), dimnames(ds[name].var), length(grid.y), length(grid.x), name)
        [A[grid.js[k], grid.is[k]] for k in eachindex(grid.js)]
    end
end

function load_daily(path)
    loaded = load_forcing_file(path; time_name="TIME", air_temperature_name="TTZ",
        wind_speed_name="UVZ", air_pressure_name="SP", mask_name="MSK", mask_threshold=CFG.mask)
    f, g = loaded.forcing, loaded.grid
    # MAR stores SP in hPa, whereas Chion uses Pa.
    p = f.air_pressure .* 100
    a = clamp.(field(path, "AL2", g, length(f.time_values)), 0, 1)
    swu = field(path, "SWU", g, length(f.time_values))
    SnowpackForcing(time_values=f.time_values, dt_days=f.dt_days, air_temperature=f.air_temperature,
        snowfall_rate=f.snowfall_rate, rainfall_rate=f.rainfall_rate, shortwave_down=f.shortwave_down,
        wind_speed=f.wind_speed, q_sw_net=f.shortwave_down .- swu, has_q_sw_net=true,
        q_lw_down=f.q_lw_down, has_q_lw_down=true, q_sh=f.q_sh, has_q_sh=true,
        q_lh=f.q_lh, has_q_lh=true, relative_humidity=f.relative_humidity,
        has_relative_humidity=f.has_relative_humidity, air_pressure=p, surface_height=f.surface_height,
        prescribed_albedo=a, has_prescribed_albedo=true, latitude_deg=f.latitude_deg), g
end

function load_native(path)
    loaded = load_forcing_file(path; x_name="X14_163", y_name="Y19_288", time_name="TIME",
        air_temperature_name="TT", wind_speed_name=nothing, relative_humidity_name=nothing,
        air_pressure_name="SP", prescribed_albedo_name="AL", mask_name="MSK", mask_threshold=CFG.mask)
    f, g = loaded.forcing, loaded.grid; n = length(f.time_values)
    u, v, q = field(path,"UU",g,n), field(path,"VV",g,n), field(path,"QQ",g,n) ./ 1000
    p = f.air_pressure .* 100
    # Specific humidity to vapour pressure and then RH over ice/water using Chion's closure.
    e = q .* p ./ (0.622 .+ 0.378 .* q)
    rh = similar(e)
    for j in axes(e,2), i in axes(e,1)
        es = Chion._bessi_ice_saturation_vapor_pressure(f.air_temperature[i,j], 273.15)
        rh[i,j] = clamp(100e[i,j] / es, 0, 100)
    end
    a = clamp.(f.prescribed_albedo, 0, 1)
    forcing = SnowpackForcing(time_values=f.time_values .- Minute(90), dt_days=f.dt_days,
        air_temperature=f.air_temperature, snowfall_rate=f.snowfall_rate .* 8,
        rainfall_rate=f.rainfall_rate .* 8, shortwave_down=f.shortwave_down,
        wind_speed=hypot.(u,v), q_sw_net=f.shortwave_down .* (1 .- a), has_q_sw_net=true,
        q_lw_down=f.q_lw_down, has_q_lw_down=true, q_sh=f.q_sh, has_q_sh=true,
        q_lh=f.q_lh, has_q_lh=true, relative_humidity=rh, has_relative_humidity=true,
        air_pressure=p, surface_height=f.surface_height, prescribed_albedo=a,
        has_prescribed_albedo=true, latitude_deg=f.latitude_deg)
    forcing, g
end

function expand_daily(f; prescribed=true, temperature_amplitude_c=1.0)
    ncol, nd = size(f.air_temperature); nt = 8nd
    mats = Dict{Symbol,Matrix{Float64}}(n => Matrix{Float64}(undef,ncol,nt) for n in
        (:air_temperature,:snowfall_rate,:rainfall_rate,:shortwave_down,:wind_speed,:q_sw_net,
         :q_lw_down,:q_sh,:q_lh,:relative_humidity,:air_pressure,:surface_height,:prescribed_albedo))
    times=DateTime[]
    for d in 1:nd, s in 1:8
        k=8(d-1)+s; h0=-pi+(s-1)*pi/4; h1=h0+pi/4
        push!(times, DateTime(Date(f.time_values[d])) + Hour(3s) - Minute(90))
        # Copy daily-constant fields once per interval, outside the cell loop.
        for n in (:snowfall_rate,:rainfall_rate,:wind_speed,:q_lw_down,:q_sh,:q_lh,:relative_humidity,:air_pressure,:surface_height,:prescribed_albedo)
            @views mats[n][:,k] .= getproperty(f,n)[:,d]
        end
        air_temperature=mats[:air_temperature]
        shortwave_down=mats[:shortwave_down]
        q_sw_net=mats[:q_sw_net]
        sol=f.solar_longitude_deg[d]
        for i in 1:ncol
            lat=f.latitude_deg[i,d]
            air_temperature[i,k]=Chion._diurnal_temperature_interval_average(
                f.air_temperature[i,d],temperature_amplitude_c,h0,h1)
            shortwave_down[i,k]=Chion._diurnal_shortwave_interval_average(f.shortwave_down[i,d],lat,sol,h0,h1)
            q_sw_net[i,k]=Chion._diurnal_shortwave_interval_average(f.q_sw_net[i,d],lat,sol,h0,h1)
        end
    end
    SnowpackForcing(; time_values=times,dt_days=fill(0.125,nt),mats...,
        has_q_sw_net=prescribed,has_q_lw_down=prescribed,has_q_sh=prescribed,has_q_lh=prescribed,
        has_relative_humidity=true,has_prescribed_albedo=prescribed,latitude_deg=f.latitude_deg[:,1])
end

"""Chion for 3-hourly forcing: the package defaults (calibrated setup) without the
internal diurnal reconstruction, because native and expanded forcing are already sub-daily."""
model(grid; kwargs...) = BESSIModel(grid; Ntot=15, diurnal_shortwave_substeps=false,
    diurnal_temperature_cycle=false, kwargs...)

budget(s)=(melt=Array(s.melt),runoff=Array(s.runoff),refreezing=Array(s.refreezing),
    vapor_mass=Array(s.vapor_mass),smb_ice=Array(s.smb_ice))

function run_branch(m, forcing, initial, name, area)
    before=budget(initial)
    sim=Simulation(m;forcing,state=deepcopy(initial),years=1,backend=CFG.backend,
        write_netcdf=true,netcdf_variables=:monthly,netcdf_path=joinpath(CFG.out,"$(name)_monthly.nc"),
        compute_year_metrics=false,name=name)
    result=run!(sim); result.status==:complete || error("$name failed")
    after=budget(sim.now)
    vals=NamedTuple{keys(before)}(Tuple(sum((getfield(after,k).-getfield(before,k)).*area)*1e-6 for k in keys(before)))
    return vals
end

function main()
    mkpath(CFG.out)
    open(joinpath(CFG.out,"RUNNING.txt"),"w") do io
        println(io,"Shared 200-year climatology spin-up followed by four 1981 branches.")
        println(io,"Started: ",now())
    end
    clim,g=timed_phase("load_climatology") do; load_daily(CFG.climatology); end
    daily,gd=timed_phase("load_daily") do; load_daily(CFG.daily); end
    native,gn=timed_phase("load_native") do; load_native(CFG.native); end
    (g.x == gd.x && g.x == gn.x && g.y == gd.y && g.y == gn.y &&
     g.is == gd.is && g.is == gn.is && g.js == gd.js && g.js == gn.js) || error("MAR grids differ")
    m=model(g); area=static_field(CFG.daily,"AREA",g)
    expanded_clim=timed_phase("expand_climatology") do; expand_daily(clim); end
    spin=timed_phase("construct_spinup") do
        Simulation(m;forcing=expanded_clim,years=CFG.spinup_years,backend=CFG.backend,
            write_netcdf=false,compute_year_metrics=false,name="shared_prescribed_climatology_spinup")
    end
    timed_phase("run_spinup_including_JIT") do
        run!(spin).status==:complete || error("spinup failed")
    end
    expanded_daily=timed_phase("expand_daily_prescribed") do; expand_daily(daily); end
    expanded_parameterized=timed_phase("expand_daily_parameterized") do; expand_daily(daily;prescribed=false); end
    cases=(("native_3hourly_prescribed",native),("daily_prescribed",daily),
           ("daily_prescribed_8x",expanded_daily),("daily_parameterized_8x",expanded_parameterized))
    open(joinpath(CFG.out,"annual_budgets_gt.csv"),"w") do io
        println(io,"case,melt,runoff,refreezing,vapor_mass,smb_ice")
        for (name,f) in cases
            b=timed_phase("branch_"*name) do
                run_branch(m,f,spin.now,name,area)
            end
            println(io,join((name,values(b)...),','))
        end
    end
    open(joinpath(CFG.out,"configuration.csv"),"w") do io
        println(io,"setting,value"); for (k,v) in pairs(CFG); println(io,"$k,$v"); end
        println(io,"seb_scheme,semix")
        println(io,"layers,15\nnear_surface_layer_max_m,0.02|0.05|0.10|0.30\nshared_spinup_substeps,8")
        println(io,"native_precipitation_conversion,mmWE_per_3h_divided_by_10800s")
    end
    rm(joinpath(CFG.out,"RUNNING.txt"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    timed_phase(main, "main_including_JIT")
end
