#!/usr/bin/env julia

"""Write monthly output for a fully prescribed 1940--2026 GrIS transient.

The model is first spun up with climatological MAR forcing, then continued one
annual MAR file at a time. This avoids loading the 47 GB consolidated forcing
file, which exceeds a GPU job's host-memory allocation. Monthly records are
staged per year and consolidated after the complete transient.
"""

include(joinpath(@__DIR__, "run_gris_fully_prescribed_transient.jl"))

using Dates: DateTime, Millisecond, month

const MONTHLY_CONFIG = (
    forcing_dir=get(ENV, "CHION_FORCING_DIR", DOMAIN_SETUP.forcing_dir),
    output_dir=get(ENV, "CHION_OUTPUT_DIR", joinpath(@__DIR__, "..", "..", "output", "gris_prescribed_transient_1940_2026")),
)

# CHION_FORCING_MODE selects prescribed MAR surface fluxes (default), the
# fully parameterized SEB with dynamic albedo, or the parameterized SEB with
# MAR albedo (AL2) prescribed (`parameterized_mar_albedo`). Output columns keep
# the `chion_prescribed_*` names in every mode; configuration.csv records it.
const FORCING_MODE = Symbol(get(ENV, "CHION_FORCING_MODE", "prescribed"))
FORCING_MODE in (:prescribed, :parameterized, :parameterized_mar_albedo) ||
    error("CHION_FORCING_MODE must be prescribed, parameterized or parameterized_mar_albedo.")
FORCING_MODE == :parameterized_mar_albedo && (ENV["CHION_ALBEDO_SCHEME"] = "prescribed")
build_mode_model(grid) = FORCING_MODE == :prescribed ? build_fully_prescribed_model(grid) : build_parameterized_model(grid)
mode_forcing(path, loaded) =
    FORCING_MODE == :prescribed ? fully_prescribed_forcing(path, loaded) :
    FORCING_MODE == :parameterized ? parameterized_forcing(loaded.forcing; wind_speed=forcing_wind(path, loaded)) :
    parameterized_forcing(fully_prescribed_forcing(path, loaded); use_prescribed_albedo=true, wind_speed=forcing_wind(path, loaded))

# The last forcing year (2026) is incomplete, so it stays in the CSV but is
# excluded from plotted series and error metrics.
complete_year_rows(rows) = length(rows) > 1 ? rows[1:end-1] : rows

function staging_paths(staging_dir, files)
    return [joinpath(staging_dir, replace(basename(path), r"\.nc$" => "_monthly.nc")) for path in files]
end

"""Consolidate staged monthly outputs without retaining all forcing in memory."""
function consolidate_monthly_outputs!(model, state, staged_files, output_file)
    options = RunOptions(
        name="gris_prescribed_monthly_1940_2026",
        netcdf_path=output_file,
        write_netcdf=true,
        netcdf_variables=:monthly,
        years=1,
    )
    output = Chion.init_state_netcdf(
        output_file,
        options,
        DateTime[],
        model.grid,
        state,
        Chion.monthly_output_variables(model);
        nlayer=state.Ntot,
    )
    first_record = 1
    try
        for staged_file in staged_files
            NCDataset(staged_file) do staged
                nrecords = size(staged["smb_ice"], 3)
                last_record = first_record + nrecords - 1
                output.dataset["t"][first_record:last_record] = staged["t"][1:nrecords]
                for key in keys(output.vars)
                    output.vars[key][:, :, first_record:last_record] = staged[String(key)][:, :, 1:nrecords]
                end
                first_record = last_record + 1
            end
        end
        Chion.close_output!(output, :complete, first_record - 1)
    catch
        Chion.close_output!(output, :failed, first_record - 1)
        rethrow()
    end
    return output_file
end

"""Write the ice-sheet-integrated monthly surface SMB accompanying the maps."""
function write_surface_smb_timeseries(path, monthly_file, grid, area_km2)
    NCDataset(monthly_file) do dataset
        haskey(dataset, "surface_smb") || return nothing
        surface_smb = dataset["surface_smb"]
        times = dataset["t"][:]
        open(path, "w") do io
            println(io, "date,year,month,surface_smb_gt_per_month")
            for record in eachindex(times)
                field = surface_smb[:, :, record]
                values = Float64[
                    coalesce(field[grid.is[column], grid.js[column]], 0.0)
                    for column in eachindex(grid.is)
                ]
                timestamp = times[record] isa DateTime ?
                    times[record] :
                    DateTime(1970, 1, 1) + Millisecond(round(Int, times[record] * 86_400_000))
                integrated_gt = sum(values .* area_km2) * 1e-6
                println(io, "$(timestamp),$(year(timestamp)),$(month(timestamp)),$(integrated_gt)")
            end
        end
    end
    return path
end

function write_streamed_budget_comparison(path, rows)
    open(path, "w") do io
        columns = String["year"]
        for name in keys(PLOT_VARIABLES)
            append!(columns, ("mar_$(name)_gt_per_year", "chion_prescribed_$(name)_gt_per_year"))
        end
        println(io, join(columns, ','))
        for row in rows
            values = Any[row.year]
            for name in keys(PLOT_VARIABLES)
                append!(values, (getfield(row.mar, name), getfield(row.chion_prescribed, name)))
            end
            println(io, join(values, ','))
        end
    end
    return path
end

function monthly_transient_main()
    CONFIG.spinup_years > 0 || error("CHION_SPINUP_YEARS must be positive.")
    isdir(MONTHLY_CONFIG.forcing_dir) || error("Forcing directory not found: $(MONTHLY_CONFIG.forcing_dir)")
    mkpath(MONTHLY_CONFIG.output_dir)
    files = annual_files(MONTHLY_CONFIG.forcing_dir)
    years = [parse(Int, only(match(r"-(\d{4})\.nc$", basename(path)).captures)) for path in files]
    first(years) == DOMAIN_SETUP.first_year && last(years) == 2026 ||
        error("Expected annual forcing files for $(DOMAIN_SETUP.first_year)--2026, got $(first(years))--$(last(years)).")
    if DOMAIN == :greenland
        # The raw Greenland 2026 file lacks fields; a prepared compatible copy replaces it.
        compatible_2026 = get(ENV, "CHION_2026_COMPATIBLE_FILE", joinpath(
            dirname(MONTHLY_CONFIG.output_dir),
            "mar_forcing",
            "MARv3.14.3-10km-daily-ERA5-2026_chion_compatible.nc",
        ))
        isfile(compatible_2026) || error("2026 compatible forcing file not found: $compatible_2026")
        files[end] = compatible_2026
    end

    output_file = get(ENV, "CHION_MONTHLY_OUTPUT_FILE", joinpath(MONTHLY_CONFIG.output_dir, "gris_prescribed_monthly_1940_2026.nc"))
    timeseries_file = get(ENV, "CHION_SURFACE_SMB_TIMESERIES_FILE", joinpath(MONTHLY_CONFIG.output_dir, "monthly_surface_smb_timeseries.csv"))
    staging_dir = joinpath(MONTHLY_CONFIG.output_dir, ".gris_prescribed_monthly_staging")

    spinup_loaded = load_mar(CONFIG.spinup_file)
    model = build_mode_model(spinup_loaded.grid)
    staged_files = staging_paths(staging_dir, files)
    if isdir(staging_dir)
        all(isfile, staged_files) || error("Incomplete staging directory: $staging_dir")
        isfile(output_file) ||
            consolidate_monthly_outputs!(model, initial_state(model), staged_files, output_file)
        write_surface_smb_timeseries(
            timeseries_file,
            output_file,
            model.grid,
            read_mar_static_columns(first(files), "AREA", model.grid),
        )
        rm(staging_dir; recursive=true)
        @info "Consolidated previously completed prescribed monthly transient" output_file timeseries_file
        return output_file
    end
    isfile(output_file) && error("Monthly output already exists: $output_file")
    spinup = Simulation(
        model;
        forcing=mode_forcing(CONFIG.spinup_file, spinup_loaded),
        years=CONFIG.spinup_years,
        backend=CONFIG.backend,
        write_netcdf=false,
        compute_year_metrics=false,
        name="gris_prescribed_monthly_spinup",
    )
    result = run!(spinup)
    result.status == :complete || error("Spin-up ended with status $(result.status)")
    @info "Spin-up complete" years=CONFIG.spinup_years

    mkpath(staging_dir)
    state = spinup.now
    area_km2 = read_mar_static_columns(first(files), "AREA", model.grid)
    annual_rows = NamedTuple[]
    for (path, year, staged_file) in zip(files, years, staged_files)
        loaded = load_mar(path)
        same_columns(model.grid, loaded.grid) || error("Forcing grid differs from spin-up in $(basename(path)).")
        forcing = mode_forcing(path, loaded)
        before = state_budget(state)
        transient = Simulation(
            model;
            forcing,
            state,
            years=1,
            backend=CONFIG.backend,
            write_netcdf=true,
            netcdf_variables=:monthly,
            netcdf_path=staged_file,
            compute_year_metrics=false,
            name="gris_prescribed_monthly_$(year)",
        )
        result = run!(transient)
        result.status == :complete || error("Transient year $year ended with status $(result.status)")
        state = transient.now
        after = state_budget(state)
        push!(annual_rows, (
            year=year,
            mar=mar_annual_budget(path, model.grid, forcing, area_km2),
            chion_prescribed=chion_annual_budget(before, after, forcing, area_km2),
            chion_parameterized=nothing,
        ))
        @info "Completed prescribed transient year" year
        GC.gc()
    end
    consolidate_monthly_outputs!(model, state, staged_files, output_file)
    write_surface_smb_timeseries(
        timeseries_file,
        output_file,
        model.grid,
        read_mar_static_columns(first(files), "AREA", model.grid),
    )
    annual_csv = write_streamed_budget_comparison(
        joinpath(MONTHLY_CONFIG.output_dir, "yearly_integrated_rates.csv"), annual_rows)
    annual_pdf = plot_results(
        joinpath(MONTHLY_CONFIG.output_dir, "yearly_integrated_rates.pdf"), complete_year_rows(annual_rows))
    open(joinpath(MONTHLY_CONFIG.output_dir, "configuration.csv"), "w") do io
        println(io, "parameter,value")
        println(io, "forcing,$(FORCING_MODE)_MAR")
        println(io, "albedo_scheme,$(get(ENV, "CHION_ALBEDO_SCHEME", FORCING_MODE == :prescribed ? "prescribed" : "dynamic"))")
        println(io, "spinup_years,$(CONFIG.spinup_years)\nsnow_depth_cap_m,22.5")
        write_model_configuration(io, model)
    end
    rm(staging_dir; recursive=true)
    @info "Wrote prescribed monthly transient" output_file timeseries_file annual_csv annual_pdf
    return output_file
end

if abspath(PROGRAM_FILE) == @__FILE__
    monthly_transient_main()
end
