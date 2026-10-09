#!/usr/bin/env julia
# Re-plot yearly_integrated_rates.pdf from the CSV of streamed daily transients,
# excluding the incomplete final forcing year from series and error metrics.
# Usage: julia --project=. replot_gris_transient_yearly_rates.jl <transient output dir>...
include(joinpath(@__DIR__, "run_gris_fully_prescribed_monthly_transient.jl"))

for dir in ARGS
    csv = joinpath(dir, "yearly_integrated_rates.csv")
    isfile(csv) || (@warn "No yearly_integrated_rates.csv" dir; continue)
    plot_results(joinpath(dir, "yearly_integrated_rates.pdf"), complete_year_rows(read_results(csv)))
    println("Re-plotted ", dir)
end
