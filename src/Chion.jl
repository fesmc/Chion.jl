module Chion

using Printf
using Dates
using Base.Threads: @threads
using Adapt: Adapt, adapt, @adapt_structure
using CUDA
using KernelAbstractions

include("constants.jl")
include("domain_interfaces.jl")
include("execution_backends.jl")
include("snowpack_domain.jl")
include("step_forcing.jl")
include("step_scratch.jl")
include("column_state_utils.jl")
include("processes/albedo.jl")
include("processes/layer_structure.jl")
include("processes/accumulation.jl")
include("processes/melt.jl")
include("processes/diurnal_shortwave.jl")
include("processes/surface_fluxes.jl")
include("processes/energy_flux.jl")
include("processes/densification.jl")
include("processes/percolation.jl")
include("processes/refreezing.jl")
include("domain_summaries.jl")
include("step_timing.jl")
include("step.jl")
include("step_field_batches.jl")
include("state_access.jl")
include("dataloaders.jl")
include("grid.jl")
include("simulation.jl")
include("api.jl")

export step!
export SnowpackDomain
export SnowpackStepForcing
export SnowpackStepFields
export SnowpackPhysicalConstants
export RunResult
export StepTimingStats, physics
export run!
export LoadedProblem, load_gris_forcing_file_problem
export read_dataset_shapes, read_hdf5_subset, read_hdf5_full
export read_timeslice_2d, read_timeslice_3d
export valid_or, mmwe_day_to_kgm2s
export read_forcing_times, choose_time_index, infer_dt_days
export extract_forcing_file_layers, populate_domain_column_from_forcing_file!
export read_full_timeseries_3d, read_first_available_timeseries_3d
export ColumnarStepWorkspace
export continuous_bottom_deplete!
export update_surface_albedo!
export get_state, print_state
export add_timing!, timing_rows, print_timing_summary

end
