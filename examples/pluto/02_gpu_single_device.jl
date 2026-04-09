### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82002
begin
	import Pkg
	const NOTEBOOK_DIR = @__DIR__
	const PROJECT_ROOT = normpath(joinpath(NOTEBOOK_DIR, "..", ".."))
	Pkg.activate(PROJECT_ROOT)
	# Uncomment on a fresh machine:
	# Pkg.instantiate()
	PROJECT_ROOT
end

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82003
begin
	try
		using PlutoUI
	catch
		global PlutoUI = (TableOfContents=() -> nothing,)
	end
	using CUDA
	using Chion
	using Statistics
	using Plots
	include(joinpath(PROJECT_ROOT, "examples", "shared", "notebook_helpers.jl"))
	using .ChionNotebookHelpers
end

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82001
md"""
# Chion GPU Single-Device Workflow

This notebook shows the same Chion workflow as the CPU notebook, but with
`backend=:gpu`.

We do the following steps:
1. activate the project
   by adding the repo to `LOAD_PATH`
2. define `physics`
3. define forcing vectors
4. build a `CaseDefinition` from those vectors
5. `build_case(...; run=RunConfig(..., backend=:gpu, cycles=1))`
6. `run_case(...)`
7. move results back to host only when you need host-side inspection

The preferred path is the high-level case API. A later section shows the lower-level device transfer steps for manual batch stepping.
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82005
md"""
## GPU Prerequisites

Make sure you have the appropiate modules installed or loaded. 
- Load the cluster modules first:
  - `module load julia/1.12.2`
  - `module load cuda/13.1.0`
  - `module load hdf5`
  - `module load netcdf-c`
- Request a GPU allocation before opening Pluto.
- Only set `NETCDF_LIB` if you enable NetCDF output.
- This notebook prepends `~/.chion-pluto/notebook-runtime/depot` to `DEPOT_PATH` automatically.

"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82004
PlutoUI.TableOfContents()

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82008
md"""
## Device And Memory Checks

Check `gpu_status.functional` before running GPU cells. `gpu_device.pool_status` gives a quick view of CUDA memory-pool state.
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82006
gpu_status = cuda_preflight()

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82007
gpu_device = device_report()

# ╔═╡ 7f7c10f1-3ae9-49d8-b729-8f8b061520bf
md"""
## Define Physics

Let's define the physics options, i.e., the albedo scheme, the high density densification scheme and the fresh snow density scheme.
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82010
physics = Chion.physics(
	albedo=:dynamic,
	densification=:bessi,
	fresh_snow_density=:constant,
)

# ╔═╡ 504961b5-e1e5-46db-9e2b-9f45524428f7
md"""
Now we define some arbitrary forcings (snowfall, rainfall, shortwave radiation and wind_speed).
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82011
single_forcing_vectors = let
	ntime = 500
	step = collect(1:ntime)
	air_temperature_c = -18.0 .+ 0.02 .* step
	air_temperature_c[120:170] .+= 4.0
	air_temperature_c[330:390] .-= 3.5
	snowfall_mm_day = fill(0.15, ntime)
	snowfall_mm_day[40:80] .+= 1.2
	snowfall_mm_day[200:250] .+= 0.8
	snowfall_mm_day[410:470] .+= 1.0
	rainfall_mm_day = zeros(Float64, ntime)
	rainfall_mm_day[260:320] .= 0.6
	rainfall_mm_day[321:380] .= 1.1
	shortwave_down = fill(140.0, ntime)
	shortwave_down[150:260] .+= 50.0
	shortwave_down[261:360] .+= 85.0
	wind_speed = 4.5 .+ 0.4 .* sin.(0.05 .* step)
	(
		dt_days=fill(1.0, ntime),
		air_temperature_c=air_temperature_c,
		snowfall_mm_day=snowfall_mm_day,
		rainfall_mm_day=rainfall_mm_day,
		shortwave_down=shortwave_down,
		wind_speed=wind_speed,
	)
end

# ╔═╡ a756c300-41f0-4817-a197-616efa8a0487
md"""
Let's build the single-column case. The maximum amount of vertical layers is 5 (ntot). One column means nx=1 and ny=1. We initalise the column with some surface mass, temperature and density. 

%%% Currently, max mass is fixed at 15 layers -> needs to change. 

"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82012
single_data = Chion.prescribed_case(
	physics=physics,
	ntot=5,
	nx=1,
	ny=1,
	initial_surface_mass=250.0,
	initial_density=320.0,
	initial_temperature_c=-12.0,
	input_label="single_column_vectors",
	; single_forcing_vectors...,
)

# ╔═╡ c29af6d1-2ca9-40b3-a51d-43e54f73224c
gpu_single_forcing_plot = plots_available() ?
	forcing_timeseries_plot(single_data.forcing; idx=1, title_prefix="GPU single-column") :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ 911569b5-812c-47a2-a8e8-f72dcde0a5dc
md"""
We tell the computer, that we want to run on the GPU and define what forcing data to use. Furthermore, we tell the computer where to save output and if we even want output. 
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82013
begin
	single_case = Chion.build_case(
		single_data;
		run=Chion.RunConfig(
			name="gpu_single_column",
			backend=:gpu,
			output_dir=output_dir_for("02_gpu_single_device"),
			write_outputs=false,
			write_netcdf=false,
			cycles=1,
			history_stride=1,
		),
	)
	gpu_single_run = gpu_status.functional ? run_case_capture(single_case) : nothing
	single_case
end

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82014
gpu_single_summary = isnothing(gpu_single_run) ? (
	message="CUDA is not functional in this session. Read the troubleshooting section below.",
) : (
	case=single_case,
	metadata=single_data.metadata,
	status=gpu_single_run.result.status,
	forcing_steps=length(single_data.forcing.time_values),
	history_records=length(gpu_single_run.result.history),
	summary_path=gpu_single_run.result.summary_path,
	history_csv_path=gpu_single_run.result.history_csv_path,
	column_summary=summarize_column(gpu_single_run.result, 1),
)

# ╔═╡ e3508184-17eb-4d2f-a909-41caac42d417
md"""
## Single-Column GPU Plots
"""

# ╔═╡ ed47e84a-08f9-4d94-ad7b-f6f9c1e38413
gpu_single_profile_plot = isnothing(gpu_single_run) ? gpu_single_summary.message :
	plots_available() ?
	column_profile_plot(gpu_single_run.result, 1; title="GPU single-column") :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ 74f2e022-fcd8-4251-831c-a566dbea0697
md""" 
Let's increase the amount of snow that falls, i.e. constant = 3 mmWE/day and run the column again.

"""

# ╔═╡ 026429de-3c0e-42de-b2ca-39116cf64bbd
begin
	single_forcing_vectors_constant_snow = let
		ntime = 500
		step = collect(1:ntime)
		air_temperature_c = -18.0 .+ 0.02 .* step
		air_temperature_c[120:170] .+= 4.0
		air_temperature_c[330:390] .-= 3.5
		snowfall_mm_day = fill(3, ntime)
		rainfall_mm_day = zeros(Float64, ntime)
		rainfall_mm_day[260:320] .= 0.6
		rainfall_mm_day[321:380] .= 1.1
		shortwave_down = fill(140.0, ntime)
		shortwave_down[150:260] .+= 50.0
		shortwave_down[261:360] .+= 85.0
		wind_speed = 4.5 .+ 0.4 .* sin.(0.05 .* step)
		(
			dt_days=fill(1.0, ntime),
			air_temperature_c=air_temperature_c,
			snowfall_mm_day=snowfall_mm_day,
			rainfall_mm_day=rainfall_mm_day,
			shortwave_down=shortwave_down,
			wind_speed=wind_speed,
		)
	end
	constant_snow = Chion.prescribed_case(
		physics=physics,
		ntot=5,
		nx=1,
		ny=1,
		initial_surface_mass=250.0,
		initial_density=320.0,
		initial_temperature_c=-12.0,
		input_label="single_column_vectors",
		; single_forcing_vectors_constant_snow...,
	)
	
end

# ╔═╡ c180d37e-0277-4c42-9464-7c06c8f7ad6a
md"""
As we see, the amount of snow increases and we have 5 active layers now with the bottom layer accumulating the snow mass. 
"""

# ╔═╡ 3c726710-1c66-42c3-82bb-eaeaebeb3920
begin
	const_snow_case = Chion.build_case(
		constant_snow;
		run=Chion.RunConfig(
			name="gpu_single_column",
			backend=:gpu,
			output_dir=output_dir_for("02_gpu_single_device"),
			write_outputs=true,
			write_netcdf=false,
			cycles=1,
			history_stride=1,
		),
	)
	gpu_snow_run = gpu_status.functional ? run_case_capture(const_snow_case) : nothing
	const_snow_case

	constant_snow_plot = isnothing(gpu_snow_run) ? "GPU run skipped because CUDA is unavailable." :
		plots_available() ?
		column_profile_plot(gpu_snow_run.result, 1; title="GPU single-column") :
		"Plots.jl is not available in this Pluto session."
end

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82015
md"""
## Multi-Column GPU Example

Let's run a multi-column example, this is where the GPU really helps. 
First, we define a simple grid with 40 times 50 cells. 
Then we initialise the cells with some mass, density and temperature again and build the case. 

As a sanity check, we only vary the snowfall spatially as function of the normalised coordinates snow = (x-0.5)^2 + (y-0.5)^2 here. This means, that there is the most snow at the edges and no snow in the center.
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82016
multi_grid = (
	nx=40,
	ny=50,
)

# ╔═╡ f5d81564-8ac7-4f12-ba2e-e2ef0b390e12
begin
	gridfield_to_columns(A) = vec(permutedims(A, (2, 1)))

	gridcube_to_columns(A) = begin
		@assert size(A, 1) == multi_grid.ny
		@assert size(A, 2) == multi_grid.nx
		reshape(
			permutedims(A, (2, 1, 3)),
			multi_grid.nx * multi_grid.ny,
			size(A, 3),
		)
	end
end

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82017
multi_data = let
	nx = multi_grid.nx
	ny = multi_grid.ny
	ntime = length(single_forcing_vectors.dt_days)

	multi_surface_mass = fill(220.0, nx * ny)

	air_temperature_c_cube = Array{Float64}(undef, ny, nx, ntime)
	snowfall_mm_day_cube = Array{Float64}(undef, ny, nx, ntime)
	rainfall_mm_day_cube = Array{Float64}(undef, ny, nx, ntime)
	shortwave_down_cube = Array{Float64}(undef, ny, nx, ntime)
	wind_speed_cube = Array{Float64}(undef, ny, nx, ntime)

	for t in 1:ntime, j in 1:ny, i in 1:nx
		xfrac = nx == 1 ? 0.0 : (i - 1) / (nx - 1)
		yfrac = ny == 1 ? 0.0 : (j - 1) / (ny - 1)

		air_temperature_c_cube[j, i, t] = (xfrac-1)
		
		snowfall_mm_day_cube[j, i, t] = (xfrac-0.5)^2 + (yfrac-0.5)^2
		rainfall_mm_day_cube[j, i, t] = 0
		shortwave_down_cube[j, i, t] = 300
		wind_speed_cube[j, i, t] = 5
	end

	Chion.prescribed_case(
		physics=physics,
		ntot=15,
		nx=nx,
		ny=ny,
		initial_surface_mass=multi_surface_mass,
		initial_density=320.0,
		initial_temperature_c=-12.0,
		input_label="multi_column_spatiotemporal_forcing",
		dt_days=single_forcing_vectors.dt_days,
		air_temperature_c=gridcube_to_columns(air_temperature_c_cube),
		snowfall_mm_day=gridcube_to_columns(snowfall_mm_day_cube),
		rainfall_mm_day=gridcube_to_columns(rainfall_mm_day_cube),
		shortwave_down=gridcube_to_columns(shortwave_down_cube),
		wind_speed=gridcube_to_columns(wind_speed_cube),
	)

end



# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82018
multi_case = Chion.build_case(
	multi_data;
	run=Chion.RunConfig(
		name="gpu_multi_column",
		backend=:gpu,
		output_dir=output_dir_for("02_gpu_single_device"),
		write_outputs=false,
		write_netcdf=false,
		cycles=10,
		history_stride=1,
	),
)

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82019
gpu_multi_run = gpu_status.functional ? run_case_capture(multi_case) : nothing

# ╔═╡ 76a4a56d-ff88-46ef-bf61-a4472eeb4ba5
gpu_multi_sample_indices = unique([
	1,
	cld(multi_data.metadata.ncol, 2),
	multi_data.metadata.ncol,
])

# ╔═╡ 5df2f2c8-6cc1-4d7a-b2f5-e77db43a9c82
gpu_multi_summary = isnothing(gpu_multi_run) ? (
	message="GPU run skipped because CUDA is unavailable.",
) : (
	case=multi_case,
	metadata=multi_data.metadata,
	status=gpu_multi_run.result.status,
	forcing_steps=length(multi_data.forcing.time_values),
	history_records=length(gpu_multi_run.result.history),
	summary_path=gpu_multi_run.result.summary_path,
	history_csv_path=gpu_multi_run.result.history_csv_path,
	sample_columns=[summarize_column(gpu_multi_run.result, idx) for idx in gpu_multi_sample_indices],
)

# ╔═╡ 3c6efec4-12ab-4a0d-8f3d-372477851cb6
gpu_multi_thickness_plot = isnothing(gpu_multi_run) ? gpu_multi_summary.message :
	plots_available() ?
	layout_heatmap_plot(
		multi_data.layout,
		domain_metric_values(gpu_multi_run.result, :thickness);
		title="GPU multi-column final thickness",
		unit="m",
		color=:ice,
	) :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ a22d9552-5a97-4a33-a98c-75d8e66ef7ac
md"""
As expected, there is no snow in the center and the snowdepth distribution follows the snowfall since we do not have any melt etc.
"""

# ╔═╡ b8b1a062-939f-444d-9c86-9f90959bc2d8
md"""
## Comparison with CPU

We can compare the runtime and results with a CPU execution instead.
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82020
md"""
## What Changes Under The Hood?

`run_case(...)` handles the device transfer for you. The lower-level batch workflow makes those steps explicit:

1. move the domain to GPU with `gpu_domain`
2. move forcing arrays to GPU storage with `adapt(CUDA.CuArray, ...)`
3. allocate a `ColumnarStepWorkspace` on the same device
4. call `step!`
5. move the domain back with `cpu_domain` before host-side inspection
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82021
manual_gpu_step = if gpu_status.functional
	case_data = Chion.prescribed_case(
		physics=physics,
		ntot=5,
		nx=4,
		ny=3,
		initial_surface_mass=240.0,
		initial_density=320.0,
		initial_temperature_c=-10.0,
		input_label="manual_gpu_step_vectors",
		dt_days=fill(1.0, 8),
		air_temperature_c=collect(range(-16.0, -6.0; length=8)),
		snowfall_mm_day=[0.0, 0.4, 1.0, 0.8, 0.2, 0.0, 0.0, 0.0],
		rainfall_mm_day=[0.0, 0.0, 0.0, 0.0, 0.2, 0.6, 0.4, 0.1],
		shortwave_down=[120.0, 130.0, 150.0, 180.0, 210.0, 220.0, 200.0, 160.0],
		wind_speed=[4.2, 4.4, 4.8, 5.0, 5.2, 5.0, 4.7, 4.5],
	)
	domain_cpu = deepcopy(case_data.domain)
	forcing_cpu = case_data.forcing
	step_fields_cpu = Chion.SnowpackModel.SnowpackStepFields(forcing_cpu)
	domain_gpu = Chion.SnowpackModel.gpu_domain(deepcopy(domain_cpu))
	step_fields_gpu = Chion.SnowpackModel.adapt(CUDA.CuArray, step_fields_cpu)
	workspace_gpu = Chion.ColumnarStepWorkspace(domain_gpu)
	Chion.step!(domain_gpu, step_fields_gpu, 1, workspace_gpu)
	CUDA.synchronize()
	host_domain = Chion.SnowpackModel.cpu_domain(domain_gpu)
	summary = Chion.SnowpackModel.summarize_domain_state(host_domain; backend=:threads)
	(
		device_domain_type=typeof(domain_gpu.mass),
		device_forcing_type=typeof(step_fields_gpu.air_temperature),
		host_domain_type=typeof(host_domain.mass),
		mean_thickness_after_one_step=mean(summary.thickness),
	)
else
	(
		message="Manual GPU stepping skipped because CUDA is unavailable.",
	)
end

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82022
md"""
## CPU vs GPU Summary
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82023
[
	(
		mode="CPU case API",
		backend=:cpu,
		domain_storage="Array",
		manual_transfer="none",
		workspace="threaded_workspaces(domain)",
		supported=true,
	),
	(
		mode="GPU case API",
		backend=:gpu,
		domain_storage="CuArray after transfer",
		manual_transfer="handled inside run_case(...)",
		workspace="ColumnarStepWorkspace(domain_gpu)",
		supported=true,
	),
	(
		mode="Manual GPU batch step",
		backend=:gpu,
		domain_storage="CuArray",
		manual_transfer="required",
		workspace="ColumnarStepWorkspace(domain_gpu)",
		supported=true,
	),
]

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82024
md"""
## Troubleshooting

- `CUDA.functional()==false`: confirm that a GPU is allocated and `cuda/13.1.0` is loaded.
- Read-only depot or failed precompile: set `JULIA_DEPOT_PATH=\$HOME/.chion-pluto/manual-depot:\$HOME/.julia`.
- Wrong device/module combination: reload Julia after changing CUDA modules so `CUDA.jl` sees the correct driver/runtime pair.
- OOM on larger problems: reduce forcing length, `ncol`, or `ntot`, or clear cached allocations with `CUDA.reclaim()`.
- Host-side inspection fails on device arrays: call `result_domain_cpu(result)` or `cpu_domain(domain_gpu)` before `get_state`.
- NetCDF write errors: only enable NetCDF when `NETCDF_LIB` points to a valid library and a grid layout is available.
"""

# ╔═╡ Cell order:
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82001
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82002
# ╟─5d9fd6b6-0c79-11ef-86e7-174c14f82005
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82003
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82004
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82008
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82006
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82007
# ╟─7f7c10f1-3ae9-49d8-b729-8f8b061520bf
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82010
# ╟─504961b5-e1e5-46db-9e2b-9f45524428f7
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82011
# ╠═c29af6d1-2ca9-40b3-a51d-43e54f73224c
# ╟─a756c300-41f0-4817-a197-616efa8a0487
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82012
# ╟─911569b5-812c-47a2-a8e8-f72dcde0a5dc
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82013
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82014
# ╟─e3508184-17eb-4d2f-a909-41caac42d417
# ╠═ed47e84a-08f9-4d94-ad7b-f6f9c1e38413
# ╟─74f2e022-fcd8-4251-831c-a566dbea0697
# ╠═026429de-3c0e-42de-b2ca-39116cf64bbd
# ╟─c180d37e-0277-4c42-9464-7c06c8f7ad6a
# ╠═3c726710-1c66-42c3-82bb-eaeaebeb3920
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82015
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82016
# ╠═f5d81564-8ac7-4f12-ba2e-e2ef0b390e12
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82017
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82018
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82019
# ╠═76a4a56d-ff88-46ef-bf61-a4472eeb4ba5
# ╠═5df2f2c8-6cc1-4d7a-b2f5-e77db43a9c82
# ╠═3c6efec4-12ab-4a0d-8f3d-372477851cb6
# ╟─a22d9552-5a97-4a33-a98c-75d8e66ef7ac
# ╠═b8b1a062-939f-444d-9c86-9f90959bc2d8
# ╟─5d9fd6b6-0c79-11ef-86e7-174c14f82020
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82021
# ╟─5d9fd6b6-0c79-11ef-86e7-174c14f82022
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82023
# ╟─5d9fd6b6-0c79-11ef-86e7-174c14f82024
