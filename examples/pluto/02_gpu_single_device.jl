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

Same workflow as the CPU notebook, but with `backend=:gpu`.

1. construct a `SnowpackGrid` and a `BESSIModel`
2. define a `SnowpackForcing`
3. build a `Simulation` with `backend=:gpu` and call `run!(simulation)`
4. move domain back to host for inspection via `cpu_domain(model.domain)`
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82005
md"""
## GPU Prerequisites

- `module load julia/1.12.2`
- `module load cuda/13.1.0`
- Request a GPU allocation before opening Pluto.
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82004
PlutoUI.TableOfContents()

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82008
md"""
## Device And Memory Checks
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82006
gpu_status = cuda_preflight()

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82007
gpu_device = device_report()

# ╔═╡ 7f7c10f1-3ae9-49d8-b729-8f8b061520bf
md"""
## Define Forcing Vectors
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82011
single_forcing_vectors = let
	ntime = 500
	step = collect(1:ntime)
	air_temperature_c = -18.0 .+ 0.02 .* step
	air_temperature_c[120:170] .+= 4.0
	air_temperature_c[330:390] .-= 3.5
	snowfall_mm_day = fill(0.15, ntime)
	snowfall_mm_day[40:80]   .+= 1.2
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
## Single-Column GPU Example
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82012
begin
	single_grid   = SnowpackGrid(CPU(), 1)
	single_model  = seed_domain!(
		BESSIModel(single_grid; Ntot=5);
		surface_mass=250.0, density=320.0, temperature_c=-12.0,
	)
	single_forcing = SnowpackForcing(
		dt_days          = single_forcing_vectors.dt_days,
		air_temperature_c = single_forcing_vectors.air_temperature_c,
		snowfall_mm_day  = single_forcing_vectors.snowfall_mm_day,
		rainfall_mm_day  = single_forcing_vectors.rainfall_mm_day,
		shortwave_down   = single_forcing_vectors.shortwave_down,
		wind_speed       = single_forcing_vectors.wind_speed,
		ncol             = 1,
	)
end

# ╔═╡ c29af6d1-2ca9-40b3-a51d-43e54f73224c
gpu_single_forcing_plot = plots_available() ?
	forcing_timeseries_plot(single_forcing; idx=1, title_prefix="GPU single-column") :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ 911569b5-812c-47a2-a8e8-f72dcde0a5dc
md"""
## Run on GPU
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82013
gpu_single_result = if gpu_status.functional
	sim = Simulation(single_model;
		forcing=single_forcing, cycles=1, backend=:gpu,
		save=Symbol[], output_dir=output_dir_for("02_gpu_single_device"),
		write_outputs=false, history_stride=1,
	)
	run!(sim)
else
	nothing
end

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82014
gpu_single_summary = isnothing(gpu_single_result) ? (
	message="CUDA is not functional in this session.",
) : (
	status          = gpu_single_result.status,
	forcing_steps   = length(single_forcing.time_values),
	history_records = length(gpu_single_result.history),
	column_summary  = summarize_column(single_model, 1),
)

# ╔═╡ e3508184-17eb-4d2f-a909-41caac42d417
md"""
## Single-Column GPU Plots
"""

# ╔═╡ ed47e84a-08f9-4d94-ad7b-f6f9c1e38413
gpu_single_profile_plot = isnothing(gpu_single_result) ? "GPU run skipped." :
	plots_available() ?
	column_profile_plot(single_model, 1; title="GPU single-column") :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ 74f2e022-fcd8-4251-831c-a566dbea0697
md"""
## Constant-Snowfall Variant
"""

# ╔═╡ 026429de-3c0e-42de-b2ca-39116cf64bbd
begin
	const_snow_vectors = let
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
		(dt_days=fill(1.0, ntime), air_temperature_c=air_temperature_c,
		 snowfall_mm_day=snowfall_mm_day, rainfall_mm_day=rainfall_mm_day,
		 shortwave_down=shortwave_down, wind_speed=wind_speed)
	end
	const_snow_grid   = SnowpackGrid(CPU(), 1)
	const_snow_model  = seed_domain!(
		BESSIModel(const_snow_grid; Ntot=5);
		surface_mass=250.0, density=320.0, temperature_c=-12.0,
	)
	const_snow_forcing = SnowpackForcing(
		dt_days          = const_snow_vectors.dt_days,
		air_temperature_c = const_snow_vectors.air_temperature_c,
		snowfall_mm_day  = const_snow_vectors.snowfall_mm_day,
		rainfall_mm_day  = const_snow_vectors.rainfall_mm_day,
		shortwave_down   = const_snow_vectors.shortwave_down,
		wind_speed       = const_snow_vectors.wind_speed,
		ncol             = 1,
	)
end

# ╔═╡ 3c726710-1c66-42c3-82bb-eaeaebeb3920
const_snow_result = if gpu_status.functional
	sim = Simulation(const_snow_model;
		forcing=const_snow_forcing, cycles=1, backend=:gpu,
		save=Symbol[], output_dir=output_dir_for("02_gpu_single_device"),
		write_outputs=false, history_stride=1,
	)
	run!(sim)
else
	nothing
end

# ╔═╡ c180d37e-0277-4c42-9464-7c06c8f7ad6a
constant_snow_plot = isnothing(const_snow_result) ? "GPU run skipped." :
	plots_available() ?
	column_profile_plot(const_snow_model, 1; title="GPU constant-snow") :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82015
md"""
## Multi-Column GPU Example
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82016
multi_dims = (nx=40, ny=50)

# ╔═╡ f5d81564-8ac7-4f12-ba2e-e2ef0b390e12
begin
	gridcube_to_columns(A, nx, ny) = reshape(
		permutedims(A, (2, 1, 3)), nx * ny, size(A, 3),
	)
end

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82017
begin
	nx, ny = multi_dims.nx, multi_dims.ny
	ntime  = length(single_forcing_vectors.dt_days)

	air_cube  = Array{Float64}(undef, ny, nx, ntime)
	snow_cube = Array{Float64}(undef, ny, nx, ntime)
	rain_cube = zeros(Float64, ny, nx, ntime)
	sw_cube   = fill(300.0, ny, nx, ntime)
	wind_cube = fill(5.0, ny, nx, ntime)

	for t in 1:ntime, j in 1:ny, i in 1:nx
		xfrac = nx == 1 ? 0.0 : (i - 1) / (nx - 1)
		yfrac = ny == 1 ? 0.0 : (j - 1) / (ny - 1)
		air_cube[j, i, t]  = xfrac - 1.0
		snow_cube[j, i, t] = (xfrac - 0.5)^2 + (yfrac - 0.5)^2
	end

	multi_grid   = regular_grid(nx, ny)
	multi_model  = seed_domain!(
		BESSIModel(multi_grid; Ntot=15);
		surface_mass=fill(220.0, nx * ny), density=320.0, temperature_c=-12.0,
	)
	multi_forcing = SnowpackForcing(
		dt_days          = single_forcing_vectors.dt_days,
		air_temperature_c = gridcube_to_columns(air_cube, nx, ny),
		snowfall_mm_day  = gridcube_to_columns(snow_cube, nx, ny),
		rainfall_mm_day  = gridcube_to_columns(rain_cube, nx, ny),
		shortwave_down   = gridcube_to_columns(sw_cube, nx, ny),
		wind_speed       = gridcube_to_columns(wind_cube, nx, ny),
		ncol             = nx * ny,
	)
end

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82019
gpu_multi_result = if gpu_status.functional
	sim = Simulation(multi_model;
		forcing=multi_forcing, cycles=10, backend=:gpu,
		save=Symbol[], output_dir=output_dir_for("02_gpu_single_device"),
		write_outputs=false, history_stride=1,
	)
	run!(sim)
else
	nothing
end

# ╔═╡ 5df2f2c8-6cc1-4d7a-b2f5-e77db43a9c82
gpu_multi_summary = isnothing(gpu_multi_result) ? (
	message="GPU run skipped because CUDA is unavailable.",
) : let
	ncol = nx * ny
	sample_idx = unique([1, cld(ncol, 2), ncol])
	(
		status          = gpu_multi_result.status,
		forcing_steps   = length(multi_forcing.time_values),
		history_records = length(gpu_multi_result.history),
		sample_columns  = [summarize_column(multi_model, i) for i in sample_idx],
	)
end

# ╔═╡ 3c6efec4-12ab-4a0d-8f3d-372477851cb6
gpu_multi_thickness_plot = isnothing(gpu_multi_result) ? "GPU run skipped." :
	plots_available() ?
	layout_heatmap_plot(
		multi_grid,
		domain_metric_values(multi_model, :thickness);
		title="GPU multi-column final thickness",
		unit="m",
		color=:ice,
	) :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82020
md"""
## Manual Low-Level GPU Stepping

`run!(simulation; backend=:gpu)` handles device transfers automatically.
The code below shows the explicit steps for manual batch stepping.
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82021
manual_gpu_step = if gpu_status.functional
	g = SnowpackGrid(CPU(), 12)
	m = seed_domain!(
		BESSIModel(g; Ntot=5);
		surface_mass=240.0, density=320.0, temperature_c=-10.0,
	)
	f_manual = SnowpackForcing(
		dt_days          = fill(1.0, 8),
		air_temperature_c = collect(range(-16.0, -6.0; length=8)),
		snowfall_mm_day  = [0.0, 0.4, 1.0, 0.8, 0.2, 0.0, 0.0, 0.0],
		rainfall_mm_day  = [0.0, 0.0, 0.0, 0.0, 0.2, 0.6, 0.4, 0.1],
		shortwave_down   = [120.0, 130.0, 150.0, 180.0, 210.0, 220.0, 200.0, 160.0],
		wind_speed       = [4.2, 4.4, 4.8, 5.0, 5.2, 5.0, 4.7, 4.5],
		ncol             = 12,
	)
	step_fields_cpu = f_manual
	domain_gpu      = Chion.gpu_domain(deepcopy(m.domain))
	step_fields_gpu = Chion.adapt(CUDA.CuArray, step_fields_cpu)
	workspace_gpu   = Chion.ColumnarStepWorkspace(domain_gpu)
	Chion.step!(domain_gpu, step_fields_gpu, 1, workspace_gpu)
	CUDA.synchronize()
	host_domain = Chion.cpu_domain(domain_gpu)
	summary = Chion.summarize_domain_state(host_domain)
	(
		device_domain_type  = typeof(domain_gpu.mass),
		device_forcing_type = typeof(step_fields_gpu.air_temperature),
		host_domain_type    = typeof(host_domain.mass),
		mean_thickness      = mean(summary.thickness),
	)
else
	(message="Manual GPU stepping skipped because CUDA is unavailable.",)
end

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82022
md"""
## CPU vs GPU Summary
"""

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82023
[
	(mode="CPU Simulation",   backend=:cpu, domain_storage="Array",            manual_transfer="none",                         supported=true),
	(mode="GPU Simulation",   backend=:gpu, domain_storage="CuArray (auto)",   manual_transfer="handled inside run!",           supported=true),
	(mode="Manual GPU step",  backend=:gpu, domain_storage="CuArray (manual)", manual_transfer="gpu_domain / cpu_domain",       supported=true),
]

# ╔═╡ 5d9fd6b6-0c79-11ef-86e7-174c14f82024
md"""
## Troubleshooting

- `CUDA.functional()==false`: confirm GPU allocation and `cuda/13.1.0` is loaded.
- Read-only depot: `export JULIA_DEPOT_PATH=\$HOME/.chion-pluto/manual-depot:\$HOME/.julia`.
- OOM: reduce forcing length, `ncol`, or `Ntot`, or call `CUDA.reclaim()`.
- Host inspection on device arrays: call `cpu_domain(model.domain)` before `get_state`.
"""

# ╔═╡ Cell order:
# ╟─5d9fd6b6-0c79-11ef-86e7-174c14f82001
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82002
# ╟─5d9fd6b6-0c79-11ef-86e7-174c14f82005
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82003
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82004
# ╟─5d9fd6b6-0c79-11ef-86e7-174c14f82008
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82006
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82007
# ╟─7f7c10f1-3ae9-49d8-b729-8f8b061520bf
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
# ╠═3c726710-1c66-42c3-82bb-eaeaebeb3920
# ╠═c180d37e-0277-4c42-9464-7c06c8f7ad6a
# ╟─5d9fd6b6-0c79-11ef-86e7-174c14f82015
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82016
# ╠═f5d81564-8ac7-4f12-ba2e-e2ef0b390e12
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82017
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82019
# ╠═5df2f2c8-6cc1-4d7a-b2f5-e77db43a9c82
# ╠═3c6efec4-12ab-4a0d-8f3d-372477851cb6
# ╟─5d9fd6b6-0c79-11ef-86e7-174c14f82020
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82021
# ╟─5d9fd6b6-0c79-11ef-86e7-174c14f82022
# ╠═5d9fd6b6-0c79-11ef-86e7-174c14f82023
# ╟─5d9fd6b6-0c79-11ef-86e7-174c14f82024
