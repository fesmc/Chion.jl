### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1002
begin
	import Pkg
	const NOTEBOOK_DIR = @__DIR__
	const PROJECT_ROOT = normpath(joinpath(NOTEBOOK_DIR, "..", ".."))
	Pkg.activate(PROJECT_ROOT)
	# Uncomment on a fresh machine:
	# Pkg.instantiate()
	PROJECT_ROOT
end

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1003
begin
	try
		using PlutoUI
	catch
		global PlutoUI = (TableOfContents=() -> nothing,)
	end
	using Chion
	include(joinpath(PROJECT_ROOT, "examples", "shared", "notebook_helpers.jl"))
	using .ChionNotebookHelpers
end

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1001
md"""
# Chion CPU Workflows

This notebook shows the standard Chion Pluto workflow on CPU using
user-supplied forcing vectors:

1. activate the project
2. define `physics`
3. define forcing vectors
4. build case data from those vectors
5. `build_case(...; backend=:cpu, run_forcing_once=true)`
6. `run_case(...)`
7. inspect the result and saved outputs

The notebooks run one pass through the provided forcing. There is no synthetic
cycle parameter to tune.
"""

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1004
PlutoUI.TableOfContents()

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1005
md"""
## Prerequisites

- Start Pluto from a Julia environment that has `Pluto` installed.
- This notebook uses the package project at `$(PROJECT_ROOT)`.
- On the cluster, load Julia first: `module load julia/1.12.2`.
- If precompilation fails because the default depot is read-only, set:
  `export JULIA_DEPOT_PATH=/tmp/chion-pluto:\$HOME/.julia`

The notebook itself depends on the project packages plus `PlutoUI`.
"""

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1006
physics = build_physics(
	albedo_scheme=:dynamic,
	densification=:bessi,
	fresh_snow_density=:constant,
)

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1007
md"""
## Single-Column Example

Replace the example vectors below with your own forcing arrays. The only hard
requirement is that all forcing vectors have the same length.
"""

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1008
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

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1009
single_data = build_prescribed_case_data(
	physics=physics,
	ntot=5,
	nx=1,
	ny=1,
	initial_surface_mass=250.0,
	initial_density=320.0,
	initial_temperature_c=-12.0,
	forcing_label="single_column_vectors",
	; single_forcing_vectors...,
)

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1010
single_case = Chion.build_case(
	single_data;
	name="cpu_single_column",
	backend=:cpu,
	out_dir=output_dir_for("01_cpu_workflows"),
	write_outputs=true,
	write_netcdf=false,
	run_forcing_once=true,
	cycle_metrics_stride=1,
)

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1011
single_run = run_case_capture(single_case)

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1012
let
	single_state = Chion.get_state(result_domain_cpu(single_run.result), 1)
	(
		case=single_case,
		metadata=single_data.metadata,
		status=single_run.result.status,
		forcing_steps=length(single_data.forcing.time_values),
		history_records=length(single_run.result.history),
		summary_path=single_run.result.summary_path,
		history_csv_path=single_run.result.history_csv_path,
		column_summary=summarize_column(single_run.result, 1),
		state_keys=collect(keys(single_state)),
	)
end

# ╔═╡ 54d75be0-65f9-4e5c-b178-8a0d4272b677
md"""
## Single-Column Plots
"""

# ╔═╡ 7012df88-b561-4962-bf36-d662f5d54bf3
single_forcing_plot = plots_available() ?
	forcing_timeseries_plot(single_data.forcing; idx=1, title_prefix="CPU single-column") :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ 4c2846e2-7426-4b87-8e9f-c5b8bcb8aec7
single_forcing_summary = (
	forcing_steps=length(single_forcing_vectors.dt_days),
	dt_days_preview=single_forcing_vectors.dt_days[1:5],
	air_temperature_c_range=(minimum(single_forcing_vectors.air_temperature_c), maximum(single_forcing_vectors.air_temperature_c)),
)

# ╔═╡ 09351cdb-5748-4935-8f0a-e5c931fd716d
single_profile_plot = plots_available() ?
	column_profile_plot(single_run.result, 1; title="CPU single-column") :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1013
md"""
## Multi-Column Example

The workflow is identical. We keep the same forcing vectors, choose an explicit
grid size, and broadcast the forcing across all columns.
"""

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1014
multi_grid = (
	nx=40,
	ny=50,
)

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1015
multi_data = let
	multi_surface_mass = [
		220.0 + 35.0 * ((i - 1) / max(multi_grid.nx - 1, 1)) + 20.0 * ((j - 1) / max(multi_grid.ny - 1, 1))
		for j in 1:multi_grid.ny for i in 1:multi_grid.nx
	]
	build_prescribed_case_data(
		physics=physics,
		ntot=5,
		nx=multi_grid.nx,
		ny=multi_grid.ny,
		initial_surface_mass=multi_surface_mass,
		initial_density=320.0,
		initial_temperature_c=-12.0,
		forcing_label="multi_column_vectors",
		; single_forcing_vectors...,
	)
end

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1016
multi_case = Chion.build_case(
	multi_data;
	name="cpu_multi_column",
	backend=:cpu,
	out_dir=output_dir_for("01_cpu_workflows"),
	write_outputs=true,
	write_netcdf=false,
	run_forcing_once=true,
	cycle_metrics_stride=1,
)

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1017
multi_run = run_case_capture(multi_case)

# ╔═╡ 6e837b58-a3bc-4d9d-b6f5-6d8c2b8d43d1
multi_sample_indices = unique([
	1,
	cld(multi_data.metadata.ncol, 2),
	multi_data.metadata.ncol,
])

# ╔═╡ 0b0dd5b7-3044-4aae-b512-cf3d4e6d5dc9
(
	case=multi_case,
	metadata=multi_data.metadata,
	status=multi_run.result.status,
	forcing_steps=length(multi_data.forcing.time_values),
	history_records=length(multi_run.result.history),
	summary_path=multi_run.result.summary_path,
	history_csv_path=multi_run.result.history_csv_path,
	sample_columns=[summarize_column(multi_run.result, idx) for idx in multi_sample_indices],
)

# ╔═╡ bbe5850d-eb9b-4c8f-a6c4-a87836f7b33d
md"""
## Multi-Column Plot
"""

# ╔═╡ aaad1caa-0e49-4d1e-87d5-2ad2d3614acf
multi_thickness_plot = plots_available() ?
	layout_heatmap_plot(
		multi_data.layout,
		domain_metric_values(multi_run.result, :thickness);
		title="CPU multi-column final thickness",
		unit="m",
		color=:ice,
	) :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1018
md"""
## Physics-Scheme Variant Comparison

These short runs keep the execution path the same while changing process-related
physics options.
"""

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1019
scheme_comparison = let
	variants = [
		(label="default", kwargs=(; albedo_scheme=:dynamic, densification=:bessi, fresh_snow_density=:constant)),
		(label="htessel_like", kwargs=(; albedo_scheme=:constant, densification=:htessel, fresh_snow_density=:parameterized)),
	]
	[
			let
				variant_physics = build_physics(; variant.kwargs...)
				case_data = build_prescribed_case_data(
					physics=variant_physics,
					ntot=5,
					nx=1,
					ny=1,
					forcing_label="comparison_" * variant.label,
					; single_forcing_vectors...,
				)
				case = Chion.build_case(
					case_data;
					name="cpu_" * variant.label,
					backend=:cpu,
					out_dir=output_dir_for("01_cpu_workflows"),
					write_outputs=false,
					write_netcdf=false,
					run_forcing_once=true,
				)
			run = run_case_capture(case)
			(
				label=variant.label,
				case=case,
				status=run.result.status,
				last_record=last(run.result.history),
				column_summary=summarize_column(run.result, 1),
			)
		end
		for variant in variants
	]
end

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1020
md"""
## Output Inspection

- `single_run.result.history` and `multi_run.result.history` each contain one summary record because the notebook runs a single pass through the forcing.
- `run_forcing_once=true` is just a user-facing shortcut for "use the provided forcing vectors once".
- `single_run.result.summary_path` and `history_csv_path` point to saved text/CSV outputs.
- `Chion.get_state(result_domain_cpu(result), idx)` gives a host-side snapshot you can inspect or plot.

The raw captured logs are also available:
"""

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1021
(
	single_log_preview=join(split(single_run.log, '\n')[1:min(8, length(split(single_run.log, '\n')))], '\n'),
	multi_log_preview=join(split(multi_run.log, '\n')[1:min(8, length(split(multi_run.log, '\n')))], '\n'),
)

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1022
md"""
## Performance Notes

- `backend=:cpu` dispatches to the threaded CPU backend.
- Interactive notebooks should keep `nx`, `ny`, forcing length, and `ntot` small enough for the machine you are using.
- Writing text outputs is cheap; enabling NetCDF adds extra I/O and requires a valid grid layout and NetCDF library.
- For quick iteration, set `write_outputs=false` and shorten the forcing vectors.
"""

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1023
md"""
## Troubleshooting

- If `using Chion` fails in a restricted shell, set `JULIA_DEPOT_PATH=/tmp/chion-pluto:\$HOME/.julia`.
- If Julia is not on `PATH`, load it first with `module load julia/1.12.2`.
- If notebook cells feel slow, shorten the forcing vectors or lower `nx`, `ny`, or `ntot`.
- If you want GPU execution, switch to `02_gpu_single_device.jl` instead of rewriting the workflow.
"""

# ╔═╡ 0d8c12b5-885f-4e5d-96d1-b8f3f18944a3
md"""
## Adding A New Forcing Format

The extension point is `Chion.AbstractForcingSource` plus one `load_forcing` method:

```julia
struct CSVForcingSource <: Chion.AbstractForcingSource
    path::String
end

function Chion.load_forcing(source::CSVForcingSource; physics=Chion.SnowpackPhysicalConstants(), ntot=20)
    domain = Chion.SnowpackDomain(ncol=1, Ntot=ntot, c=physics)
    forcing = Chion.EquilibriumForcing(
        dt_days=[1.0],
        air_temperature=fill(physics.T0 - 10.0, 1, 1),
        snowfall_rate=zeros(1, 1),
        rainfall_rate=zeros(1, 1),
        shortwave_down=zeros(1, 1),
    )
    return Chion.SnowpackCaseData(domain, forcing; forcing_label=source.path)
end
```

After that, the rest of the workflow does not change:

```julia
case = Chion.build_case(name="csv_demo", forcing=CSVForcingSource("forcing.csv"), backend=:cpu)
result = Chion.run_case(case)
```
"""

# ╔═╡ Cell order:
# ╟─2c9d12d4-0c79-11ef-9f26-6b694f0c1001
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1002
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1003
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1004
# ╟─2c9d12d4-0c79-11ef-9f26-6b694f0c1005
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1006
# ╟─2c9d12d4-0c79-11ef-9f26-6b694f0c1007
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1008
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1009
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1010
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1011
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1012
# ╟─54d75be0-65f9-4e5c-b178-8a0d4272b677
# ╠═7012df88-b561-4962-bf36-d662f5d54bf3
# ╠═4c2846e2-7426-4b87-8e9f-c5b8bcb8aec7
# ╠═09351cdb-5748-4935-8f0a-e5c931fd716d
# ╟─2c9d12d4-0c79-11ef-9f26-6b694f0c1013
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1014
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1015
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1016
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1017
# ╠═6e837b58-a3bc-4d9d-b6f5-6d8c2b8d43d1
# ╠═0b0dd5b7-3044-4aae-b512-cf3d4e6d5dc9
# ╟─bbe5850d-eb9b-4c8f-a6c4-a87836f7b33d
# ╠═aaad1caa-0e49-4d1e-87d5-2ad2d3614acf
# ╟─2c9d12d4-0c79-11ef-9f26-6b694f0c1018
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1019
# ╟─2c9d12d4-0c79-11ef-9f26-6b694f0c1020
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1021
# ╟─2c9d12d4-0c79-11ef-9f26-6b694f0c1022
# ╟─2c9d12d4-0c79-11ef-9f26-6b694f0c1023
# ╟─0d8c12b5-885f-4e5d-96d1-b8f3f18944a3
