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

This notebook shows the standard Chion workflow on CPU:

1. activate the project
2. construct a `SnowpackGrid` and a `BESSIModel`
3. define a `SnowpackForcing`
4. build a `Simulation` and call `run!(simulation)`
5. inspect the result via `Chion.get_state(simulation.model.domain, idx)`
"""

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1004
PlutoUI.TableOfContents()

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1006
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

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1007
md"""
## Single-Column Example
"""

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1009
begin
	single_grid  = SnowpackGrid(CPU(), 1)
	single_model = seed_domain!(
		BESSIModel(single_grid;
			albedo=DynamicAlbedo(), densification=BESSIDensification(),
			fresh_snow_density=ConstantFreshSnowDensity(), Ntot=5);
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
	single_simulation = Simulation(single_model;
		forcing        = single_forcing,
		cycles         = 1,
		backend        = :cpu,
		save           = Symbol[],
		output_dir     = output_dir_for("01_cpu_workflows"),
		write_outputs  = true,
		history_stride = 1,
	)
end

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1010
single_result = run!(single_simulation)

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1012
let
	(
		status          = single_result.status,
		forcing_steps   = length(single_forcing.time_values),
		history_records = length(single_result.history),
		summary_path    = single_result.summary_path,
		column_summary  = summarize_column(single_model, 1),
	)
end

# ╔═╡ 54d75be0-65f9-4e5c-b178-8a0d4272b677
md"""
## Single-Column Plots
"""

# ╔═╡ 7012df88-b561-4962-bf36-d662f5d54bf3
single_forcing_plot = plots_available() ?
	forcing_timeseries_plot(single_forcing; idx=1, title_prefix="CPU single-column") :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ 4c2846e2-7426-4b87-8e9f-c5b8bcb8aec7
single_forcing_summary = (
	forcing_steps            = length(single_forcing_vectors.dt_days),
	dt_days_preview          = single_forcing_vectors.dt_days[1:5],
	air_temperature_c_range  = (minimum(single_forcing_vectors.air_temperature_c), maximum(single_forcing_vectors.air_temperature_c)),
)

# ╔═╡ 09351cdb-5748-4935-8f0a-e5c931fd716d
single_profile_plot = plots_available() ?
	column_profile_plot(single_model, 1; title="CPU single-column") :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1013
md"""
## Multi-Column Example
"""

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1014
multi_dims = (nx=40, ny=50)

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1015
begin
	multi_surface_mass = [
		220.0 + 35.0 * ((i - 1) / max(multi_dims.nx - 1, 1)) + 20.0 * ((j - 1) / max(multi_dims.ny - 1, 1))
		for j in 1:multi_dims.ny for i in 1:multi_dims.nx
	]
	multi_grid  = regular_grid(multi_dims.nx, multi_dims.ny)
	multi_model = seed_domain!(
		BESSIModel(multi_grid; Ntot=5);
		surface_mass=multi_surface_mass, density=320.0, temperature_c=-12.0,
	)
	multi_forcing = SnowpackForcing(
		dt_days          = single_forcing_vectors.dt_days,
		air_temperature_c = single_forcing_vectors.air_temperature_c,
		snowfall_mm_day  = single_forcing_vectors.snowfall_mm_day,
		rainfall_mm_day  = single_forcing_vectors.rainfall_mm_day,
		shortwave_down   = single_forcing_vectors.shortwave_down,
		wind_speed       = single_forcing_vectors.wind_speed,
		ncol             = multi_dims.nx * multi_dims.ny,
	)
	multi_simulation = Simulation(multi_model;
		forcing        = multi_forcing,
		cycles         = 1,
		backend        = :cpu,
		save           = Symbol[],
		output_dir     = output_dir_for("01_cpu_workflows"),
		write_outputs  = true,
		history_stride = 1,
	)
end

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1017
multi_result = run!(multi_simulation)

# ╔═╡ 0b0dd5b7-3044-4aae-b512-cf3d4e6d5dc9
let
	ncol = multi_dims.nx * multi_dims.ny
	sample_indices = unique([1, cld(ncol, 2), ncol])
	(
		status          = multi_result.status,
		forcing_steps   = length(multi_forcing.time_values),
		history_records = length(multi_result.history),
		sample_columns  = [summarize_column(multi_model, idx) for idx in sample_indices],
	)
end

# ╔═╡ bbe5850d-eb9b-4c8f-a6c4-a87836f7b33d
md"""
## Multi-Column Plot
"""

# ╔═╡ aaad1caa-0e49-4d1e-87d5-2ad2d3614acf
multi_thickness_plot = plots_available() ?
	layout_heatmap_plot(
		multi_grid,
		domain_metric_values(multi_model, :thickness);
		title="CPU multi-column final thickness",
		unit="m",
		color=:ice,
	) :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1018
md"""
## Physics-Scheme Variant Comparison
"""

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1019
scheme_comparison = let
	variants = [
		(label="default",      albedo=DynamicAlbedo(),  densification=BESSIDensification(),   fresh_snow=ConstantFreshSnowDensity()),
		(label="htessel_like", albedo=ConstantAlbedo(), densification=HTESSELDensification(), fresh_snow=ParameterizedFreshSnowDensity()),
	]
	[
		let
			g = SnowpackGrid(CPU(), 1)
			m = seed_domain!(
				BESSIModel(g; albedo=v.albedo, densification=v.densification,
					fresh_snow_density=v.fresh_snow, Ntot=5);
				surface_mass=250.0, density=320.0, temperature_c=-12.0,
			)
			f = SnowpackForcing(
				dt_days          = single_forcing_vectors.dt_days,
				air_temperature_c = single_forcing_vectors.air_temperature_c,
				snowfall_mm_day  = single_forcing_vectors.snowfall_mm_day,
				rainfall_mm_day  = single_forcing_vectors.rainfall_mm_day,
				shortwave_down   = single_forcing_vectors.shortwave_down,
				wind_speed       = single_forcing_vectors.wind_speed,
				ncol             = 1,
			)
			sim = Simulation(m; forcing=f, cycles=1, backend=:cpu,
				save=Symbol[], write_outputs=false, history_stride=1)
			r = run!(sim)
			(label=v.label, status=r.status, last_record=last(r.history),
			 column_summary=summarize_column(m, 1))
		end
		for v in variants
	]
end

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1020
md"""
## Output Inspection

- `result.history` contains one summary record per recorded cycle.
- `Chion.get_state(simulation.model.domain, idx)` returns a copied column snapshot.
- `simulation.model.domain` gives lower-level access to the internal state arrays.
- `result.summary_path` and `result.history_csv_path` point to saved text/CSV outputs.
- `Chion.get_state(simulation.model.domain, idx)` gives a named snapshot for a single column.
"""

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1022
md"""
## Performance Notes

- `backend=:cpu` dispatches to the threaded CPU backend.
- Keep `nx`, `ny`, forcing length, and `Ntot` small for interactive notebooks.
- Writing text outputs is cheap; enabling NetCDF adds extra I/O and requires
  a grid with spatial coordinates and a NetCDF library.
- For quick iteration, set `write_outputs=false` and shorten the forcing vectors.
"""

# ╔═╡ 2c9d12d4-0c79-11ef-9f26-6b694f0c1023
md"""
## Troubleshooting

- If `using Chion` fails in a restricted shell, set `JULIA_DEPOT_PATH=/tmp/chion-pluto:\$HOME/.julia`.
- If Julia is not on `PATH`, load it first with `module load julia/1.12.2`.
- If notebook cells feel slow, shorten the forcing vectors or lower `nx`, `ny`, or `Ntot`.
- For GPU execution, switch to `02_gpu_single_device.jl`.
"""

# ╔═╡ 0d8c12b5-885f-4e5d-96d1-b8f3f18944a3
md"""
## Minimal Example

```julia
using Chion

grid  = SnowpackGrid(CPU(), 1)
model = BESSIModel(grid; albedo=DynamicAlbedo(), Ntot=5)

forcing = SnowpackForcing(
    dt_days=[1.0, 1.0],
    air_temperature_c=[-12.0, -10.0],
    snowfall_mm_day=[0.3, 0.0],
    rainfall_mm_day=[0.0, 0.1],
    shortwave_down=[120.0, 180.0],
)

simulation = Simulation(model; forcing=forcing, cycles=1)
run!(simulation)
```
"""

# ╔═╡ Cell order:
# ╟─2c9d12d4-0c79-11ef-9f26-6b694f0c1001
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1002
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1003
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1004
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1006
# ╟─2c9d12d4-0c79-11ef-9f26-6b694f0c1007
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1009
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1010
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1012
# ╟─54d75be0-65f9-4e5c-b178-8a0d4272b677
# ╠═7012df88-b561-4962-bf36-d662f5d54bf3
# ╠═4c2846e2-7426-4b87-8e9f-c5b8bcb8aec7
# ╠═09351cdb-5748-4935-8f0a-e5c931fd716d
# ╟─2c9d12d4-0c79-11ef-9f26-6b694f0c1013
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1014
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1015
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1017
# ╠═0b0dd5b7-3044-4aae-b512-cf3d4e6d5dc9
# ╟─bbe5850d-eb9b-4c8f-a6c4-a87836f7b33d
# ╠═aaad1caa-0e49-4d1e-87d5-2ad2d3614acf
# ╟─2c9d12d4-0c79-11ef-9f26-6b694f0c1018
# ╠═2c9d12d4-0c79-11ef-9f26-6b694f0c1019
# ╟─2c9d12d4-0c79-11ef-9f26-6b694f0c1020
# ╟─2c9d12d4-0c79-11ef-9f26-6b694f0c1022
# ╟─2c9d12d4-0c79-11ef-9f26-6b694f0c1023
# ╟─0d8c12b5-885f-4e5d-96d1-b8f3f18944a3
