### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4002
begin
	import Pkg
	const NOTEBOOK_DIR = @__DIR__
	const PROJECT_ROOT = normpath(joinpath(NOTEBOOK_DIR, "..", ".."))
	Pkg.activate(PROJECT_ROOT)
	PROJECT_ROOT
end

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4003
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

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4001
md"""
# Chion External Forcing-File Scaffold

This notebook is optional. It reuses the same public case workflow as the synthetic notebooks, but the domain and forcing come from an external prepared HDF5/NetCDF forcing file.
"""

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4004
PlutoUI.TableOfContents()

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4005
md"""
## Prerequisites

- Load the cluster modules before launching Pluto:
  - `module load julia/1.12.2`
  - `module load hdf5`
  - `module load netcdf-c`
  - `module load cuda/13.1.0` only if you want GPU execution
- Point `forcing_file_path` below at a readable prepared forcing file.
- Keep `enable_netcdf=false` unless `NETCDF_LIB` is set and you want file output beyond the text summary/history CSV.
- The file-backed case is large. This notebook uses explicit `load_case_data` and `run_case_now` gates so Pluto does not automatically load and run the whole dataset on open.
"""

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4006
default_forcing_file_path = let candidates = [
	"/p/projects/ou/labs/ai/Nils/MARv3.14.3-10km-daily-ERA5-2025.nc",
	"/Users/niboch001/Downloads/MARv3.14.3-10km-daily-ERA5-2026.nc",
]
	found = ""
	for candidate in candidates
		if isfile(candidate)
			found = candidate
			break
		end
	end
	found
end

# ╔═╡ 1f51ba02-55a5-455d-aa59-1840d430d149
forcing_file_path = default_forcing_file_path

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4007
backend_choice = :gpu

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4008
enable_netcdf = false

# ╔═╡ 4ad0f4a4-c267-4551-a9f3-6e3996fcd65e
load_case_data = true

# ╔═╡ 8cfefc4a-aa0d-4211-ba99-d516f37603a4
run_case_now = true

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4009
physics = Chion.physics(
	albedo=:dynamic,
	densification=:bessi,
	fresh_snow_density=:constant,
)

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4010
forcing_file_status = if isempty(strip(forcing_file_path))
	(
		ready=false,
		message="Set `forcing_file_path` to a prepared forcing file path before running the loader cells.",
	)
elseif !isfile(strip(forcing_file_path))
	(
		ready=false,
		message="The configured `forcing_file_path` does not exist or is not readable: $(strip(forcing_file_path))",
	)
elseif enable_netcdf && isempty(strip(get(ENV, "NETCDF_LIB", "")))
	(
		ready=false,
		message="NetCDF output was requested, but `NETCDF_LIB` is not set. Disable NetCDF or configure the library path first.",
	)
elseif !load_case_data
	(
		ready=false,
		message="Path looks valid. Set `load_case_data=true` when you want to read the forcing file into memory.",
	)
else
	(
		ready=true,
		message="Forcing-file path and output settings look valid. The loader is enabled.",
	)
end

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4011
md"""
## Path Status And Loader Gate
"""

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4012
forcing_file_status

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4013
forcing_file_case = forcing_file_status.ready ? Chion.prescribed_case(
	forcing_file=strip(forcing_file_path),
	physics=physics,
	ntot=15,
	run=Chion.RunConfig(
		name="forcing_file_scaffold",
		backend=backend_choice,
		output_dir=output_dir_for("90_forcing_file_external_scaffold"),
		write_outputs=true,
		write_netcdf=enable_netcdf,
		cycles=100,
		history_stride=1,
	),
) : nothing

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4014
forcing_file_data = isnothing(forcing_file_case) ? nothing : forcing_file_case.definition

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4015
forcing_file_run = if isnothing(forcing_file_case)
	nothing
elseif !run_case_now
	nothing
else
	run_case_capture(forcing_file_case)
end

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4016
forcing_file_summary = isnothing(forcing_file_run) ? (
	message=isnothing(forcing_file_case) ? forcing_file_status.message : "Case is loaded. Set `run_case_now=true` to execute one equilibrium cycle.",
	case=isnothing(forcing_file_case) ? nothing : forcing_file_case,
	metadata=isnothing(forcing_file_data) ? nothing : forcing_file_data.metadata,
	loader_notes=isnothing(forcing_file_data) ? nothing : forcing_file_data.notes,
) : (
	case=forcing_file_case,
	metadata=forcing_file_data.metadata,
	loader_notes=forcing_file_data.notes,
	status=forcing_file_run.result.status,
	cycles_completed=length(forcing_file_run.result.history),
	summary_path=forcing_file_run.result.summary_path,
	history_csv_path=forcing_file_run.result.history_csv_path,
	nc_path=forcing_file_run.result.nc_path,
	column_1=summarize_column(forcing_file_run.result, 1),
)

# ╔═╡ c5ed28d4-e3e6-43bc-a48e-28fc42be111c
md"""
## Spatial Plots
"""

# ╔═╡ f9557bf3-6977-4979-b5f0-38ec3898bc46
forcing_file_initial_thickness_plot = isnothing(forcing_file_data) ? forcing_file_status.message :
	plots_available() ?
	layout_heatmap_plot(
		forcing_file_data.layout,
		domain_metric_values(forcing_file_data.domain, :thickness);
		title="Initial snow thickness",
		unit="m",
		color=:ice,
	) :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ 819ad742-b7da-45ab-b13c-d0afc1454078
forcing_file_initial_density_plot = isnothing(forcing_file_data) ? forcing_file_status.message :
	plots_available() ?
	layout_heatmap_plot(
		forcing_file_data.layout,
		domain_metric_values(forcing_file_data.domain, :bulk_density);
		title="Initial bulk density",
		unit="kg/m^3",
		color=:dense,
	) :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ 6bc7f31e-7124-4c80-8f77-bb312117f6db
md"""
## Run Plots
"""

# ╔═╡ 8cdb0f3c-4f3a-4108-99f1-49f54b985462
forcing_file_history_plot = isnothing(forcing_file_run) ? forcing_file_summary.message :
	plots_available() ?
	history_plot(forcing_file_run.result.history; title="Equilibrium cycle history") :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ d27e0310-0bfe-4889-8397-1d7891dc004e
forcing_file_final_thickness_plot = isnothing(forcing_file_run) ? forcing_file_summary.message :
	plots_available() ?
	layout_heatmap_plot(
		forcing_file_data.layout,
		domain_metric_values(forcing_file_run.result, :thickness);
		title="Final snow thickness",
		unit="m",
		color=:ice,
	) :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ 0f6e6f13-ac87-48cd-9dc4-b4e5f4f26ced
forcing_file_delta_thickness_plot = isnothing(forcing_file_run) ? forcing_file_summary.message :
	plots_available() ?
	layout_heatmap_plot(
		forcing_file_data.layout,
		domain_metric_values(forcing_file_run.result, :thickness) .- domain_metric_values(forcing_file_data.domain, :thickness);
		title="Thickness change after run",
		unit="m",
		symmetric=true,
	) :
	"Plots.jl is not available in this Pluto session."

# ╔═╡ a41f2b9a-0c79-11ef-9f26-f7ea663b4017
md"""
## Notes

- This notebook is intentionally not the primary onboarding path; use the synthetic notebooks first.
- `prescribed_case(; forcing_file=...)` expects a prepared forcing file; do any spatial masking outside Chion before loading.
- `load_case_data=true` reads the full external forcing into memory.
- `run_case_now=true` executes the model after the case is loaded.
- Keep `write_netcdf=false` unless you explicitly need NetCDF output and have a valid `NETCDF_LIB`.
- If you switch `backend_choice` to `:gpu`, make sure `CUDA.functional()` is true first.
"""

# ╔═╡ Cell order:
# ╟─a41f2b9a-0c79-11ef-9f26-f7ea663b4001
# ╠═a41f2b9a-0c79-11ef-9f26-f7ea663b4002
# ╠═a41f2b9a-0c79-11ef-9f26-f7ea663b4003
# ╠═a41f2b9a-0c79-11ef-9f26-f7ea663b4004
# ╟─a41f2b9a-0c79-11ef-9f26-f7ea663b4005
# ╠═a41f2b9a-0c79-11ef-9f26-f7ea663b4006
# ╠═1f51ba02-55a5-455d-aa59-1840d430d149
# ╠═a41f2b9a-0c79-11ef-9f26-f7ea663b4007
# ╠═a41f2b9a-0c79-11ef-9f26-f7ea663b4008
# ╠═4ad0f4a4-c267-4551-a9f3-6e3996fcd65e
# ╠═8cfefc4a-aa0d-4211-ba99-d516f37603a4
# ╠═a41f2b9a-0c79-11ef-9f26-f7ea663b4009
# ╠═a41f2b9a-0c79-11ef-9f26-f7ea663b4010
# ╟─a41f2b9a-0c79-11ef-9f26-f7ea663b4011
# ╠═a41f2b9a-0c79-11ef-9f26-f7ea663b4012
# ╠═a41f2b9a-0c79-11ef-9f26-f7ea663b4013
# ╠═a41f2b9a-0c79-11ef-9f26-f7ea663b4014
# ╠═a41f2b9a-0c79-11ef-9f26-f7ea663b4015
# ╠═a41f2b9a-0c79-11ef-9f26-f7ea663b4016
# ╟─c5ed28d4-e3e6-43bc-a48e-28fc42be111c
# ╠═f9557bf3-6977-4979-b5f0-38ec3898bc46
# ╠═819ad742-b7da-45ab-b13c-d0afc1454078
# ╟─6bc7f31e-7124-4c80-8f77-bb312117f6db
# ╠═8cdb0f3c-4f3a-4108-99f1-49f54b985462
# ╠═d27e0310-0bfe-4889-8397-1d7891dc004e
# ╠═0f6e6f13-ac87-48cd-9dc4-b4e5f4f26ced
# ╟─a41f2b9a-0c79-11ef-9f26-f7ea663b4017
