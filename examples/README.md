# Examples

The configurations and scripts behind the figures of the Chion model description paper. All model
runs are in `run_paper_experiments.sh` in the repository root.

| Path | What it does |
| --- | --- |
| `configs/greenland_parameterized.toml` | Greenland 1940–2026, surface energy balance and albedo from Chion |
| `configs/greenland_prescribed.toml` | the same with the MAR surface energy fluxes and albedo prescribed |
| `configs/greenland_firn.toml` | firn density and temperature profiles every six months (SUMup comparison) |
| `configs/antarctica_parameterized.toml`, `antarctica_prescribed.toml` | Antarctica 1979–2026 (Antarctic surface calibration) |
| `configs/antarctica_firn.toml` | Antarctic firn profiles |
| `runs/nh_era5_snow.jl` | Northern Hemisphere seasonal snow from daily ERA5 |
| `runs/snowmip_sites.jl` | all ten ESM-SnowMIP sites |
| `figures/` | the paper figures; `compare_with_mar.jl <output dir>` compares a MAR-forced run with MAR |
| `research/` | native 3-hourly forcing, longwave fit, benchmarks; see `research/README.md` |

Run a configuration with (on the PIK cluster: `sbatch run_config_gpu.sh <config.toml>`)

```bash
julia --project=. -m Chion examples/configs/greenland_parameterized.toml --backend=gpu
```

`--backend=threads` runs on the CPU (start Julia with `--threads N`), and `--output-dir=DIR`
changes where the results go. The data paths in the files point to the PIK cluster; change
`[forcing] path` for other locations.

## Paper figures

| Figure | Script | Inputs |
| --- | --- | --- |
| Model schematic | `figures/model_schematic.py` | none |
| Verification (step response, periodic forcing, convergence) | `figures/verification_suite.jl` | none |
| Diurnal shortwave reconstruction | `figures/plot_diurnal_shortwave_cycle.jl` | none |
| L1 cloud-proxy longwave | `research/fit_cloud_proxy_longwave.jl` | MAR daily 1980–1989 and 2010–2019 |
| E1–E3 Greenland against MAR | `figures/evaluation_figures.jl` | groups `greenland` and `native` |
| A1, A3 Antarctica against MAR | `figures/antarctica_figures.jl` | group `antarctica` |
| D1–D3 firn against SUMup | `figures/firn_evaluation_figures.jl` (`CHION_DOMAIN=antarctica` for the AIS) | groups `greenland`, `antarctica` |
| S2, S3 ESM-SnowMIP | `figures/snowmip_figures.jl` | group `snow` |
| N1, N2 Northern Hemisphere snow | `figures/nh_snow_figures.jl` | group `snow` |
| P performance | `figures/performance_figures.jl` | `research/performance_benchmark.jl`, `research/performance_greenland.jl` |
