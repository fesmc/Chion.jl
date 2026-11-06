# Research scripts

Scripts behind parts of the paper that do not fit a configuration file. They use the low-level
`Simulation` API and read their settings from `CHION_*` environment variables. The model runs behind
the paper figures are in `run_paper_experiments.sh` (configurations in `examples/configs`, scripts in
`examples/runs`).

- `run_gris_fully_prescribed_transient.jl`: shared MAR helpers (`load_mar`, `read_mar_columns`,
  `read_mar_static_columns`, MAR budgets) used by several figure scripts.
- `run_gris_native_3hourly_transient.jl` (with `run_gris_1981_three_hourly_comparison.jl`): Greenland
  1980–2017 with the native 3-hourly MAR forcing (`run_paper_experiments.sh native`; Figs. E1–E3).
- `fit_cloud_proxy_longwave.jl`: fit of the cloud-proxy longwave parameterisation to MAR (Fig. L1).
- `performance_benchmark.jl`, `performance_greenland.jl`: CPU and GPU timings (Fig. P).
- Earlier studies, kept for reference: `calibrate_gris_surface_parameterization.jl`,
  `run_gris_1981_*`, `run_gris_native_3hourly_cycled_equilibrium.jl`,
  `run_gris_minimal_parameterized_comparison.jl`, `run_gris_prescribed_parameterized_comparison.jl`,
  `benchmark_snowpack_columns.jl`, `single_column_conservation_experiment.jl`.
