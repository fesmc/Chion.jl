# Chion.jl — Agent Rules

## Project Overview

Chion.jl is a Julia package for fast snowpack and firn mass and energy balance on CPUs and GPUs. It
computes the surface mass balance of ice sheets and snow-covered land from atmospheric forcing. The
main model, `BESSIModel`, is a layered column model: a surface energy balance (SEMIX-type turbulent
fluxes, parameterized longwave, Robin surface boundary) over snow and firn layers with accumulation,
melt, densification, percolation and refreezing, plus a thermal ice substrate. `PDDModel` and
`ITMModel` are bulk alternatives. Columns are independent, so one GPU thread or CPU task steps one
column. Chion is designed to be modular in the sense that (almost) all components should be
exchangeable with alternative implementations. Chion is also intended to become fully
differentiable.

## Language & Environment

- **Julia 1.10+** (`[compat]` admits 1.9 to 1.12; the manifests are resolved with 1.12.2) | CPU and GPU
- On the cluster Julia is not on `PATH`: `module load julia/1.12.2 hdf5 netcdf-c`. GPU jobs go through
  `sbatch` (see the `run_*_gpu.sh` scripts for the SLURM header and `JULIA_DEPOT_PATH`)
- **Key packages**: KernelAbstractions.jl, Adapt.jl, NCDatasets.jl, CUDA.jl (weak dependency, loaded
  by `using CUDA`), Enzyme.jl (forward mode on the column step), Reactant.jl (reverse-mode prototype,
  see Differentiability)

## Critical Rules

### Kernel Functions (GPU compatibility)

- Use `@kernel` / `@index` (KernelAbstractions.jl)
- Kernels must be **type-stable** and **allocation-free**. This includes everything they call: the
  whole column step (`column_step!` → `column_step_core!` → the `processes/` functions) is inlined into
  the GPU kernel
- Use `ifelse` — never short-circuiting `if`/`else` in kernels. Existing column code still has
  `if`/`return` branches (layer split/merge, bare-ice branch); new code uses `ifelse`, and converting
  old branches is a separate change, not a side effect of another one
- No error messages
- Mark functions called inside kernels with `@inline` (or `Base.@propagate_inbounds`)
- **Never loop over grid points outside kernels**: columns are stepped only through
  `_launch_step_columns!` (`src/step.jl`), which runs the KA kernel on the GPU and a precompiled
  threaded loop on the CPU. Loops over layers *inside* a column step are fine

### Type Stability & Memory

- All structs must be concretely typed
- Type annotations are for **dispatch**, not documentation
- Minimize allocation; favor inline computation
- **Never hardcode Float64**: no float literal such as `0.0` or `1.0` in kernels or in code that builds
  arrays. Use `zero(NF)`, `one(NF)`, `NF(0.5)`, `oftype(x, 0.5)`, `convert(NF, 1//2)`, or rational
  literals. `NF` is the model's number type (`Params(; NF)`, Float64 or Float32 today); write code so that it may
  also be a ForwardDiff `Dual` or a Reactant traced number. Literal defaults are fine in host-side
  constructors whose values are converted to `NF` afterwards (`BESSIModel` keywords,
  the `@kwdef` defaults of the schemes and `Layering`, `_param_defaults`, example configs)
- Absolute thresholds such as `EPS_TINY = 1e-12` are below one Float32 ulp of typical layer masses.
  A `while` loop or convergence test needs a relative floor (`_roundoff_floor`), or it can hang in
  Float32

### Parameterizations

- Each exchangeable scheme is a concrete subtype of an abstract scheme type
  (`AbstractAlbedoScheme`, `AbstractDensificationScheme`, `AbstractTurbulentFluxScheme`,
  `AbstractLongwaveScheme`, `AbstractDiurnalCycle`, `AbstractPDDMethod`, …; all in `src/schemes.jl`)
  and owns its coefficients as documented `@kwdef` fields (`AgingAlbedo(; cold_timescale_days = 20)`),
  as in Laddie. It is passed as an instance in a typed slot of `Params`
  (`Params(albedo_scheme = SemixAlbedo())`); constants shared by several schemes are flat `Params`
  fields. Names (`"semix"`) exist only in TOML and are converted once, by dispatch, in `config.jl`
  (`_SCHEME_TABLE`)
- A new scheme is a new type, its `validate` method if it has coefficients, its methods in
  `src/processes/`, an entry in `_SCHEME_TABLE`, and an export. Never a flag, a symbol comparison or an
  `isa` test in the model code (`test/code_quality.jl` checks this)
- Selecting a whole process (which column-step branch runs) dispatches on the scheme where the
  choice is made, e.g. `_column_step_diurnal!(::DiurnalSubsteps, …)`; a local operation dispatches
  inline (`_net_longwave(params.surface_energy_balance_scheme, …)`)
- Kernels receive arrays (`fields`, a NamedTuple from `get_fields(model.state)`), the forcing fields
  (`model.forcing.fields`), the workspace, and the isbits `Params` and `Layering`. Never a model,
  `Simulation`, `Forcing` or another container
- When a scheme changes the kernel arithmetic, prefer scalar coefficients computed at the launch site,
  and keep the existing default bit-identical
- Numbers read inside kernels, including those of scheme structs, go through `_val`, so that Reactant
  can trace them. Chion has no `_val` yet: introduce it with the first Reactant path in `src/`

### Reactant

- Reactant-specific code lives in `ext/ChionReactantExt.jl` (not yet created; the current prototype is
  the standalone port in `examples/research/differentiable_column/` on branch `fixed-grid-ad`, with
  its own `Project.toml`). The extension will need both Reactant and CUDA loaded (also on CPU)
- **No reachable throw path in kernels**: no `@assert`, error, un-elided bounds check, `InexactError`
  from `convert`/`round(Int, …)` or integer division by zero. These lower to traps that fail
  compilation even if never hit. Validate arguments in host-side constructors instead (as
  `BESSIModel` does)
- `ifelse` evaluates both branches: select the *index*, not the loads (`x[ifelse(c, i, n)]`)
- Functions traced by Reactant outside kernels must be branch-free on values: use `ifelse`, `min`,
  `max`, `clamp` instead of `if`, `? :`, `&&`, `||`
- Closures compiled into kernels capture only `isbits` values (no `Type`s such as `NF`)
- Do not run Reactant tests with `--check-bounds=yes`
- Known limits: the step kernel takes `time_index::Int`, so a `@trace` loop over time fails;
  `to_rarray` cannot rebuild the compact forcing (needs dense forcing); Reactant lacks `exp10` and
  `cbrt`; loop-carried `@trace for` values are written back into input buffers (copy the state first)

### Differentiability

- Supported AD paths (`test/test_autodiff.jl`, checked against finite differences): Enzyme forward
  mode through `column_step!` on the CPU and inside the KA kernel on the GPU (`autodiff_deferred`),
  and ForwardDiff through a whole model built at `NF = Dual` (`Params(; NF, alpha_dry = Dual(…))`,
  `integrate!` or a `Simulation`; `NF <: Real`). Both must keep working. Enzyme reverse mode through
  Reactant is the target; it exists only in the fixed-grid port (`fixed-grid-ad`). Do not add code that
  blocks any of these paths
- Host code that leaves the number type (outputs, logs, metadata) converts with `_primal` (the value of
  a `Dual`, from `ChionForwardDiffExt`), never with `Float64(x)` on model numbers
- Every field of a scheme, boundary condition or `Layering` declared with the type parameter `NF` is
  converted to the model's number type by `_promote_param` (from the struct definition); declare
  integer fields (counts) as `Int`
- Never call `Enzyme.API.strictAliasing!(false)`: it silently gave wrong tangents and broke reverse mode
- No global or hidden mutable state; every input of a computation is passed explicitly. The exceptions
  are launch configuration that never changes results: `CHION_CPU_TIME_BLOCK_STEPS` /
  `CHION_GPU_TIME_BLOCK_STEPS` (time steps per launch; columns are independent) and the cached
  multiprocessor count in `ext/ChionCUDAExt.jl` (block size)
- Guard non-differentiable points instead of introducing new kinks. Powers: double-where
  `ifelse(x > 0, ifelse(x > 0, x, one(x))^p, zero(x))` (not `max(x, floatmin)^p`, whose derivative
  overflows). Divisions: `ifelse(b > ε, a / max(b, ε), zero(a))` (the double-where form of a division
  is folded away by LLVM/XLA and gives NaN adjoints). Add shared helpers (`_safe_pow`, `_safe_div`,
  `_safe_sqrt`) rather than repeating the pattern. A new `abs`/`max`/`min` in the dynamics needs a
  comment on why the kink is acceptable
- Layer split/merge, snow↔ice switches and remeshing make annual SMB piecewise smooth (jumps of
  0.2–7 kg m⁻² per 0.01 K): AD gives the pathwise slope. Do not add new discrete switches to the
  column step where a continuous formulation exists
- A change to the dynamics or to a scheme needs a gradient check: extend `test/test_autodiff.jl`
  (Enzyme forward and ForwardDiff against finite differences on a short column run), and the FD tests in `examples/research/differentiable_column/run_tests.jl` when the port is
  touched

### Physics

- Every process cites its governing equation (and its source) in the docstring or in comments
- Mass and enthalpy are not conserved per layer: snowfall, rain, melt, runoff, sublimation and the
  surface and basal heat fluxes exchange them with the atmosphere and the ice below. Test budgets
  instead: the change of a quantity equals its fluxes. Mass closes through `smb_ice` (bare-ice melt is
  subtracted there, not from `mass_base`); enthalpy closes through `column_enthalpy` and
  `workspace.boundary_enthalpy` over a full step. A process that does conserve a quantity (remeshing,
  split/merge conserve mass and enthalpy) has a test that verifies it
- Prefer continuous-time formulations (tendencies, relaxations with a time scale). Instantaneous
  adjustments (layer split/merge/remesh, snow-depth cap, bucket percolation, removing trace surface
  layers, clamps to the melting point) are exceptions: document why they are needed, and do not copy
  them as a pattern for new schemes

### Docstrings

- Use DocStringExtensions.jl: `$(TYPEDEF)` and `$(TYPEDFIELDS)` for structs (every field gets a
  docstring with unit and default, as in `Params`), `$(TYPEDSIGNATURES)` when the signature does not
  include default values for args and/or kwargs
- Use `jldoctest` blocks. Plain `julia` blocks are only exceptionally allowed since they are not tested
  by CI
- Include `# output` with verifiable output; prefer `show` methods over boolean comparisons
- Use unicode for math (`Δt`, `ρ`, `T₀`), not LaTeX — LaTeX doesn't render in the REPL

## Naming Conventions

- **Files**: snake_case matching the type or process they define (`processes/albedo.jl`)
- **Types/Constructors**: PascalCase **only for true constructors**
- **Functions**: snake_case — `column_step!`, `yearly_step!`; functions that return values are not
  PascalCase, except named setups that return a configured object like a constructor
- **Kernels**: underscore prefix and `_kernel!` suffix — `_step_columns_kernel!`
- **Variables**: English long name or readable unicode math notation — never mix abbreviated and
  full forms (e.g., `cond` vs `condition`) to imply a difference; be specific. Older names
  (`mass_w`, `Tsrf`, `Ntot`, `c` for the constants) stay until a dedicated rename; new names follow
  this rule, with units in the name where ambiguous (`dt_days`, `_mm_day`, `_c` for °C)

## Module Structure

`src/` holds the infrastructure (schemes, params, grid, forcing, state, model, step, outputs,
simulation, TOML config); `src/processes/` holds one file per process (`albedo.jl`,
`densification.jl`, `energy_flux.jl`, `percolation.jl`, …) with the kernel-side methods of its
schemes. The API follows Laddie.jl: `Grid` → model (`BESSIModel`, `PDDModel`, `ITMModel`) →
`Simulation`. The object structure of a run (abbreviated):

```
Simulation
├── model :: BESSIModel                 # mutable; knows nothing about time beyond the loaded year
│   ├── grid :: Grid{NF}                # ncol, optional SpatialLayout (x, y, is, js, mask),
│   │   └── layering :: Layering{NF}    # Ntot, mass_max/split/min, near-surface thicknesses, substrate
│   ├── state :: BESSIState             # prognostic arrays only:
│   │                                   #   layers (layer, column), TransposedLayerMatrix on a GPU:
│   │                                   #     mass, mass_w, density, temperature, ice_temperature
│   │                                   #   per column: N, mass_base, smb_ice, runoff, melt, refreezing,
│   │                                   #     vapor_mass, sublimation, *_sum SEB accumulators, Tsrf,
│   │                                   #     albedo, snow_age_days, w_snow_max
│   ├── cache :: BESSICache             # workspace (heat solve, boundary_enthalpy), diagnostics
│   │                                   # (thickness, wet_mass, bulk_density, liquid_water), active columns
│   ├── params :: Params{NF, schemes…}  # physical constants + typed scheme slots, each scheme
│   │                                   # carrying its coefficients (isbits, passed to kernels)
│   └── forcing :: AtmosphericForcing   # source::Forcing (files / in memory / function), loaded year:
│                                       #   data (host SnowpackForcing), fields (kernel view on device)
├── clock :: Clock                      # year_index, step, iteration, time
├── years, spinup                       # NoSpinup / CycledSpinup / ClimatologySpinup
├── output :: OutputConfig              # run directory (log.txt, chion_run.toml, restarts) and the
│                                       # field outputs (one NetCDF or CSV file each, calendar periods)
└── io :: IOState                       # open recorders, spun_up, wall time
```

Flat forwarding: `model.mass` (state), `model.thickness` (cache), `model.alpha_dry` (params),
`model.Ntot` (grid layering), `model.year`; names must be unique across them (checked at
construction). Print any object's structure:

```julia
using Chion, Dates
function fieldtree(x, depth = 0; maxdepth = 2)
    for f in fieldnames(typeof(x))
        v = getfield(x, f)
        println("  "^depth, f, " :: ", nameof(typeof(v)))
        depth < maxdepth - 1 && isstructtype(typeof(v)) && !(v isa AbstractArray) &&
            fieldtree(v, depth + 1; maxdepth)
    end
end
forcing = SnowpackForcing(dt_days = [1.0, 1.0], ncol = 4, air_temperature_c = [-15.0, -14.0],
    snowfall_mm_day = 1.0, rainfall_mm_day = 0.0, shortwave_down = [120.0, 140.0],
    wind_speed = 3.0, latitude_deg = 70.0,
    time_values = [DateTime(2001, 1, 1, 12), DateTime(2001, 1, 2, 12)])
fieldtree(Simulation(BESSIModel(Grid(4); forcing)); maxdepth = 3)
```

List the implementations of a scheme with `using InteractiveUtils; subtypes(Chion.AbstractAlbedoScheme)`.

## Typical Usage

```julia
using Chion
forcing = Forcing("/path/MARv3.14.3-10km-daily-ERA5-{year}.nc"; surface_fluxes = :parameterized)
grid = Grid(forcing; layering = Layering(Ntot = 15))
model = BESSIModel(grid; forcing, params = Params())   # ParamsAISv1() for Antarctica; backend = CUDABackend()
sim = Simulation(model; years = 1980:1989, spinup = CycledSpinup(1960:1979; cycles = 2),
    output = OutputConfig(dir = "output/x", restart_every = 10,
        outputs = [Output("monthly.nc"; variables = :monthly, frequency = :monthly)]))
run!(sim)                  # run!(sim; years = 1) continues year by year; time_step!(sim) one step
```

The same as a TOML file: `julia --project -m Chion examples/configs/greenland_parameterized.toml`
(`[grid]` = `Layering`, `[model]` = `Params` fields with schemes by name or as tables with `type`).
In-memory forcing: `BESSIModel(Grid(n); forcing = SnowpackForcing(…))`. Coupling: `yearly_step!(sim)`,
`set_active_mask!`, `set_input!`. Without a `Simulation` (AD): `integrate!(model, steps)`. For a GPU,
`using CUDA` and `backend = CUDABackend()` (or `to_backend(model, CUDABackend())`); `cpu_state(model)`
copies the state back.

## Time-Step Flow

`run!`/`time_step!` (`src/simulation.jl`) spin up once, then advance the clock year by year:
`update_forcing!` loads the year, `_advance!` steps up to the next end of an output period and records
the outputs. `_step!(model, …)` (`src/model.jl`) dispatches on the model type; for BESSI
`_step_columns!` (`src/step.jl`) splits the range into blocks and calls `_launch_step_columns!`,
which steps every active column over the block:

```
column_step!(fields, params, layering, idx, forcing, workspace)
├── solar geometry once per step; parameterized longwave (cloud proxy)
└── _column_step_diurnal!               NoDiurnalCycle: one call; DiurnalSubsteps: up to
    └── column_step_core!               max_substeps, each with its own SW / air-temperature average
        ├── accumulation                snowfall into the surface layer, trace layers removed,
        │                               near-surface remesh
        ├── bare ice?  → _step_bare_ice_substrate!  (SEB on ice, melt to runoff) and return
        ├── albedo                      scheme update (dynamic/aging/semix/prescribed/constant)
        ├── densification
        ├── energy flux                 implicit heat solve (Thomas), Robin surface boundary,
        │                               ice substrate; SEB accumulators
        ├── vapour flux                 sublimation/deposition on the surface layer
        ├── melt                        surface layers, then the ice below
        ├── percolation                 bucket scheme, runoff
        ├── liquid-water compaction
        ├── refreezing
        └── near-surface remesh         restores the layer geometry for the next solve
```

- Columns are independent: the step order across columns and the time blocking never change results.
  A refactor of the step must be bitwise identical on the CPU (`test/regression/bitwise.jl`). On the
  GPU, bitwise identity holds only while the control flow is unchanged: `ptxas` fuses multiply-adds
  into FMAs per basic block, and GPU bits also change with the CUDA.jl/LLVM version. A refactor that
  removes branches (e.g. a flag test replaced by dispatch) must instead keep the GPU differences at
  roundoff level, below the difference between the unchanged code's own GPU and CPU results
  (compare baseline and change with the same CUDA environment on the same GPU model)
- `ColumnarStepWorkspace` (`model.cache.workspace`) holds per-column scratch (`energy` heat-solve
  arrays, `liquid_water_before_energy`), overwritten by every step: do not read it across steps. The
  exception is `boundary_enthalpy`, a cumulative diagnostic for the enthalpy budget
- A new process usually enters as a function in `src/processes/` called from `column_step_core!` at
  the place its inputs are ready, and as a term in the energy or mass budget it changes

## Common Pitfalls

1. **Type instability** in kernels ruins GPU performance (and spills registers: the step kernel is
   already at ~255 registers per thread)
2. **Overconstraining types**: use annotations for dispatch, not documentation
3. **Missing imports**: tests will catch this — add to `using` statements
4. **Subtle bugs from missing method imports**, especially in extensions
5. **Expecting unexported names**: consider exporting them rather than changing user scripts
6. **"Type is not callable" errors**: variable name shadows a function — rename or qualify
7. **Quick fixes that break correctness**: if a test fails after a change, revisit the original edit
8. **Commented-out code**: delete it. Git is the journal — don't leave commented code, debugging
   artifacts, or stale copy-paste remnants
9. **Hardcoded Float64**: never use `0.0`, `1.0` in kernels, use `zero(NF)`, `NF(x)` or `oftype(x, y)`
10. **Changing dependencies**: never add, remove or change `[deps]`/`[weakdeps]` in `Project.toml`
    unless the task requires it; touch `[compat]` only when asked
11. **Scope creep**: keep a change to a single concern; unrelated cleanup goes in a separate commit or PR
12. **Closures over mutable state**: they break AD — pass the values as explicit arguments
13. **GPU scripts without `using CUDA`**: `CUDABackend()` and `backend = :gpu` need CUDA loaded
14. **Concurrent compilation on CPU threads**: the CPU path precompiles the column step before
    `Threads.@threads`; keep that when changing `_cpu_step_columns!` (otherwise every process call
    allocates)

## Design Principles

- **Dispatch over conditionals**: use Julia's type system and multiple dispatch instead of
  `if`/`else` branching. Backend-specific code goes in `ext/` extensions, not `if` branches in `src/`
- **Use `to_backend(model, backend)` (or `backend = …` at construction) and `cpu_state(model)` to move
  data between devices** — never manual `CuArray()` calls. Copying device data back to the host for
  I/O, output or restarts with `Array(x)` is fine
- **Defaults serve the common case**: avoid `nothing` defaults when a concrete default (like
  `backend = :threads`) covers 80% of usage. Minimize boilerplate for the typical user.
- **Keyword argument names must be consistent** across related types and constructors
- **Always use explicit `return`** in functions longer than one expression
- **One operation per line** as default; break long expressions across lines

## Testing

- Full suite: `julia --project=. -e 'using Pkg; Pkg.test()'`. `test/test_gpu.jl` runs only when a CUDA
  device is functional (submit through SLURM for GPU coverage)
- Single area, from the package root (the test environment is `test/Project.toml`, with Aqua and
  CUDA; `test/setup.jl` loads the shared packages):
  ```
  julia --project=test -e 'include("test/setup.jl"); include("test/test_simulation.jl")'
  ```
- Regression gate of refactors: `julia --project=. --threads 16 test/regression/bitwise.jl compare
  output/regression/base cpu` (27 cases incl. two GrIS years with output files; `save` writes a
  baseline)
- Julia stack traces are long: write test output to a file and `grep` it
- A failing test is a signal about the change, not about the test. Do not loosen a tolerance to make it
  pass
- **Stale manifests**: `Manifest.toml` is resolved per Julia version. If a *stdlib* fails to
  precompile, compare `grep -m1 julia_version Manifest.toml` with `julia --version`, and re-resolve
  that environment (`julia --project=<env> -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'`) instead
  of debugging the code
- `examples/research/differentiable_column/` (branch `fixed-grid-ad`) has its own environment and
  tests: `run_tests.jl` (FD-checked gradients), `regression.jl` (port vs its committed version)

## Documentation

- Build: `julia --project=docs docs/make.jl`. Build the docs after changing public API, docstrings or
  `docs/src/`: a stale `@docs` entry or a missing docstring of an export fails the build
  (`checkdocs = :exports`)
- `docs/Manifest.toml` must point to the checkout being documented (`path` under `[[deps.Chion]]`; it
  is the absolute path of the `dev_nils` checkout, so in the `Chion.jl-fixedgrid` worktree it
  documents `dev_nils`). Fix with
  `julia --project=docs -e 'using Pkg; Pkg.develop(path = "."); Pkg.instantiate()'`
- Cite with DocumenterCitations (`[Key](@cite)`) and add the BibTeX entry to
  `docs/src/assets/chion.bib`. Neither exists yet: set them up (docs environment only) with the first
  citation instead of writing references as plain text
- Write units in parentheses, (kg m⁻² yr⁻¹), not brackets, which clash with Markdown links

## Benchmarks

- Benchmarks live in `examples/research/` (`performance_benchmark*.jl`, `benchmark_snowpack_columns.jl`,
  `gpu_profile.jl`), launched by `run_performance_benchmark_{cpu,gpu}.sh`. Run them only when asked or
  after a major change to the compute path (column step, kernel launch, state layout, number type)
- Compare against the most recent stored results (`output/performance/*.csv`): deviations within ±20%
  are noise; report larger regressions

## Calibrated defaults are a regression contract

The defaults of `Params` (= `ParamsGrISv1()`) and `Layering` are the Greenland calibration against MAR
v3.14.3 (daily forcing); the Antarctic runs use `ParamsAISv1()` (`setup = "AISv1"` in TOML). The paper runs
(`run_paper_experiments.sh`) and the manuscript depend on them. `PDDModel` and `ITMModel` are pinned to
the Fortran budgets by tests in `test/test_case_api.jl`.

- Never change what the defaults produce as a side effect of another change: a refactor is bitwise
  identical on the CPU, in Float64 and Float32, and within the GPU roundoff envelope on the GPU (see
  Time-Step Flow)
- A new scheme or option keeps the existing default; a config or experiment pins any setting it relies
  on explicitly
- Changing a model-wide default requires checking the GrIS SMB, melt and runoff against MAR (and the
  Antarctic setup when it shares the parameter) and saying so in the plan. If a default changes on
  purpose, update the calibration notes and the affected test values in the same commit and say why

## Agent Behavior

- Prioritize type stability and GPU compatibility
- Maintain AD-compatibility (Enzyme forward now, Reactant reverse as the target)
- Follow established patterns in existing code
- Add tests for new functionality; update exports when adding public API
- Reference physics equations in comments when implementing dynamics

## Implementation plans

Major features, bug fixes with wide reach, and refactors that change substantial existing code need an
**implementation plan** that the human has **explicitly approved** before any code is written. A change
of direction (new bucket, different design) needs approval again.

Plans live in the internal `chion-roadmap/`, one Markdown file per topic (see below). Skeleton:

```md
# Title

Status YYYY-MM-DD: **planned / approved / in progress / done**. One sentence on where it stands.
Branch, base revision, code, tests, docs: pointers only.

## Originating prompt
## Rationale            physics/problem, design decisions, not in scope
## Tasks                buckets, each small enough to review and test on its own; approval per bucket
## Results              dated findings, numbers copied in (runs, jobs, figures)
## Open questions
```

- Update the status line whenever the state changes
- Record a change of direction as a new dated entry rather than rewriting history
- When done, move the plan to `chion-roadmap/archived/` and write the devnote

## Roadmaps and devnotes

- `chion-roadmap/` is **internal and not committed** (listed in `.git/info/exclude`): implementation
  plans and working notes, written while a change is developed. Do not cite it from code, docs or
  commits
- `devnotes/` holds the settled findings, summarised from the roadmaps: `improvements/` (changes in
  main), `deadends/` (ideas dropped, with the reason), `open-questions/` (unresolved). Read the
  relevant devnote before reworking an area, and do not retry a dead end without addressing why it
  failed
- Write devnotes following `devnotes/rules.md` (create it with the first devnote: one file per finding,
  dated, with the numbers, the run or job that produced them, and pointers to code and commits)
