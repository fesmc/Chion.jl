#=
Numeric tolerances and default values used across Chion. The physical constants are fields of
`Params`, the scheme choices types in `schemes.jl`.
=#

# ---------------------------------------------------------------------------
# Numeric tolerances
# ---------------------------------------------------------------------------
const EPS_TINY = 1.0e-12
const EPS_EMPTY_LAYER = 1.0e-10

# ---------------------------------------------------------------------------
# Time conversion defaults
# ---------------------------------------------------------------------------
const DEFAULT_SECONDS_PER_DAY = 86_400.0

# ---------------------------------------------------------------------------
# Atmospheric defaults
# ---------------------------------------------------------------------------
const DEFAULT_SEA_LEVEL_AIR_PRESSURE = 101_325.0
const DEFAULT_GRAVITY = 9.80665
const DEFAULT_MOLAR_MASS_DRY_AIR = 0.0289644
const DEFAULT_UNIVERSAL_GAS_CONSTANT = 8.31446261815324

# ---------------------------------------------------------------------------
# Snowpack column defaults
# ---------------------------------------------------------------------------
const DEFAULT_NTOT = 15
const DEFAULT_MASS_MAX = 500.0
const DEFAULT_MASS_SPLIT = 300.0
const DEFAULT_MASS_MIN = 100.0
const DEFAULT_DENSITY_INIT = 300.0
const DEFAULT_TEMPERATURE_INIT = 273.0

# Legacy BESSI reference uses a 15-layer column when capping total column depth.
const BESSI_REFERENCE_LAYER_COUNT = 15
const BESSI_REFERENCE_DEPTH_DENSITY = 300.0
const BESSI_REFERENCE_SNOW_DEPTH_M = 22.5
