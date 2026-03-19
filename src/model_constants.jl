"""
Shared non-physical constants used across snowpack-process modules.
"""

# Numeric tolerances
const EPS_TINY = 1.0e-12
const EPS_EMPTY_LAYER = 1.0e-10

# Time conversion defaults
const DEFAULT_SECONDS_PER_DAY = 86_400.0
const DEFAULT_SECONDS_PER_MONTH = DEFAULT_SECONDS_PER_DAY * 30.0
const DEFAULT_SECONDS_PER_YEAR = DEFAULT_SECONDS_PER_MONTH * 12.0

# Snowpack column default setup
const DEFAULT_NTOT = 7
const DEFAULT_N_ACTIVE = 0
const DEFAULT_MASS_MAX = 500.0
const DEFAULT_MASS_SPLIT = 300.0
const DEFAULT_MASS_MIN = 100.0
const DEFAULT_RHO_MAX = 900.0
const DEFAULT_F_BASE_MAX = 0.6
const DEFAULT_DENSITY_INIT = 300.0
const DEFAULT_TEMPERATURE_INIT = 273.0

# Legacy BESSI reference uses a 15-layer column when capping total column mass.
const BESSI_REFERENCE_LAYER_COUNT = 15
