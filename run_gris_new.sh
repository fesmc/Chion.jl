#!/bin/bash

#SBATCH --qos=gpushort
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
##SBATCH --time=2-1:00:00
#SBATCH --time=0-2:30:00
#SBATCH --job-name=gris_mar_case
#SBATCH --gres=gpu:1

# Don't change anything below this line:
#SBATCH --output=/p/projects/ou/labs/ai/Nils/Chion.jl/logs/gris-mar-%j.log
#SBATCH --error=/p/projects/ou/labs/ai/Nils/Chion.jl/logs/gris-mar-%j.err

set -euo pipefail

module load julia
module load hdf5
module load netcdf-c

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export JULIA_EXCLUSIVE=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export H5DUMP_BIN="$(which h5dump)"
export H5LS_BIN="$(which h5ls)"
export NETCDF_LIB="${NETCDF_LIB:-libnetcdf.so}"

PROJECT_DIR=/p/projects/ou/labs/ai/Nils/Chion.jl
SCRIPT_PATH="${PROJECT_DIR}/examples/scripts/run_gris_forcing_file_case.jl"
LOG_DIR="${PROJECT_DIR}/logs"
DEFAULT_FORCING_PATH="/p/projects/ou/labs/ai/Nils/MARv3.14.3-10km-daily-ERA5-2025.nc"

# Defaults can be overridden with exported env vars or script args.
FORCING_PATH="${FORCING_PATH:-${DEFAULT_FORCING_PATH}}"
OUTPUT_DIR="${OUTPUT_DIR:-}"
NETCDF_PATH="${NETCDF_PATH:-}"
BACKEND="${BACKEND:-}"
YEARS="${YEARS:-}"
HISTORY_YEAR_STRIDE="${HISTORY_YEAR_STRIDE:-}"
THREAD_COUNT="${THREAD_COUNT:-8}"
WRITE_OUTPUTS="${WRITE_OUTPUTS:-}"
WRITE_NETCDF="${WRITE_NETCDF:-}"
NETCDF_VARIABLES="${NETCDF_VARIABLES:-}"
MODEL="${MODEL:-}"
MASK_THRESHOLD="${MASK_THRESHOLD:-}"
NTOT="${NTOT:-}"
PDD_DDF_SNOW="${PDD_DDF_SNOW:-}"
PDD_DDF_ICE="${PDD_DDF_ICE:-}"
PDD_REFREEZING_FRACTION="${PDD_REFREEZING_FRACTION:-}"
TURBULENT_FLUX_SIGN="${TURBULENT_FLUX_SIGN:-}"
ALBEDO="${ALBEDO:-}"
DENSIFICATION="${DENSIFICATION:-}"
FRESH_SNOW_DENSITY="${FRESH_SNOW_DENSITY:-}"
EXTRA_ARGS="${EXTRA_ARGS:-}"

mkdir -p "${LOG_DIR}"
cd "${PROJECT_DIR}"

has_arg_prefix() {
    local prefix="$1"
    local arg
    for arg in "${EXTRA_ARGS_ARR[@]}"; do
        if [[ "${arg}" == "${prefix}"* ]]; then
            return 0
        fi
    done
    return 1
}

append_arg() {
    local arg="$1"
    EXTRA_ARGS_ARR+=("${arg}")
}

EXTRA_ARGS_ARR=()
if [[ $# -gt 0 ]]; then
    EXTRA_ARGS_ARR=("$@")
elif [[ -n "${EXTRA_ARGS}" ]]; then
    read -r -a EXTRA_ARGS_ARR <<< "${EXTRA_ARGS}"
fi

if [[ -n "${FORCING_PATH}" ]] && ! has_arg_prefix "--forcing-file="; then
    append_arg "--forcing-file=${FORCING_PATH}"
fi
if [[ -n "${OUTPUT_DIR}" ]] && ! has_arg_prefix "--output-dir="; then
    append_arg "--output-dir=${OUTPUT_DIR}"
fi
if [[ -n "${NETCDF_PATH}" ]] && ! has_arg_prefix "--netcdf-path="; then
    append_arg "--netcdf-path=${NETCDF_PATH}"
fi
if [[ -n "${BACKEND}" ]] && ! has_arg_prefix "--backend="; then
    append_arg "--backend=${BACKEND}"
fi
if [[ -n "${YEARS}" ]] && ! has_arg_prefix "--years="; then
    append_arg "--years=${YEARS}"
fi
if [[ -n "${NETCDF_VARIABLES}" ]] && ! has_arg_prefix "--netcdf-vars="; then
    append_arg "--netcdf-vars=${NETCDF_VARIABLES}"
fi
if [[ -n "${MODEL}" ]] && ! has_arg_prefix "--model="; then
    append_arg "--model=${MODEL}"
fi
if [[ -n "${MASK_THRESHOLD}" ]] && ! has_arg_prefix "--mask-threshold="; then
    append_arg "--mask-threshold=${MASK_THRESHOLD}"
fi
if [[ -n "${NTOT}" ]] && ! has_arg_prefix "--ntot="; then
    append_arg "--ntot=${NTOT}"
fi
if [[ -n "${PDD_DDF_SNOW}" ]] && ! has_arg_prefix "--pdd-ddf-snow="; then
    append_arg "--pdd-ddf-snow=${PDD_DDF_SNOW}"
fi
if [[ -n "${PDD_DDF_ICE}" ]] && ! has_arg_prefix "--pdd-ddf-ice="; then
    append_arg "--pdd-ddf-ice=${PDD_DDF_ICE}"
fi
if [[ -n "${PDD_REFREEZING_FRACTION}" ]] && ! has_arg_prefix "--pdd-refreezing-fraction="; then
    append_arg "--pdd-refreezing-fraction=${PDD_REFREEZING_FRACTION}"
fi
if [[ -n "${ALBEDO}" ]] && ! has_arg_prefix "--albedo="; then
    append_arg "--albedo=${ALBEDO}"
fi
if [[ -n "${DENSIFICATION}" ]] && ! has_arg_prefix "--densification="; then
    append_arg "--densification=${DENSIFICATION}"
fi
if [[ -n "${FRESH_SNOW_DENSITY}" ]] && ! has_arg_prefix "--fresh-snow-density="; then
    append_arg "--fresh-snow-density=${FRESH_SNOW_DENSITY}"
fi

case "${WRITE_OUTPUTS,,}" in
    false|0|no|off)
        has_arg_prefix "--no-output" || append_arg "--no-output"
        ;;
esac
case "${WRITE_NETCDF,,}" in
    false|0|no|off)
        has_arg_prefix "--no-nc" || append_arg "--no-nc"
        ;;
esac

REQUESTED_BACKEND="gpu"
for arg in "${EXTRA_ARGS_ARR[@]}"; do
    if [[ "${arg}" == --backend=* ]]; then
        REQUESTED_BACKEND="${arg#--backend=}"
        break
    fi
done

SRUN_ARGS=()
if [[ "${REQUESTED_BACKEND}" == "threads" || "${REQUESTED_BACKEND}" == "cpu" ]]; then
    SRUN_ARGS+=(--cpu-bind=cores --distribution=block:block --mem-bind=local)
fi

srun "${SRUN_ARGS[@]}" env \
    JULIA_NUM_THREADS="${THREAD_COUNT}" \
    FORCING_PATH="${FORCING_PATH}" \
    OUTPUT_DIR="${OUTPUT_DIR}" \
    NETCDF_PATH="${NETCDF_PATH}" \
    BACKEND="${BACKEND}" \
    YEARS="${YEARS}" \
    HISTORY_YEAR_STRIDE="${HISTORY_YEAR_STRIDE}" \
    WRITE_OUTPUTS="${WRITE_OUTPUTS}" \
    WRITE_NETCDF="${WRITE_NETCDF}" \
    NETCDF_VARIABLES="${NETCDF_VARIABLES}" \
    MODEL="${MODEL}" \
    MASK_THRESHOLD="${MASK_THRESHOLD}" \
    NTOT="${NTOT}" \
    PDD_DDF_SNOW="${PDD_DDF_SNOW}" \
    PDD_DDF_ICE="${PDD_DDF_ICE}" \
    PDD_REFREEZING_FRACTION="${PDD_REFREEZING_FRACTION}" \
    TURBULENT_FLUX_SIGN="${TURBULENT_FLUX_SIGN}" \
    ALBEDO="${ALBEDO}" \
    DENSIFICATION="${DENSIFICATION}" \
    FRESH_SNOW_DENSITY="${FRESH_SNOW_DENSITY}" \
    julia -O3 --check-bounds=no --math-mode=fast --threads "${THREAD_COUNT}" "${SCRIPT_PATH}" \
    --forcing-file="${FORCING_PATH}" \
    "${EXTRA_ARGS_ARR[@]}"

wait
