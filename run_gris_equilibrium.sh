#!/bin/bash

#SBATCH --qos=gpushort
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
##SBATCH --time=2-1:00:00
#SBATCH --time=0-0:30:00
#SBATCH --job-name=gris_equilibrium
#SBATCH --gres=gpu:1

# Don't change anything below this line:
#SBATCH --output=/p/projects/ou/labs/ai/Nils/Chion.jl/logs/gris-%j.log
#SBATCH --error=/p/projects/ou/labs/ai/Nils/Chion.jl/logs/gris-%j.err

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
SCRIPT_PATH="${PROJECT_DIR}/examples/scripts/run_gris_equilibrium.jl"
LOG_DIR="${PROJECT_DIR}/logs"

# Optional: set this to your MAR forcing file before submission, or pass it
# inline with:
# sbatch --export=ALL,NC_PATH=/path/to/file.nc run_gris_equilibrium.sh
NC_PATH="${NC_PATH:-/p/projects/ou/labs/ai/Nils/MARv3.14.3-10km-daily-ERA5-2025.nc}"
EXTRA_ARGS="${EXTRA_ARGS:-}"

mkdir -p "${LOG_DIR}"
cd "${PROJECT_DIR}"

EXTRA_ARGS_ARR=()
HAS_BACKEND_ARG=0
if [[ -n "${EXTRA_ARGS}" ]]; then
    read -r -a EXTRA_ARGS_ARR <<< "${EXTRA_ARGS}"
    for arg in "${EXTRA_ARGS_ARR[@]}"; do
        if [[ "${arg}" == --backend=* ]]; then
            HAS_BACKEND_ARG=1
            break
        fi
    done
fi

if [[ "${HAS_BACKEND_ARG}" -eq 0 ]]; then
    EXTRA_ARGS_ARR+=(--backend=gpu)
fi

if [[ -n "${NC_PATH}" ]]; then
    srun julia -O3 --check-bounds=no --math-mode=fast --threads 8 "${SCRIPT_PATH}" --nc="${NC_PATH}" "${EXTRA_ARGS_ARR[@]}"
else
    srun julia -O3 --check-bounds=no --math-mode=fast --threads 8 "${SCRIPT_PATH}" "${EXTRA_ARGS_ARR[@]}"
fi

wait
