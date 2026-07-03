#!/bin/bash

#SBATCH --qos=priority
##SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0-5:30:00
#SBATCH --job-name=gris_mar_case
##SBATCH --gres=gpu:1
#SBATCH --output=/p/projects/ou/labs/ai/Nils/Chion.jl/logs/gris-mar-%j.log
#SBATCH --error=/p/projects/ou/labs/ai/Nils/Chion.jl/logs/gris-mar-%j.err

set -euo pipefail

module load julia
module load hdf5
module load netcdf-c

PROJECT_DIR=/p/projects/ou/labs/ai/Nils/Chion.jl
SCRIPT_PATH="${PROJECT_DIR}/examples/scripts/run_gris_forcing_file_case.jl"
THREAD_COUNT="${SLURM_CPUS_PER_TASK}"

export JULIA_NUM_THREADS="${THREAD_COUNT}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export H5DUMP_BIN="$(which h5dump)"
export H5LS_BIN="$(which h5ls)"
export NETCDF_LIB="${NETCDF_LIB:-libnetcdf.so}"

mkdir -p "${PROJECT_DIR}/logs"
cd "${PROJECT_DIR}"

srun julia -O3 --check-bounds=no --math-mode=fast \
    --threads "${THREAD_COUNT}" \
    "${SCRIPT_PATH}"
