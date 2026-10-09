#!/bin/bash
#SBATCH --qos=gpushort
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0-12:00:00
#SBATCH --job-name=chion_column_bench
#SBATCH --gres=gpu:1
#SBATCH --output=/p/projects/ou/labs/ai/Nils/Chion.jl/logs/snowpack-column-benchmark-%j.log
#SBATCH --error=/p/projects/ou/labs/ai/Nils/Chion.jl/logs/snowpack-column-benchmark-%j.err

set -euo pipefail

module load julia

PROJECT_DIR=/p/projects/ou/labs/ai/Nils/Chion.jl
SCRIPT_PATH="${PROJECT_DIR}/examples/scripts/benchmark_snowpack_columns.jl"
COLUMNS="${COLUMNS:-128,1024,16384,65536,524288,1048576}"
NTOTS="${NTOTS:-4,8,12,20}"
BACKENDS="${BACKENDS:-gpu}"
THREADS="${THREADS:-8}"
YEARS="${YEARS:-100}"
STEPS="${STEPS:-365}"
REPETITIONS="${REPETITIONS:-3}"
OUTPUT="${OUTPUT:-${PROJECT_DIR}/logs/snowpack-columns-cpu8.csv}"

export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-${PROJECT_DIR}/.julia-depot:${HOME}/.julia:}"
export JULIA_PKG_PRECOMPILE_AUTO=0
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
mkdir -p "${PROJECT_DIR}/logs"
cd "${PROJECT_DIR}"

srun julia -O3 --project="${PROJECT_DIR}" "${SCRIPT_PATH}" \
    --columns="${COLUMNS}" \
    --ntots="${NTOTS}" \
    --backends="${BACKENDS}" \
    --threads="${THREADS}" \
    --years="${YEARS}" \
    --steps="${STEPS}" \
    --repetitions="${REPETITIONS}" \
    --output="${OUTPUT}"
