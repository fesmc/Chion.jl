#!/bin/bash
# Run an experiment configuration on one GPU:
#   sbatch run_config_gpu.sh examples/configs/greenland_parameterized.toml
#SBATCH --qos=gpushort
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=0-08:00:00
#SBATCH --job-name=chion
#SBATCH --gres=gpu:1
#SBATCH --output=/p/projects/ou/labs/ai/Nils/Chion.jl/logs/chion-%j.log
set -euo pipefail
module load julia/1.12.2 hdf5 netcdf-c
PROJECT_DIR=/p/projects/ou/labs/ai/Nils/Chion.jl
export JULIA_DEPOT_PATH="/home/bochow/.julia:${PROJECT_DIR}/.julia-depot:" JULIA_PKG_PRECOMPILE_AUTO=0 GKSwstype=100
cd "${PROJECT_DIR}"
julia --project=. --startup-file=no -O3 --check-bounds=no --threads "${SLURM_CPUS_PER_TASK}" \
    -m Chion "$@" --backend=gpu
