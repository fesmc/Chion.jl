#!/usr/bin/env bash
# All model runs behind the paper figures, with the configurations in examples/configs and the
# scripts in examples/runs and examples/research. One job runs one group; the runs of a group share
# the node and its GPU (NPAR at a time). Each run logs to logs/paper/<name>.log; finished runs are skipped.
#
#   sbatch run_paper_experiments.sh greenland    # transients and firn: E1-E3, D1-D3 (5 runs)
#   sbatch run_paper_experiments.sh native       # transient with native 3-hourly MAR: E1-E3 (1 run)
#   sbatch run_paper_experiments.sh antarctica   # transients and firn: A1, A3, D1-D3 (5 runs)
#   sbatch --partition=standard --qos=short --gres=gpu:0 run_paper_experiments.sh snow   # SnowMIP sites, NH ERA5 snow
#
# Not in this script: the performance benchmarks (examples/research/performance_benchmark.jl and
# performance_greenland.jl, timed on their own nodes), the longwave fit (fit_cloud_proxy_longwave.jl)
# and the figures without model runs (verification_suite.jl, plot_diurnal_shortwave_cycle.jl,
# model_schematic.py). The NH run needs the daily ERA5 forcing (CHION_ERA5_NH_DIR).
#SBATCH --qos=gpushort
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=chion_paper
#SBATCH --output=/p/projects/ou/labs/ai/Nils/Chion.jl/logs/paper-%x-%j.log
set -uo pipefail
module load julia/1.12.2 hdf5 netcdf-c
PROJECT_DIR=/p/projects/ou/labs/ai/Nils/Chion.jl
export JULIA_DEPOT_PATH="/home/bochow/.julia:${PROJECT_DIR}/.julia-depot:" JULIA_PKG_PRECOMPILE_AUTO=0 GKSwstype=100
cd "${PROJECT_DIR}"
mkdir -p logs/paper
GROUP=$1
C=examples/configs
CONFIG="-m Chion"
PRESCRIBED="--set forcing.surface_fluxes=prescribed --set model.albedo=prescribed"

if [[ $GROUP == native ]]; then
    # Spin-up: 5 cycles of the daily MAR years 1940-1969, then 1940-1979 on daily MAR, then the
    # native 3-hourly MAR 1980-2017 (all with prescribed surface fluxes and albedo).
    export CHION_BACKEND=gpu CHION_SPINUP_MODE=cycled CHION_SPINUP_FIRST=1940 CHION_SPINUP_LAST=1969 CHION_SPINUP_CYCLES=5
    export CHION_DIURNAL_TEMPERATURE_CONSTANT_AMPLITUDE_C=1.0
    export CHION_OUTPUT_DIR=output/gris_substrate5_suite/transient_native_3hourly
    export CHION_NATIVE_MONTHLY_DIR=${CHION_OUTPUT_DIR}/monthly
    mkdir -p "${CHION_NATIVE_MONTHLY_DIR}"
fi

# name | marker file of a finished run | arguments after `julia --project=. ...`
runs() {
    case $GROUP in
    greenland)
        echo "greenland_prescribed_steps_8|output/greenland_prescribed_steps_8/budget_yearly.csv|$CONFIG $C/greenland_prescribed.toml --backend=gpu --output-dir=output/greenland_prescribed_steps_8 --set model.diurnal_shortwave_max_substeps=8"
        echo "greenland_parameterized|output/greenland_parameterized/budget_yearly.csv|$CONFIG $C/greenland_parameterized.toml --backend=gpu"
        echo "greenland_firn_parameterized|output/greenland_firn_parameterized/firn_profiles.nc|$CONFIG $C/greenland_firn.toml --backend=gpu --output-dir=output/greenland_firn_parameterized"
        echo "greenland_firn_parameterized_ntot30|output/greenland_firn_parameterized_ntot30/firn_profiles.nc|$CONFIG $C/greenland_firn.toml --backend=gpu --output-dir=output/greenland_firn_parameterized_ntot30 --set model.Ntot=30"
        echo "greenland_firn_prescribed|output/greenland_firn_prescribed/firn_profiles.nc|$CONFIG $C/greenland_firn.toml --backend=gpu --output-dir=output/greenland_firn_prescribed $PRESCRIBED"
        ;;
    native)
        echo "greenland_native_3hourly|output/gris_substrate5_suite/transient_native_3hourly/native_3hourly_transient_1980_2017_budgets.csv|examples/research/run_gris_native_3hourly_transient.jl"
        ;;
    antarctica)
        for m in parameterized prescribed; do
            set_m=""; [[ $m == prescribed ]] && set_m=$PRESCRIBED
            echo "antarctica_$m|output/antarctica_$m/budget_yearly.csv|$CONFIG $C/antarctica_$m.toml --backend=gpu"
            echo "antarctica_firn_$m|output/antarctica_firn_$m/firn_profiles.nc|$CONFIG $C/antarctica_firn.toml --backend=gpu --output-dir=output/antarctica_firn_$m $set_m"
        done
        echo "antarctica_firn_parameterized_ntot30|output/antarctica_firn_parameterized_ntot30/firn_profiles.nc|$CONFIG $C/antarctica_firn.toml --backend=gpu --output-dir=output/antarctica_firn_parameterized_ntot30 --set model.Ntot=30"
        ;;
    snow)
        echo "snowmip_sites|output/snowmip/runs/cdp.csv|examples/runs/snowmip_sites.jl"
        echo "nh_era5_snow|output/nh_era5_snow/monthly_snow.nc|examples/runs/nh_era5_snow.jl"
        ;;
    *) echo "unknown group $GROUP"; exit 1 ;;
    esac
}

case $GROUP in native) NPAR=1 ;; snow) NPAR=2 ;; *) NPAR=5 ;; esac
THREADS=$(( SLURM_CPUS_PER_TASK / NPAR > 1 ? SLURM_CPUS_PER_TASK / NPAR : 2 ))
run_one() {
    IFS='|' read -r name marker args <<< "$1"
    [[ -f $marker ]] && { echo "skip $name"; return 0; }
    julia --project=. --startup-file=no -O3 --check-bounds=no --threads "$THREADS" $args > "logs/paper/$name.log" 2>&1
    echo "$name exit $?"
}
export -f run_one
export THREADS
runs | xargs -d '\n' -P "$NPAR" -I{} bash -c 'run_one "$@"' _ {}
