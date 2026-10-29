#!/bin/bash
#SBATCH --job-name=gusto_held_suarez
#SBATCH --output=logs/gusto_held_suarez_%j.out
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=120
#SBATCH --mem=0
#SBATCH --distribution=block:block
#SBATCH --exclusive

set -euo pipefail

example='test_scripts/held_suarez.py'
rundir='/data/users/thomas.bendall/gusto_gmd_paper_results'
procs=120
n_cycles=12

###############################################################################
# Queue the remaining cycles. Only the first job (PICKUP unset) does this,
# chaining each cycle to the previous with afterok, so one sbatch is enough.
###############################################################################
if [[ -z "${PICKUP:-}" ]]; then
    script=$(scontrol show job "${SLURM_JOB_ID}" | grep -oP 'Command=\K\S+')
    jid=${SLURM_JOB_ID}
    for i in $(seq 2 "${n_cycles}"); do
        jid=$(sbatch --parsable --export=ALL,PICKUP=1 \
              --dependency=afterok:${jid} --kill-on-invalid-dep=yes "${script}")
        echo "Queued cycle ${i}: job ${jid}"
    done
fi

###############################################################################
# Modules and activate environment
###############################################################################

module purge
ml use ~lfricadmin/lmod
ml lfric
export LD_LIBRARY_PATH=/data/apps/spack/releases/2026.7.1/install/linux-zen2/mpich-4.3.2-milbtch3bj6p7e363urmulyl6aslipy3/lib
export OMP_NUM_THREADS=1

###############################################################################
# Firedrake / PyOP2 caches (portable, AZ‑Spice safe)
###############################################################################
if [[ -n "${SLURM_TMPDIR:-}" ]]; then
    CACHE_BASE="${SLURM_TMPDIR}/firedrake_cache"
elif [[ -n "${TMPDIR:-}" ]]; then
    CACHE_BASE="${TMPDIR}/firedrake_cache"
elif [[ -d "/scratch/${USER}" ]]; then
    CACHE_BASE="/scratch/${USER}/firedrake_cache"
else
    CACHE_BASE="${HOME}/.firedrake_cache"
fi

export PYOP2_CACHE_DIR="${CACHE_BASE}/pyop2"
export FIREDRAKE_TSFC_KERNEL_CACHE_DIR="${CACHE_BASE}/tsfc"
export XDG_CACHE_HOME="${CACHE_BASE}/xdg"

mkdir -p \
  "${PYOP2_CACHE_DIR}" \
  "${FIREDRAKE_TSFC_KERNEL_CACHE_DIR}" \
  "${XDG_CACHE_HOME}"

echo "Using Firedrake cache base: ${CACHE_BASE}"

###############################################################################
# Run model
###############################################################################
cd ${rundir}
source /data/users/thomas.bendall/firedrakevenv/bin/activate

# PICKUP=1 is set by cycle 1 for every cycle after the first
pickup_flag=""
if [[ "${PICKUP:-0}" == "1" ]]; then
    pickup_flag="--pickup"
fi
mpiexec -n ${procs} python3 ${example} ${pickup_flag}
