#!/bin/bash
#SBATCH --job-name=gusto_straka_0
#SBATCH --output=logs/gusto_straka_0_%j.out
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=0
#SBATCH --distribution=block:block
#SBATCH --exclusive

set -euo pipefail

rundir='/data/users/thomas.bendall/gusto_gmd_paper_results'
make_command=straka_degree_0

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
make ${make_command}
