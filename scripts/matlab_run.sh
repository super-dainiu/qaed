#!/bin/bash
#SBATCH --job-name=qaed_matlab
#SBATCH --output=matlab_run.txt
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=week
#
# Usage: sbatch scripts/matlab_run.sh <test script in tests/, without .m>
# e.g.   sbatch scripts/matlab_run.sh test_time_fullrand
module purge
module load MATLAB/2023b

ROOT="${SLURM_SUBMIT_DIR:-/nfs/roberts/project/pi_mg269/ys792/qaed}"
[ -d "$ROOT/cpp" ] || ROOT=/nfs/roberts/project/pi_mg269/ys792/qaed
cd "$ROOT"

# MATLAB bundles an old libstdc++; the mex file needs the system one.
export LD_PRELOAD=/usr/lib64/libstdc++.so.6

matlab -nodisplay -nosplash -nodesktop -batch "run('$ROOT/tests/${1}.m');"
