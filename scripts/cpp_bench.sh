#!/bin/bash
#SBATCH --job-name=qaed_cpp
#SBATCH --output=cpp_bench.txt
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=day
module purge
module load GCCcore/13.3.0

# sbatch runs a spooled copy of this script, so $0 is useless for locating the
# repo; fall back to the hardcoded path when not submitted from the repo root.
ROOT="${SLURM_SUBMIT_DIR:-/nfs/roberts/project/pi_mg269/ys792/qaed}"
[ -d "$ROOT/cpp" ] || ROOT=/nfs/roberts/project/pi_mg269/ys792/qaed
cd "$ROOT/cpp"
make clean && make

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-10}
./qaed_bench --selftest || exit 1

for n in 64 128 256 512 1024 2048; do
    echo "===== n = $n (aed) ====="
    ./qaed_bench --alg aed --n "$n"
done
for n in 64 128 256 512 1024; do
    echo "===== n = $n (iqr) ====="
    ./qaed_bench --alg iqr --n "$n"
done
