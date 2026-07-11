#!/bin/bash
#SBATCH --job-name=qaed_skew
#SBATCH --output=skew_bench.txt
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=devel
module purge
module load GCCcore/13.3.0

ROOT="${SLURM_SUBMIT_DIR:-/nfs/roberts/project/pi_mg269/ys792/qaed}"
[ -d "$ROOT/cpp" ] || ROOT=/nfs/roberts/project/pi_mg269/ys792/qaed
cd "$ROOT/cpp"
make clean && make

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
lscpu | grep 'Model name'
./qaed_bench --selftest || exit 1

for n in 64 128 256 512 1024 2048 4096 8192; do
    for alg in skew_aed skew_iqr; do
        echo "===== type=skew n=$n $alg ====="
        ./qaed_bench --alg "$alg" --n "$n"
    done
done
