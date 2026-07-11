#!/bin/bash
#SBATCH --job-name=qaed_paper
#SBATCH --output=paper_bench.txt
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

for type in full hess; do
    for n in 64 128 256 512 1024 2048 4096; do
        echo "===== type=$type n=$n aed ====="
        ./qaed_bench --alg aed --type "$type" --n "$n"
    done
    for n in 64 128 256 512 1024; do
        echo "===== type=$type n=$n iqr ====="
        ./qaed_bench --alg iqr --type "$type" --n "$n"
    done
done
