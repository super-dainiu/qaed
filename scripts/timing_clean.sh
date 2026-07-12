#!/bin/bash
#SBATCH --job-name=qtime
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --time=23:00:00
#SBATCH --constraint=cpumodel:8562Y+
module purge
module load GCCcore/13.3.0
ROOT="${SLURM_SUBMIT_DIR:-/nfs/roberts/project/pi_mg269/ys792/qaed}"
[ -d "$ROOT/cpp" ] || ROOT=/nfs/roberts/project/pi_mg269/ys792/qaed
cd "$ROOT/cpp"
export OMP_NUM_THREADS=4
lscpu | grep 'Model name'
TYPE=$1
for n in 64 128 256 512 1024 2048 4096 8192; do
  echo "===== type=$TYPE n=$n aed ====="
  ./qaed_bench --alg aed --type "$TYPE" --n "$n"
done
for n in 64 128 256 512 1024 2048 4096 8192; do
  echo "===== type=$TYPE n=$n iqr ====="
  ./qaed_bench --alg iqr --type "$TYPE" --n "$n"
done
