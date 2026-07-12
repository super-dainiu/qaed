#!/bin/bash
#SBATCH --job-name=qtime3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=23:00:00
#SBATCH --constraint=cpumodel:8562Y+
# Non-exclusive best-of-3 timing: schedules immediately; min-of-3 rejects
# transient node contention. Large n (>=4096) run once (compute dominates).
module purge
module load GCCcore/13.3.0
ROOT="${SLURM_SUBMIT_DIR:-/nfs/roberts/project/pi_mg269/ys792/qaed}"
[ -d "$ROOT/cpp" ] || ROOT=/nfs/roberts/project/pi_mg269/ys792/qaed
cd "$ROOT/cpp"
export OMP_NUM_THREADS=4
lscpu | grep 'Model name'
TYPE=$1
for alg in aed iqr; do
  for n in 64 128 256 512 1024 2048 4096 8192; do
    reps=3; [ "$n" -ge 4096 ] && reps=1
    for r in $(seq 1 $reps); do
      echo "===== type=$TYPE n=$n $alg rep=$r ====="
      ./qaed_bench --alg $alg --type "$TYPE" --n "$n"
    done
  done
done
