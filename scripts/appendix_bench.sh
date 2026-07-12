#!/bin/bash
#SBATCH --job-name=qaed_appx
#SBATCH --output=appendix_bench.txt
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=week
#SBATCH --constraint=cpumodel:8562Y+
module purge
module load MATLAB/2023b

ROOT="${SLURM_SUBMIT_DIR:-/nfs/roberts/project/pi_mg269/ys792/qaed}"
[ -d "$ROOT/cpp" ] || ROOT=/nfs/roberts/project/pi_mg269/ys792/qaed
cd "$ROOT/cpp"

# MATLAB bundles an old libstdc++; the mex file needs the system one.
export LD_PRELOAD=/usr/lib64/libstdc++.so.6

matlab -nodisplay -nosplash -nodesktop -batch "\
mex('-R2018a','-DQAED_ILP64','CXXFLAGS=\$CXXFLAGS -O3 -std=c++17 -fopenmp','LDFLAGS=\$LDFLAGS -fopenmp','qaed_mex.cpp','-lmwlapack','-lmwblas'); \
disp('mex build ok'); \
maxNumCompThreads(4); \
run('$ROOT/tests/test_appendix.m');"
