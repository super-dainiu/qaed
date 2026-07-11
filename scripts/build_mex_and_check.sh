#!/bin/bash
#SBATCH --job-name=qaed_mex
#SBATCH --output=mex_check.txt
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=day
module purge
module load MATLAB/2023b

# sbatch runs a spooled copy of this script, so $0 is useless for locating the
# repo; fall back to the hardcoded path when not submitted from the repo root.
ROOT="${SLURM_SUBMIT_DIR:-/nfs/roberts/project/pi_mg269/ys792/qaed}"
[ -d "$ROOT/cpp" ] || ROOT=/nfs/roberts/project/pi_mg269/ys792/qaed
cd "$ROOT/cpp"

# MATLAB bundles an old libstdc++ that lacks the GLIBCXX version the mex file
# (built with the system g++) needs; preload the system one instead.
export LD_PRELOAD=/usr/lib64/libstdc++.so.6

matlab -nodisplay -nosplash -nodesktop -batch "\
mex('-R2018a','-DQAED_ILP64','CXXFLAGS=\$CXXFLAGS -O3 -std=c++17','qaed_mex.cpp','-lmwlapack','-lmwblas'); \
disp('mex build ok'); \
run('$ROOT/test_cpp_vs_matlab.m');"
