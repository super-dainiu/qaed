#!/bin/bash
#SBATCH --job-name=skewrand_cpu
#SBATCH --output=skewrand_err.txt
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=20G
#SBATCH --partition=day
module purge
module load MATLAB/2023a
matlab -nodisplay -nosplash -nodesktop -r "run('/home/ys792/project/qaed/test_relerr_skewrand.m');exit;"
