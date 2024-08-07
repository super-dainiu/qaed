#!/bin/bash
#SBATCH --job-name=skewrand_cpu
#SBATCH --output=skewrand.txt
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=20G
#SBATCH --partition=week
module purge
module load MATLAB/2023a
matlab -nodisplay -nosplash -nodesktop -r "run('/home/ys792/project/qaed/test_time_skewrand.m');exit;"
