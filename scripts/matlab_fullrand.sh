#!/bin/bash
#SBATCH --job-name=fullrand_cpu
#SBATCH --output=fullrand.txt
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=20G
#SBATCH --partition=week
module purge
module load MATLAB/2023a
matlab -nodisplay -nosplash -nodesktop -r "run('/home/ys792/qaed/test_time_fullrand.m');exit;"
