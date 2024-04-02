#!/bin/bash
#SBATCH --job-name=fullrand_cpu
#SBATCH --output=fullrand.txt
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=32
module purge
module load MATLAB/R2023b
matlab -nodisplay -nosplash -nodesktop -r "run('/home/yjshao/MATLAB/test_time_fullrand.m');exit;"
