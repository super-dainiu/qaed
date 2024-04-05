#!/bin/bash
#SBATCH --job-name=hessrand_cpu
#SBATCH --output=hessrand_err.txt
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=20G
#SBATCH --partition=week
module purge
module load MATLAB/2023a
matlab -nodisplay -nosplash -nodesktop -r "run('/home/ys792/qaed/test_relerr_hessrand.m');exit;"
