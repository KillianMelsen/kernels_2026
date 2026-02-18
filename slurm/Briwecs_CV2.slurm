#!/bin/bash

# -----------------------------Name of the job-------------------------
#SBATCH --job-name=kernels_2026_Briwecs_CV2

#-----------------------------Output files-----------------------------
#SBATCH --output=output/output_Briwecs_CV2_%A.%a.txt
#SBATCH --error=output/errors_Briwecs_CV2_%A.%a.txt

#-----------------------------Other information------------------------
#SBATCH --comment='kernels_2026 Briwecs CV2'

#-----------------------------Required resources-----------------------
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1-25

#-----------------------------Environment, Operations and Job steps----
# Modules
module reset
module load 2024
module load R/4.4.2
module load R-bundle-CRAN/2024

# Converting to unix and running R script
dos2unix ./Briwecs/scripts/Briwecs_CV2_${SLURM_ARRAY_TASK_ID}.R
chmod +x ./Briwecs/scripts/Briwecs_CV2_${SLURM_ARRAY_TASK_ID}.R
./Briwecs/scripts/Briwecs_CV2_${SLURM_ARRAY_TASK_ID}.R


