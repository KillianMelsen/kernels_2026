#!/bin/bash

# -----------------------------Name of the job-------------------------
#SBATCH --job-name=MB_kernels_2026_Briwecs_CV2_withweights

#-----------------------------Output files-----------------------------
#SBATCH --output=output/output_Briwecs_MB_CV2_%A.%a.txt
#SBATCH --error=output/errors_Briwecs_MB_CV2_%A.%a.txt

#-----------------------------Other information------------------------
#SBATCH --comment='MB kernels_2026 Briwecs CV2 with weights'

#-----------------------------Required resources-----------------------
#SBATCH --time=4-00:00:00
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
dos2unix ./Briwecs/scripts/Briwecs_CV2_MB_${SLURM_ARRAY_TASK_ID}.R
chmod +x ./Briwecs/scripts/Briwecs_CV2_MB_${SLURM_ARRAY_TASK_ID}.R
./Briwecs/scripts/Briwecs_CV2_MB_${SLURM_ARRAY_TASK_ID}.R


