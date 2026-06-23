#!/bin/bash

# -----------------------------Name of the job-------------------------
#SBATCH --job-name=MB_kernels_2026_DROPS_CV2_withweights

#-----------------------------Output files-----------------------------
#SBATCH --output=output/output_DROPS_MB_CV2_%A.%a.txt
#SBATCH --error=output/errors_DROPS_MB_CV2_%A.%a.txt

#-----------------------------Other information------------------------
#SBATCH --comment='MB kernels_2026 DROPS CV2 with weights'

#-----------------------------Required resources-----------------------
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1

#-----------------------------Environment, Operations and Job steps----
# Modules
module reset
module load 2024
module load R/4.4.2
module load R-bundle-CRAN/2024

# Converting to unix and running R script
dos2unix ./DROPS/scripts/DROPS_CV2_MB_${SLURM_ARRAY_TASK_ID}.R
chmod +x ./DROPS/scripts/DROPS_CV2_MB_${SLURM_ARRAY_TASK_ID}.R
./DROPS/scripts/DROPS_CV2_MB_${SLURM_ARRAY_TASK_ID}.R


