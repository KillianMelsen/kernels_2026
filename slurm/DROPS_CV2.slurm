#!/bin/bash

# -----------------------------Name of the job-------------------------
#SBATCH --job-name=kernels_2026_DROPS_CV2

#-----------------------------Output files-----------------------------
#SBATCH --output=output/output_DROPS_CV2_%A.%a.txt
#SBATCH --error=output/errors_DROPS_CV2_%A.%a.txt

#-----------------------------Other information------------------------
#SBATCH --comment='kernels_2026 DROPS CV2'

#-----------------------------Required resources-----------------------
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1-15

#-----------------------------Environment, Operations and Job steps----
# Modules
module reset
module load 2024
module load R/4.4.2
module load R-bundle-CRAN/2024

# Converting to unix and running R script
dos2unix ./DROPS/scripts/DROPS_CV2_${SLURM_ARRAY_TASK_ID}.R
chmod +x ./DROPS/scripts/DROPS_CV2_${SLURM_ARRAY_TASK_ID}.R
./DROPS/scripts/DROPS_CV2_${SLURM_ARRAY_TASK_ID}.R


