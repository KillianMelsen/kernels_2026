This repository contains all data and code required to reproduce the results from the REML kernels paper.

# Instructions for reproducing the results
## BRIWECS
The following steps can be followed to reproduce all results for the BRIWECS dataset:
* Run Briwecs/raw_data/create_K.R to compute the kinship matrix.
* Run Briwecs/scripts/BRIWECS_create_data.R to do the single-trial analyses and compute the kernels.
* Run all Briwecs/scripts/Briwecs_CV2_{1-25}.R and Briwecs/scripts/Briwecs_CV2_MB_{1-25}.R scripts to fit the different models to the different training/test sets for the different levels of sparsity.
* Run Briwecs/scripts/Briwecs_CV2_MB_merge_bandwidths.R, Briwecs/scripts/Briwecs_CV2_MB_merge.R, and Briwecs/scripts/Briwecs_CV2_merge.R to merge all results into single files.
* Run Briwecs/scripts/Briwecs_CV2_Plotting.R to create the figures from the paper.
* Run BRIWECS_LOF/BRIWECS_LOF.R to reproduce all results and figures related to the variance partitioning for the BRIWECS dataset.

## DROPS
The following steps can be followed to reproduce all results for the DROPS dataset:
* Run DROPS/scripts/DROPS_create_data.R to compute the kinship matrix, do the single-trial analyses, and compute the kernels.
* Run all DROPS/scripts/DROPS_CV2_{1-15}.R and DROPS/scripts/DROPS_CV2_MB_{1-15}.R scripts to fit the different models to the different training/test sets for the different levels of sparsity.
* Run DROPS/scripts/DROPS_CV2_MB_merge_bandwidths.R, DROPS/scripts/DROPS_CV2_MB_merge.R, and DROPS/scripts/DROPS_CV2_merge.R to merge all results into single files.
* Run DROPS/scripts/DROPS_CV2_Plotting.R to create the figures from the paper.
* Run DROPS_LOF/DROPS_LOF.R to reproduce all results and figures related to the variance partitioning for the DROPS dataset.

## Other
* misc_plotting.R can be used to generate the bandwidth visualization figure as well as the maps from the paper.
* The slurm folder contains slurm scripts for running the genomic prediction scripts in parallel on an HPC cluster.
* The output folder will contain output from these slurm jobs.
* The plots folder contains all generated plots and figures.
