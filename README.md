This repository contains all data and code required to reproduce the results from the REML kernels paper.
# Folder/file structure
* BRIWECS_LOF contains a script and results for the variance partitioning for the Briwecs dataset.
* Briwecs/raw_data contains all raw Briwecs data as well as a script for creating the kinship matrix.
* Briwecs/scripts contains all scripts for pre-processing the dataset, fitting the models in parallel using SLURM, merging results, and plotting results.
* Briwecs/data contains the pre-processed Briwecs data.
* Briwecs/results/CV2 contains all result files for the Briwecs data.
* DROPS/raw_data contains all raw DROPS data.
* DROPS/scripts contains all scripts for pre-processing the dataset, fitting the models in parallel using SLURM, merging results, and plotting results.
* DROPS/data contains the pre-processed DROPS data.
* DROPS/results/CV2 contains all result files for the DROPS data.
* DROPS_LOF contains a script and results for the variance partitioning for the DROPS dataset.
* output is a folder with SLURM output.
* slurm contains SLURM scripts for submitting the jobs.
* misc_plotting.R produces some additional plots.
