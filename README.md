## Author:
Nhung Pham, 14-03-2023

## Data description

## Set up
All scripts were run on HPC with a customized conda environment. Similar environment can be created from environment.yaml

## Scripts description

The cut and run data were analyzed based on published protocol with modification

## Steps:

1. FastQC: check sequences quality. Run time with the current setup in the script: 3-4 hours for 18 samples. Script: before_trimming_fastqc_all_samples.sh 
2. Trim: remove adapter with trim_galore  
