## Author:
Nhung Pham, 14-03-2023

## Data description

## Set up
All scripts were run on HPC with slurm in a conda environment. A similar environment can be created from environment.yaml. The file is arranged to have trim_galore installed before other packages to prevent version conflict. 

```
conda env create -f environment.yaml
```
## Steps description

The cut and run data were analyzed based on published protocol with modification [add ref]. The running time was calculated with the current setup in each script. 

|Steps | Description | Run time (hours)| Script name|
|------|-------------|-----------------|------------|
|1. Quality checking  | check reads quality with fastQC | 3-4| before_trimming_fastqc_all_samples.sh|
|2. Trimming| remove adapter and conduct fastQC after with trim_galore | | |
|3. Alignment| | | |
|4. Remove duplicate | | | |
|5. Call peak| | | |
|6. Peak annotation | | | |
|7. Find motif | | | |
|8. Transform bam file to bigwig | | | |
|9. Analyze motif | | | | 




