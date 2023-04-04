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

|Steps | Description | Run time (with current set up, 16 threads, parallel)| Script name|
|------|-------------|-----------------|------------|
|1. Quality checking  | check reads quality with fastQC | ~ 30 mins (1 sample) | 1-qualityCheck.sh [before_trimming_fastqc_all_samples.sh]|
|2. Trimming| remove adapter and conduct fastQC after with trim_galore | ~40 mins (1 sample) |2-trimming.sh |
|3. Alignment| map reads to human genome g38| ~ 30 mins - 2 hours (1 sample), 8 hours (18 samples) | 3-alignment.sh |
|4. Remove duplicate | remove duplicate, reads < 20 |~ 20 mins (1 sample) | 4-filtering.sh |
|5. Call peak| broad and narrow peak calling with and without control | ~ 1 hour (18 samples) | 5-peakCalling.sh |
|6. Peak analysis| identify differential binding peaks between groups| | 6-DiffBind.R  |
|7. Peak annotation | | | |
|8. Transform bam file to bigwig | merge bam files of the same condition and convert to bigwig to use for heatmap generation | ~ 2-3 hours |8-bam2bigwig.sh|
|9. Prepare for motif finding| merge narrow peak files from the same condition and extract the fasta sequence for STREME meme suit | a few seconds |8-prepareMotifAnalysis.sh |
|10. Find motif | | | |
|11. Analyze motif | | | | 




