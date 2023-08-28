# Author:
Nhung Pham, 14-03-2023

# Data description

# Set up
## All steps 
Module 1, 2 and 4 were run on HPC with slurm in a conda environment.  A similar environment can be created from environment.yaml. Note: if there is packages conflict, have trim_galore installed before other packages. 

```
conda env create -f environment_export.yaml
```
Activate conda env by

```
conda activate cutnrun_env
```
Module 3,5,6 were run in a local machine in R.
## Data normalization step

To install s3norm. cd to the location where you want to save s3norm and then clone 

```
git clone https://github.com/guanjue/S3norm.git
```

# Other data
fasta genome and bowtie2index for hg38 are included in data folder

# How to run
The pipeline is not completely automatic. Adapting the steps to your experiment before running. The first 6 steps in module 1 are automatic. In module 4, a sample sheet need to be generated. Example sample sheet can be found in the data folder.

The main script is run_pipeline.sh. Read the script carefully and modify wherever needed as mentioned in the script. 

Comment out steps that are not desirable. 

Before run run_pipeline.sh, setup the Config.sh to suit your experiment.
# Runtime overview
Steps that are not mentioned here have runtime less than 30 minutes. 

The running time was recorded with --mem=90G, --cpu-per-task=8, and --time=96:0:0. Samples are run in parallel. These steps composed different modules. Check each module script for more detail. 

|Steps | Run time | Script name|
|------|-----------------|------------|
|1. Quality checking | ~ 30 mins (1 sample) | 1-qualityCheck.sh |
|2. Trimming|  ~40 mins (1 sample) |2-trimming.sh |
|3. Alignment|  ~ 30 mins - 2 hours (1 sample), 8 hours (18 samples) | 3-alignment.sh |
|4. Remove duplicate | ~ 20 mins (1 sample) | 4-filtering.sh |
|5. Convert bam file to bigwig| | 5-bam2bigwig.sh|
|6. Calculate FRiP| ~ 40 mins (24 samples) | 6-Calculate_FRiP.sh|
|7. Data normalization| ~ 4 hours (24 samples) | 7-run_s3norm.sh  |
|5. Call peak|  ~ 1 hour (18 samples) | 5-peakCalling.sh |
|6. Peak analysis| ~ 30 minutes | 6-DiffBind.R  |
|7. Peak annotation | |  |
|8. Transform bam file to bigwig | ~ 2-3 hours| |8-bam2bigwig.sh|
|9. Heatmap generation| ~30 minutes (4 samples) |9-heatmap.sh |
|10. Motif analysis| ~ 30 minutes - 1h ||
|11. Find motif location | ~ 30 minutes|  |



# Steps description

Report ideas were adapted from https://nf-co.re/cutandrun/dev/output#4--alignment-post-processing [https://github.com/nf-core/cutandrun]
## Module 1. Data processing: quality checking, trimming, alignment, replicates correlation 

<p>
<img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/M1.png" alt>
</p>
<p>
    <em>Module 1. Data processing and alignment<em>.
        </p>
        
### 1. Quality checking 
Reads quality for each sequence in each sample is checked with fastQC in order to identify poor quality sequencing sample(s). 

### 2. Trimming
Remove adapter and conduct fastQC after with trim_galore

<p>
<img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/QC_figure.png" width="400" height="300" alt>
</p>
<p>
    <em>Sequence quality before and after trimming. Before trimming all samples already have good quality sequences with most in q30 region (< 0.01% probablility for error). After trimming the quality is more even among samples<em>.
        </p>

### 3. Alignment: map reads to human genome g38
Reads were aligned to hg38 genome.
 <p>
<img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/alignment_report_w_new_samples.png" width="400" height="300" alt>
</p>
<p>
    <em>Alignment report<em>.
        </p>  
        
  <img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/alignment_length_report.png" width="400" height="300" alt>
</p>
<p>
    <em>Alignment length report<em>.
        </p>   
        
### 4. Remove duplicate: 
 Sequence duplicates and reads < 20bp were removed
<p>
        </p>
<img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/picard_deduplication.png" width="400" height="300">
<p>
    <em> Duplication report<em>
        </p>
        
<p>
<img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/duplication_rate_report.png" width="400" height = "300" alt>
    </p>
    <p> 
    <em>Duplication rate report</em>
</p>

### 5. Convert bam files to bigwig 
Aligned data after removing duplicates and short reads were converted to bigwig to be used in IGV for visualization. Replicates from the same condition are merged.

### 6. Calculate fraction of read in peak (FRiP) to check signal to noise ratio among samples
FRiP is calculated as the divison of reads in peaks to total reads in a sample.

<p>
<img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/FRIP.png" width="400" height = "300" alt>
    </p>
    <p> 
    <em>FRiP and total reads. The total reads range from 726,760 to 101,240,888 reads and FRiP score from 0.02 to 0.4. Data have a large variability in both sequencing depth and signal-to-noise ratio. These data need to be normalized before conducting any comparison or peak calling.</em>
</p>

### 7. Check replicates correlation
To check the correlation between replicates from the same group. A good experiment should have high correlation between replicates.
<p>
    <img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/heatmap_PearsonCorr_readCounts.png" width="400", height="400" alt>
</p>
<p>
    <em>Replicate correlation: overall most replicates from the same condition are highly correlated indicate that data from replicates are reliable.</em>
</p>

## Module 2. Data normalization and peak calling 
        
<p>
    <img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/M2.drawio.png" alt>
</p>
<p>
    <em> Module 2. Data normalization and peak calling. Steps, inputs and outputs.</em>
</p>

## Module 3. Peak analysis: differential analysis, occupancy analysis and annotation 

<p>
    <img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/M3.drawio.png" alt>
</p>
<p>
    <em> Module 3. Peak analysis: differential analysis, occupancy analysis and annotation. Steps, inputs and outputs.</em>
</p>


### Broad and narrow peak calling with and without control
<p>
    <img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/peaks_number_with_nocontrol_report.png" width="400", height="400" alt>
</p>
<p>
    <em>Number of peak per sample without control</em>
</p>

 ### Peak reproducibility
 Compare peaks from replicates in the same group. In a good experiment replicates from the same group should generate the same peaks. 
 <p>
    <img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/H3k4me3_samples_peak_counts.png" width="500", height="500" alt>
</p>
<p>
    <em>Peaks reproducibility in H3k4me3 samples. These samples have the same antibody (anti-H3K4me3) however the seq depth are varied among them hence some samples generate more peaks than others. Most peaks from less depth samples are present in the depth samples, therefore the difference in peaks are not a problem</em>
</p>
        <p>
    <img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/tfe3_fusion_luc_samples_peak_counts.png" width="400", height="400" alt>
</p>
<p>
    <em>Peaks reproducibility in Tfe3, luc and fusion samples. Samples in the same group have the same antibody but different seq depth. However, only less than 40% of peak are overlapped among replicates. </em>
</p>
       <p>
    <img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/peak_counts_after_s3norm_normalize.png" width="400", height="400" alt>
</p>
<p>
    <em>Peaks reproducibility in Tfe3, luc and fusion samples. Peak counts among replicates after normalizing data with s3norm are better than that from raw data.  </em>
</p>


## Module 4. Motif analysis  
<p>
    <img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/M4.drawio.png" alt>
</p>
<p>
    <em> Module 4. Motif analysis. Steps, inputs and outputs.</em>
</p>

## Module 5. Enrichment analysis: GO bioprocess, GO molecular function, GSEA
<p>
    <img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/M5.drawio.png" alt>
</p>
<p>
    <em> Module 5. Enrichment analysis. Steps, inputs and outputs.</em>
</p>

## Module 6. Data integration
Expression data such as bulkRNA or single cell RNA are integrated to study the effect of protein bindings on genes in the identified peaks.  

<p>
    <img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/M6.drawio.png" alt>
</p>
<p>
    <em> Module 6. Data integration. Steps, inputs and outputs.</em>
</p>

## Module 7. Result visualization
<p>
    <img src="https://github.com/nhungpham1707/CUTnRUN/blob/main/Figures/M7.png" alt>
</p>
<p>
    <em> Module 7. Result visualization. Steps, inputs and outputs.</em>
</p>
