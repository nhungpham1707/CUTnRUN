#!/bin/bash
#SBATCH --job-name=align
#SBATCH --output=trimming_and_alignment.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=90G
#SBATCH --cpus-per-task=8
#SBATCH --gres=tmpspace:30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl

#==================================================================
#   This is a master script to analyze cut and run data. 
#   The script can also be adapted for ATAC, ChiPSeq or 
#   data that required similar steps. 
#
#   More infomation can be found at 
#   https://github.com/nhungpham1707/CUTnRUN
#   Author: Nhung, 01 05 2023
#==================================================================

#==================================================================
#                   TOOLS & FLAGS EXPLAINATION
#
#   Quality check:
#   - fastqc v0.12.1
#
#   Trimming: 
#   - trim_galore ver 0.6.10

#   Alignment: 
#	- bowtie ver 2.5.1, input flags:
#			-end-to-end: for trimmed sequence
#       (used local if seq is not trimmed before this step, sample_ID path will need to be modified to run if no trim was done)
#			-p 16: 16 threads to speed up the process
#			-samtools -bS: to save input sam file as bam file

#==================================================================

#==================================================================
#                       INPUT & OUTPUT DATA
#
#   Quality check 
#       Input: R1 and R2 fastq sequences for each sample
#       Output:
#       - output from each sample will be saved in a subfolder with the respective sample ID
#       - For each sample a zip file and a html file
#
#   Trimming 
#       Input: R1 and R2 fastq sequences for each sample
#       Output:
#           - 2 trimming reports for R1 and R2.
#           - 2 val_ files for R1 and R2 (inputs for the next step)
#           - 2 fastqc html for R1 and R2
#           - 2 fastqc.zip files for R1 and R2.
#
#   Alignment
#       Input: R1 and R2 outputs from trimming for each sample
#       Output:
#           - 1 bam file per sample stored in its respective folder 

#==================================================================

#==================================================================
#                   DEFINE GLOBAL VARIABLES
# Directory where the sample fastq are stored

data_dir=/hpc/tmp/for_Irene

# Directory to store result 
res_dir=/hpc/pmc_drost/...

qc_dir=${res_dir}/qualityCheck
mkdir -p $qc_dir

trim_dir=${res_dir}/trimmed
mkdir -p $trim_dir

align_dir=${res_dir}/alignment
mkdir -p $align_dir

figure_dir=${res_dir}/figures
mkdir -p $figure_dir

# tool dir

bowtie2Index=/hpc/pmc_drost/SOURCES/Genomes/human/bowtie2/human_gencode37_hg38

picardTool=/hpc/pmc_drost/PROJECTS/swang/software/picard.jar
new_tmp_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/tmp # to solve the out of space with the temporary output from picard

homer_dir=/hpc/pmc_drost/nhung/anaconda3/envs/cutnrun_trimgalore/bin/homer
findMotif_dir=/hpc/pmc_drost/nhung/anaconda3/envs/cutnrun_trimgalore/bin/findMotifsGenome.pl

bamCoverage_dir=/hpc/pmc_drost/nhung/anaconda3/envs/cutnrun_trimgalore/bin/bamCoverage # to merge 
effectiveGenomeSize=2913022398 
# fasta genome directory for motif finding. It has to be the same with the reference genome that was used for alignment
fasta_genome_dir=/hpc/pmc_drost/SOURCES/Genomes/human/gencode37_GRCh38_primary_assembly_genome.fa

# hg38 genes list directory to generate heatmap
hg38_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/hg38_gene_2.bed 
 
# sample_Ids: all the name and prefix as in the fastq file

sample_IDs=( "JD-103-Luc-1_AACLJGYM5_S13_L001"\
            "JD-103-Luc-2_AACLJGYM5_S14_L001"\
            "JD-103-Luc-3_AACLJGYM5_S15_L001"\
            "SCC-bulkChIC-PMC-DRO-002"\
            "SCC-bulkChIC-PMC-DRO-005"\
            "SCC-bulkChIC-PMC-DRO-008"\
            "SCC-ChIC-PMC-DRO-L5"\
            "SCC-ChIC-PMC-DRO-LH"\
            "SCC-ChIC-PMC-DRO-F1"\
            "SCC-ChIC-PMC-DRO-F5"\
            "SCC-ChIC-PMC-DRO-T1"\
            "SCC-ChIC-PMC-DRO-T5"\
            "SCC-ChIC-PMC-DRO-TH"\ 
            "SCC-ChIC-PMC-DRO-L1"\ 
            "SCC-bulkChIC-PMC-DRO-020"\
            "SCC-bulkChIC-PMC-DRO-021"\
            "SCC-bulkChIC-PMC-DRO-022"\
            "SCC-bulkChIC-PMC-DRO-023"\
            "SCC-bulkChIC-PMC-DRO-024"\
            "SCC-bulkChIC-PMC-DRO-025")
            
total_sample=${#sample_IDs[@]}

#==================================================================

#==================================================================
#                           RUNNING ANALYSIS
# turn of core file generation to not overload home dir
 ulimit -c 0 # to disable coredump  https://community.hpe.com/t5/system-administration/how-to-disable-or-restrict-core-dumps/td-p/4276097#.ZDzoIpNByDU

echo "start running cut and run analysis at $(date)"
#==================================================================

#==================================================================
#                          BEFORE PEAK CALLING
## step 1. sequencing quality check 
echo "------------------step1. running quality check--------------"
. ./1-qualityCheck.sh

# # # step 2. adapter and bad reads trimming 
echo "-------------------step 2. running trimming-----------------"
. ./2-trimming_Irene.sh 

# # # step 3. alignment- map to hg38 genome 
echo "-------------------step 3. running alignment----------------"
. ./3-alignment.sh


echo "finish cut and run analysis at $(date)"

