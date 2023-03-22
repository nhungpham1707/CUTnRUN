#!/bin/bash
#SBATCH --job-name=motif
#SBATCH --output=motif_findingg.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=90G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl


# This is a master script to analyze cut and run data 
# Nhung, 20 03 2023

################## Tools ###########################
# Quality check:
# - fastqc v0.12.1

# Trimming:
# - trim_galore ver 0.6.10

# Alignment: 
#	- bowtie ver 2.5.1, input flags:
#			-end-to-end: for trimmed sequence
#			-p 16: 16 threads
#			-samtools -bS: to save input sam file as bam file

# Remove dup and filtering for short reads:
# - picard ver
# - samtools, ver 1.16.1, input flags:
# 	(https://wikis.utexas.edu/display/CoreNGSTools/Filtering+with+)
# 	(http://www.htslib.org/doc/samtools-view.html)
# 		-b to get output in bam file
# 		-q 20 skip alignment smaller than 20
# 		-F do not output alignment with any bit set in flag (1024 mean from pcr or optical dups)
# 		-h keep header
# 		-f only output alignment with all bits set in flag (2 = read map in proper pair https://www.samformat.info/sam-format-flag)

# Peak calling
# - macs2 ver 2.2.7.1, input flags:
#        -t: The IP data file (this is the only REQUIRED parameter for MACS)
#        -c: The control or mock data file
#        -f: format of input file; Default is “AUTO” which will allow MACS to decide the format automatically. In BAMPE mode, there is no need to estimate fragment lengths since the actual insertion length for each read pair will be considered. That's why there is no model.R generated. https://github.com/macs3-project/MACS/issues/281
#        -g: mappable genome size which is defined as the genome size which can be sequenced; some precompiled values provided.
#       Output arguments
#       --outdir: MACS2 will save all output files into speficied folder for this option
#       -n: The prefix string for output files

################### Input and output data ############################

# Quality check 
#   Input: R1 and R2 fastq sequences for each sample
#   Output:
#     - output from each sample will be saved in a subfolder with the respective sample ID
#     - For each sample a zip file and a html file

# Trimming 
#   Input: R1 and R2 fastq sequences for each sample
#   Output:
#       - 2 trimming reports for R1 and R2.
#       - 2 val_ files for R1 and R2 (inputs for next step)
#       - 2 fastqc html for R1 and R2
#       - 2 fastqc.zip files for R1 and R2.


# Alignment
#   Input: R1 and R2 outputs from trimming for each sample
#   Output:
#       - 1 bam file per sample stored in its respective folder 

# Remove dup and filtering for short reads: 
#		Input data: 1 bam file for each sample (output from alignment)
# 	Output:
# 		- 1 mark and remove duplicate bam file
# 		- 1 filter bam file (input for next step peakcalling)
# 		- 1 index bam.bai file
# 		- 1 metric text file

# Peak calling
#   Input: 
#       - per group: 1 bam file for the test, 1 bam file for the control (after removing duplicate and filtering)
#   Output:
# -  broad: .broadPeak (most important), .gappedPeak, .xls 
# -  narrow: .narrowPeak (most important), .xls, summits.bed 
# summits.bed has information of a peak where the score (protein binding intensity) is maximum in the peak whereas the narrowPeak file has complete peak boundary. The summits.bed reports a peak with width/distance 1bp. That means this is the region (of peak) where the protein binding intensity rich it's maximum. Therefore every peak should have its own summit. ref. https://www.biostars.org/p/9470414/
# more about summit and extraction of seq can be found here https://notebook.community/ssjunnebo/pathogen-informatics-training/Notebooks/ChIP-Seq/motif-analysis
######### Define global variables ################

data_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/SCC_ChIC-PMC-DRO_plates_20210520_run1

# result dir
res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test

qc_dir=${res_dir}/qualityCheck
mkdir -p $qc_dir

trim_dir=${res_dir}/trimmed
mkdir -p $trim_dir

align_dir=${res_dir}/alignment
mkdir -p $align_dir

rm_dup_dir=$res_dir/rm_dup
mkdir -p $rm_dup_dir

peak_dir=${res_dir}/peakCalling
mkdir -p ${peak_dir}

peak_no_control_dir=${res_dir}/peakCalling_nocontrol
mkdir -p ${peak_no_control_dir}

motif_dir= ${res_dir}/motif
mkdir -p ${motif_dir}
# tool dir


bowtie2Index=/hpc/pmc_drost/SOURCES/Genomes/human/bowtie2/human_gencode37_hg38

picardTool=/hpc/pmc_drost/PROJECTS/swang/software/picard.jar
new_tmp_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/tmp # to solve the out of space with the temporary output from picard

homer_dir=/hpc/pmc_drost/nhung/anaconda3/envs/cutnrun_trimgalore/bin/homer
findMotif_dir=/hpc/pmc_drost/nhung/anaconda3/envs/cutnrun_trimgalore/bin/findMotifsGenome.pl

# sample_Ids was generated as text file from ls filename.txt from the data_dir

sample_IDs=( "bulkChIC-PMC-DRO-011" \
            "bulkChIC-PMC-DRO-012"\
            "bulkChIC-PMC-DRO-013"\
            "bulkChIC-PMC-DRO-014"\
            "bulkChIC-PMC-DRO-015"\
            "bulkChIC-PMC-DRO-016"\
            "SCC-bulkChIC-PMC-DRO-002"\
            "SCC-bulkChIC-PMC-DRO-005"\
            "SCC-bulkChIC-PMC-DRO-008"\
            "SCC-ChIC-PMC-DRO-L5"\
            "SCC-ChIC-PMC-DRO-LH"\
            "SCC-ChIC-PMC-DRO-F1"\
            "SCC-ChIC-PMC-DRO-F5"\
            "SCC-ChIC-PMC-DRO-FH"\
            "SCC-ChIC-PMC-DRO-T1"\
            "SCC-ChIC-PMC-DRO-T5"\
            "SCC-ChIC-PMC-DRO-TH"\ 
            "SCC-ChIC-PMC-DRO-L1")

total_sample=${#sample_IDs[@]}

# classify sample for peakcalling
tfe3=("SCC-ChIC-PMC-DRO-T1" "SCC-ChIC-PMC-DRO-T5" "bulkChIC-PMC-DRO-016" "SCC-bulkChIC-PMC-DRO-008")
luciferase=("SCC-ChIC-PMC-DRO-L1 SCC-ChIC-PMC-DRO-L5" "bulkChIC-PMC-DRO-014" "SCC-bulkChIC-PMC-DRO-005")
fusion=("SCC-ChIC-PMC-DRO-F1" "SCC-ChIC-PMC-DRO-F5" "bulkChIC-PMC-DRO-015" "SCC-bulkChIC-PMC-DRO-002")
#allT="$tfe3 $luciferase $fusion"
allT=( "SCC-ChIC-PMC-DRO-F1")
tfe3C=$res_dir/rm_dup/bulkChIC-PMC-DRO-013/bulkChIC-PMC-DRO-013_rmdup_filt.bam
luciferaseC=$res_dir/rm_dup/bulkChIC-PMC-DRO-011/bulkChIC-PMC-DRO-011_rmdup_filt.bam
fusionC=$res_dir/rm_dup/bulkChIC-PMC-DRO-012/bulkChIC-PMC-DRO-012_rmdup_filt.bam

############### steps #######################
# step 1. quality check: inspect sequencing quality with fastqc
#echo "------------------step1. running quality check----------------------"
#. ./1-qualityCheck.sh

# step 2. adapter and bad reads trimming 
#echo "-------------------step 2. running trimming--------------------------"
#. ./2-trimming.sh 

# step 3. alignment- map to hg38 genome 
#echo "-------------------step 3. running alignment-------------------------"
#. ./3-alignment.sh

# step 4. filtering: remove duplciates and reads < 20bp
# echo "-------------------step 4. running filtering-------------------------"
# . ./4-filtering.sh

# step 5. peak calling with macs2
#echo "-------------------step 5. running peak calling----------------------"
#. ./5-peakCalling.sh

# step 6. motif finding 
echo "-------------------step 6. running motif finding----------------------"
. ./6-motifFinding.sh 

# step 6. differential peak cutnrun_analysis

#Rscript diffBind.r

# step 7. motif finding 
# step 7. super enhancer finding 
