#!/bin/bash
#SBATCH --job-name=peak_callingg
#SBATCH --output=peakcalling.out
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
#        -f: format of input file; Default is “AUTO” which will allow MACS to decide the format automatically.
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
# -  narrow: .narrowPeak (most important), .xls, .bed

######### Define global variables ################

sh ./variable.sh

############### steps #######################
# step 1. quality check: inspect sequencing quality with fastqc
echo "-------------------step 1. running quality check-------------------------"
sh ./1-qualityCheck.sh

# step 2. adapter and bad reads trimming 

#trimming.sh 

# step 3. alignment- map to hg38 genome 
# echo "-------------------step 3. running alignment-------------------------"
# sh ./3-alignment.sh

# step 4. filtering: remove duplciates and reads < 20bp
#echo "-------------------step 4. running filtering-------------------------"
#sh ./4-filtering.sh

# step 5. peak calling with macs2
#echo "-------------------step 5. running peak calling----------------------"
#sh ./5-peakCalling.sh

# step 5. differential peak cutnrun_analysis

#Rscript diffBind.r

# step 6. motif finding 
# step 7. super enhancer finding 
      

