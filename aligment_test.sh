#!/bin/bash
#SBATCH --job-name=check_alignment
#SBATCH --output=alignment_test.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=90G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl


# this script is to test alignment
# test again why the remove_dup resulted in small files
# Nhung, 20 03 2023

# Tool: 
# - bowtie ver , input flags:
#			-end-to-end: for trimmed sequence
#			-p 16: 16 threads
#			-samtools -bS: to save input sam file as bam file
# - picard ver
# - samtools, ver 1.16.1, input flags:
# 	(https://wikis.utexas.edu/display/CoreNGSTools/Filtering+with+)
# 	(http://www.htslib.org/doc/samtools-view.html)
# 		-b to get output in bam file
# 		-q 20 skip alignment smaller than 20
# 		-F do not output alignment with any bit set in flag (1024 mean from pcr or optical dups)
# 		-h keep header
# 		-f only output alignment with all bits set in flag (2 = read map in proper pair https://www.samformat.info/sam-format-flag)

# Input data: 1 bam file for each sample (output from alignment w bowtie2)

# Output:
# -  1 bam file with duplicate remove

######### define global variables ################

data_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/SCC_ChIC-PMC-DRO_plates_20210520_run1

# result dir
res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test
align_dir=${res_dir}/alignment_test
mkdir -p $align_dir
#rm_dup_dir=$res_dir/rm_dup_test
#mkdir -p $rm_dup_dir
# tool dir

#genomelib=/hpc/pmc_gen/references/RNA-Seq/starfusion/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/
#gtf=/hpc/pmc_drost/SOURCES/Genomes/human/gencode37_GRCh38_annotation.gtf
bowtie2Index=/hpc/pmc_drost/SOURCES/Genomes/human/bowtie2/human_gencode37_hg38

#picardTool=/hpc/pmc_drost/PROJECTS/swang/software/picard.jar
#new_tmp_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/tmp # to solve the out of space with the temporary output from picard

# sample_Ids was generated as text file from ls filename.txt from the data_dir

sample_ID=( "bulkChIC-PMC-DRO-011")
      
########### Make task function ##############
# alignment 
trim_IDs=( $(find ${res_dir}/trimmed/${sample_ID} -name "*.fq.gz") )

declare -p trim_IDs

len=${#trim_IDs[@]}
  
echo $len

sample_align_dir=${align_dir}/${sample_ID}
mkdir -p $sample_align_dir

 echo "-----start alignment for $sample_ID at $(date)-------"

# Alignment to both target 

bowtie2 --end-to-end -p 16 -x $bowtie2Index -1 ${trim_IDs[0]} -2 ${trim_IDs[1]} | samtools view -bS - > ${sample_align_dir}/${sample_ID}.bam

echo "-------finish alignment for $sample_ID at $(date)--------"
     
