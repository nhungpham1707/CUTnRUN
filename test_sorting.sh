#!/bin/bash
#SBATCH --job-name=sorting_test
#SBATCH --output=test_sorting.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=90G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl


# this script is to sort and remove duplicate
# test again why the remove_dup resulted in small files
# Nhung, 17 03 2023

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

rm_dup_dir=$res_dir/rm_dup_test
mkdir -p $rm_dup_dir
# tool dir


picardTool=/hpc/pmc_drost/PROJECTS/swang/software/picard.jar
new_tmp_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/tmp # to solve the out of space with the temporary output from picard

# sample_Ids was generated as text file from ls filename.txt from the data_dir

sample_ID=( "bulkChIC-PMC-DRO-011")
      
########### Make task function ##############

  echo $sample_ID
    
  bam_ID=( $(find ${res_dir}/alignment_test/${sample_ID} -name "*.bam") )
  
  len=${#bam_ID[@]}
  echo $len
  sample_rmdup=$rm_dup_dir/${sample_ID}
  mkdir -p $sample_rmdup

  echo "start sorting and remove dup for $sample_ID at $(date)"
# sort
  java -Xmx2g -Djava.io.tmpdir=$new_tmp_dir -jar $picardTool SortSam \
	-I ${bam_ID} \
	-O ${sample_rmdup}/${sample_ID}_sorted.bam \
	-VALIDATION_STRINGENCY LENIENT \
	-SORT_ORDER coordinate \
	-TMP_DIR $new_tmp_dir
	  echo "finsh sorting $sample_ID at $(date)"
