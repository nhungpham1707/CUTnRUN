#!/bin/bash
#SBATCH --job-name=check_rmdup
#SBATCH --output=20-03-rmdup_test_error.out
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
# - picard ver
# - samtools, ver 1.16.1

# Input: 1 bam file for each sample (output from bowtie2)

# Output:
# -  1 bam file with duplicate remove

######### define global variables ################

data_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/SCC_ChIC-PMC-DRO_plates_20210520_run1
# result dir
res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test
rm_dup_dir=$res_dir/rm_dup_test

# tool dir
picardTool=/hpc/pmc_drost/PROJECTS/swang/software/picard.jar
new_tmp_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/tmp # to solve the out of space with the temporary output from picard


# sample_Ids was generated as text file from ls filename.txt from the data_dir

sample_ID=( "bulkChIC-PMC-DRO-011")
      
########### Make task function ##############
task () {
	sample_ID=$str
  echo $sample_ID
    
  bam_ID=( $(find ${res_dir}/alignment/${sample_ID} -name "*.bam") )
  
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
	

  	  echo "finish picard sort for $sample_ID at $(date)" ;
}

for str in ${sample_IDs[@]}; do
  
  task '$str' & 

done

wait 
echo "all done"

