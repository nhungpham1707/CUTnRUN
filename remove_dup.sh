#!/bin/bash
#SBATCH --job-name=rmdup_test
#SBATCH --output=rmdup_test.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=90G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl


# this script is to sort and remove duplicate

# Nhung, 14 03 2023

# Tool: 
# - picard ver
# - samtools, ver 1.16.1

# Input: 1 bam file for each sample (output from bowtie2)

# Output:
# -  1 bam file with duplicate remove

#########
# data dir
data_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/SCC_ChIC-PMC-DRO_plates_20210520_run1


res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test
picardTool=/hpc/pmc_drost/PROJECTS/swang/software/picard.jar

# sample_Ids was generated as text file from ls filename.txt from the data_dir
sample_IDs=( "bulkChIC-PMC-DRO-011" \
      "bulkChIC-PMC-DRO-013"\
      "bulkChIC-PMC-DRO-016"\
      "SCC-bulkChIC-PMC-DRO-005"\
      )


task () {
	sample_ID=$str
  echo $sample_ID
    
  bam_ID=( $(find ${res_dir}/alignment_parallel/${sample_ID} -name "*.bam") )
  
  len=${#fastq_IDs[@]}
  echo $len
  mkdir -p ${res_dir}/rm_dup/${sample_ID}
  sample_rmdup=${res_dir}/rm_dup/${sample_ID}
  
  echo "start sorting and remove dup for $sample_ID at $(date)"
# sort
  java -jar $picardTool SortSam \
	I=${bam_ID} \
	O=${sample_rmdup}/${sample_ID}_sorted.bam \
	VALIDATION_STRINGENCY=LENIENT \
	SORT_ORDER=coordinate
  	  echo "start mark dup for $sample_ID at $(date)"

# Mark duplicates
java -jar $picardTool MarkDuplicates \
	I=${sample_rmdup}/${sample_ID}_sorted.bam  \
	O=${sample_rmdup}/${sample_ID}_mkdup.bam \
	VALIDATION_STRINGENCY=LENIENT \
	M=${sample_rmdup}/${sample_ID}/marked_dup_metrics.txt
  echo "start remove dup for $sample_ID at $(date)"

# Remove duplicates
java -jar $picardTool SortSam \
	I=${sample_rmdup}/${sample_ID}_sorted.bam\
	O=${sample_rmdup}/${sample_ID}_rmdup.bam \
	SORT_ORDER=coordinate
  echo "start sorting with samtools for $sample_ID at $(date)"

# ##Filter with samtools # 1024 for pcr or optical dups https://wikis.utexas.edu/display/CoreNGSTools/Filtering+with+
# http://www.htslib.org/doc/samtools-view.html
# -b to get output in bam file
# -q 20 skip alignment smaller than 20
# -F do not output alignment with any bit set in flag (1024 mean from pcr or optical dups)
# -h keep header
# -f only output alignment with all bits set in flag (2 = read map in proper pair https://www.samformat.info/sam-format-flag)
 samtools view -F 1024 -f 2 -q 20 -h -b ${sample_rmdup}/${sample_ID}_rmdup.bam -o ${sample_rmdup}/${sample_ID}_rmdup_filt.bam
# ##Index with samtools
# -@ n thread
  echo "start indexing for $sample_ID at $(date)"

 samtools index -@ 16 ${sample_rmdup}/${sample_ID}_rmdup_filt.bam ${sample_rmdup}/${sample_ID}_rmdup_filt.bam.bai

  echo "finish indexing for $sample_ID at $(date)" ;
}

for str in ${sample_IDs[@]}; do
  
  task '$str' & 

done

wait 
echo "all done"












