#!/bin/bash

# this script is to conduct fastqc prior to trimming for all samples. The code is designed to run many samples in parallel

# Nhung, 21 03 2023

# Tool: fastqc v0.12.1

# Input: R1 and R2 fastq sequences for each sample

# Output:
# - output from each sample will be saved in a subfolder with the respective sample ID
# - For each sample a zip file and a html file


task () {

  fastq_IDs=( $(find ${data_dir}/$sample_ID -name "*.fastq.gz") )
  #declare -p fastq_IDs

  len=${#fastq_IDs[@]}
  echo $len
  sample_fastqc=${qc_dir}/${sample_ID}
  mkdir -p $sample_fastqc
  
  echo "---------start fastqc for $sample_ID at $(date)---------"

  fastqc ${fastq_IDs[0]} -o $sample_fastqc
  fastqc ${fastq_IDs[1]} -o $sample_fastqc

  echo "----------finish fastqc for $sample_ID at $(date)---------" ;
}

n=0
for sample_ID in ${sample_IDs[@]}; do
  n=$((n+1))
  echo "-----------running $n out of $total_sample samples---------------------- "
  task "$sample_ID" &
done

wait

# generate report
multiqc ${qc_dir}
echo "all done quality check"
