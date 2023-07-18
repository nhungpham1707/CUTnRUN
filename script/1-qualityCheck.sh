#!/bin/bash

# this script is to conduct fastqc prior to trimming for all samples. The code is designed to run many samples in parallel

# Nhung, 21 03 2023

# Tool: fastqc v0.12.1

# Input: R1 and R2 fastq sequences for each sample

# Output:
# - output from each sample will be saved in a subfolder with the respective sample ID
# - For each sample a zip file and a html file
mkdir -p $qc_dir

task () {

  fastq_IDs=( $(find ${data_dir}/$sample_ID -maxdepth 1 -name "*.$input_extension") )
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

qc_file=${qc_dir}/quality_check_report.html 

if [ -s "$qc_file" ] 
then
    echo "$qc_file exists. Skip qc step!"
else 
  n=0
  for sample_ID in ${sample_IDs[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample samples---------------------- "
    task "$sample_ID" &
  done

  wait
  # generate report
  multiqc ${qc_dir} -n quality_check_report -o ${qc_dir}
  echo "all done quality check"
fi