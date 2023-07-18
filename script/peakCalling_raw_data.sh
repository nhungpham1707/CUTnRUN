#!/bin/bash


# this script is to call peak per group with the approriate control and everything together without control 
# Nhung, 21 03 2023

# Tool: 
# - macs2 ver 2.2.7.1


# Input file options

# -t: The IP data file (this is the only REQUIRED parameter for MACS)
# -c: The control or mock data file
# -f: format of input file; Default is “AUTO” which will allow MACS to decide the format automatically.
# -g: mappable genome size which is defined as the genome size which can be sequenced; some precompiled values provided.
# Output arguments
# --outdir: MACS2 will save all output files into speficied folder for this option
# -n: The prefix string for output files

# Input before data normalization: 
# - per group: 1 bam file for the test, 1 bam file for the control (after removing duplicate and filtering)
# Input after normalization:
# - bedgraph files in S3norm_NBP_bedgraph 

# Output:
# -  broad: .broadPeak (most important), .gappedPeak, .xls 
# -  narrow: .narrowPeak (most important), .xls, .bed

mkdir -p ${peak_dir}

##### Run per group using appropriate controls 
task () {
  output_dir=${peak_dir}/${sample_ID}
  mkdir -p ${output_dir}

  sample_dir=( $(find ${rm_dup_dir}/${sample_ID} -name "*rmdup_filt.bam") )
  
  echo "----------------start-narrow-peak-calling-for $sample_ID at $(date)----------"
  macs2 callpeak --tempdir $new_tmp_dir -t ${sample_dir} -c ${control} -f BAMPE -g hs -n ${sample_ID}_paired_control -q 0.01 --outdir $output_dir/narrow

  echo "----------------start-broad-peak-calling-for $sample_ID at $(date)----------"

  macs2 callpeak --tempdir $new_tmp_dir -t ${sample_dir} -c ${control} -f BAMPE -g hs -n ${sample_ID}_paired_control -q 0.01 --broad --outdir $output_dir/broad

  echo "----------------finish-broad-peak-calling-for $sample_ID at $(date)----------" ;

}

control="$1"
shift 1
sample_list=("$@") 

total_sample=${#sample_list[@]}
n=0

for sample_ID in $sample_list; do
  n=$((n+1))
  echo "-----------running $n out of $total_sample samples---------------"
    peak_file=( $(find ${peak_dir}/${sample_ID}  -name "*.narrowPeak") )
  if [ -s "$peak_file" ] 
  then
    echo "$peak_file exists. Skip peakcalling w control step for $sample_ID!"
  else 
  task "$sample_ID" &
  fi
done
wait

# generate report

#############################################
multiqc_file=$peak_dir/peak_control_report.html

if [ -s "$multiqc_file" ] 
  then
    echo "$multiqc_file exists. Skip generating qc report for peakcalling w control step!"
  else 
multiqc ${peak_dir} -n peak_control_report -o ${peak_dir}
fi

echo "all done for peak calling with control"