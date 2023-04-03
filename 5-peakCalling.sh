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

# Input: 
# - per group: 1 bam file for the test, 1 bam file for the control (after removing duplicate and filtering)
# 
# Output:
# -  broad: .broadPeak (most important), .gappedPeak, .xls 
# -  narrow: .narrowPeak (most important), .xls, .bed

# Run per group using appropriate controls 
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

# run for tfe3
total_sample=${#tfe3[@]}
n=0

for sample_ID in $tfe3; do
  n=$((n+1))
  echo "-----------running $n out of $total_sample tfe3 samples---------------"
  
  control=$tfe3C

  task "$sample_ID" &
  
done
wait
echo "all done for tfe3"

# run for luciferase
total_sample=${#luciferase[@]}
n=0
for sample_ID in $luciferase;do
  n=$((n+1))
  echo "-----------running $n out of $total_sample luciferase samples---------------"

  control=$luciferaseC

  task "$sample_ID" &

done

wait

echo "all done for luciferase"

# run for fusion

total_sample=${#fusion[@]}
n=0
for sample_ID in $fusion;do
  n=$((n+1))
  echo "-----------running $n out of $total_sample fusion samples---------------"

  control=$fusionC

  task "$sample_ID" &

done

wait

echo "all done for fusion"

# Run for all samples without control

task () {
  output_dir=${peak_no_control_dir}/${sample_ID}
  mkdir -p ${output_dir}

  sample_dir=( $(find ${rm_dup_dir}/${sample_ID} -name "*rmdup_filt.bam") )
  
  echo "----------------start-narrow-peak-calling-no-control-for $sample_ID at $(date)----------"
  macs2 callpeak --tempdir $new_tmp_dir -t ${sample_dir} -f BAMPE -g hs -n ${sample_ID}_paired -q 0.01 --outdir $output_dir/narrow

  echo "----------------start-broad-peak-calling-for $sample_ID at $(date)----------"

  macs2 callpeak --tempdir $new_tmp_dir -t ${sample_dir} -f BAMPE -g hs -n ${sample_ID}_paired -q 0.01 --broad --outdir $output_dir/broad

  echo "----------------finish-broad-peak-calling-no-control-for $sample_ID at $(date)----------" ;

}


total_sample=${#allT[@]}
n=0
for sample_ID in $allT;do
  n=$((n+1))
  echo "-----------running $n out of $total_sample all samples---------------"

  control=$fusionC

  task "$sample_ID" &

done

wait

echo "all done for no control"


