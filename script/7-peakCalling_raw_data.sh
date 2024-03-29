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

# run for tfe3
# total_sample=${#tfe3[@]}
# n=0

# for sample_ID in $tfe3; do
#   n=$((n+1))
#   echo "-----------running $n out of $total_sample tfe3 samples---------------"
  
#   control=$tfe3C

#   task "$sample_ID" &
  
# done
# wait
# echo "all done for tfe3"

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

#### Run for all samples without control

# task () {
#   output_dir=${peak_no_control_dir}/${sample_ID}
#   mkdir -p ${output_dir}

#   sample_dir=( $(find ${rm_dup_dir}/${sample_ID} -name "*rmdup_filt.bam") )
  
#   echo "----------------start-narrow-peak-calling-no-control-for $sample_ID at $(date)----------"
#   macs2 callpeak --tempdir $new_tmp_dir -t ${sample_dir} -f BAMPE -g hs -n ${sample_ID}_paired -q 0.01 --outdir $output_dir/narrow

#   echo "----------------start-broad-peak-calling-for $sample_ID at $(date)----------"

#   macs2 callpeak --tempdir $new_tmp_dir -t ${sample_dir} -f BAMPE -g hs -n ${sample_ID}_paired -q 0.01 --broad --outdir $output_dir/broad

#   echo "----------------finish-broad-peak-calling-no-control-for $sample_ID at $(date)----------" ;

# }

# task () {
#   output_dir=${peak_no_control_dir}/${sample_ID}
#   mkdir -p ${output_dir}

#   sample_dir=( $(find ${rm_dup_dir}/${sample_ID} -name "*rmdup_filt.bam") )
  
#   echo "----------------start-narrow-peak-calling-no-control-for $sample_ID at $(date)----------"
#   macs2 callpeak -t ${sample_dir} -f BAMPE -g hs -n ${sample_ID}_paired -q 0.01 --outdir $output_dir/narrow

#   echo "----------------start-broad-peak-calling-for $sample_ID at $(date)----------"

#   macs2 callpeak -t ${sample_dir} -f BAMPE -g hs -n ${sample_ID}_paired -q 0.01 --broad --outdir $output_dir/broad

#   echo "----------------finish-broad-peak-calling-no-control-for $sample_ID at $(date)----------" ;

# }


# total_sample=${#allT[@]}
# total_sample=${#sample_IDs[@]}

# n=0
# for sample_ID in ${sample_IDs[@]};do
#   n=$((n+1))
#   echo "-----------running $n out of $total_sample all samples---------------"

#   #control=$fusionC

#   task "$sample_ID" &

# done

# wait

# echo "all done for no control"

# generate report

#############################################
## run with normalized data

task () {
  output_dir=${normalize_peak_dir}/${sample_ID}
  mkdir -p ${output_dir}

  sample_dir=( ${s3norm_working_directory}/S3norm_NBP_bedgraph/${sample_ID}_nozeroes_augmented.bedgraph.s3norm.bedgraph)
  
  echo "----------------start-narrow-peak-calling-for $sample_ID at $(date)----------"
  macs2 callpeak --tempdir $new_tmp_dir -t ${sample_dir} -c ${control} -f BAMPE -g hs -n ${sample_ID}_paired_control -q 0.01 --outdir $output_dir/narrow

  echo "----------------start-broad-peak-calling-for $sample_ID at $(date)----------"

  macs2 callpeak --tempdir $new_tmp_dir -t ${sample_dir} -c ${control} -f BAMPE -g hs -n ${sample_ID}_paired_control -q 0.01 --broad --outdir $output_dir/broad

  echo "----------------finish-broad-peak-calling-for $sample_ID at $(date)----------" ;

}

# run for tfe3
# total_sample=${#tfe3[@]}
# n=0

# for sample_ID in $tfe3; do
#   n=$((n+1))
#   echo "-----------running $n out of $total_sample tfe3 samples---------------"
  
#   control=$tfe3C

#   task "$sample_ID" &
  
# done
# wait
# echo "all done for tfe3"

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




multiqc ${peak_dir} -n peak_control_report -o ${peak_dir}
multiqc ${peak_no_control_dir} -n peak_no_control_report -o ${peak_no_control_dir}

echo ("all done for peak calling")