#!/bin/bash

# This script is to find motif from macs2 narrow peak 
# Nhung 22 03 2023


task () {

out_dir=${motif_dir}/${sample_ID}
mkdir -p ${out_dir}

sample_dir=( $(find ${peak_dir}/${sample_ID}/narrow -name "*paired_control_summits.bed") ) 

echo "----------------start-motif-finding-for $sample_ID at $(date)----------"

$findMotif_dir $sample_dir hg38 $_out_dir

echo "----------------finish-motif-finding-for $sample_ID at $(date)----------" ;
}


total_sample=${#allT[@]}
n=0
for sample_ID in $allT; do

n=$((n+1))
  echo "-----------running $n out of $total_sample  samples---------------"

  task "$sample_ID" &

done

wait

echo "all done"
