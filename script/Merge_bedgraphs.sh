#!/bin/bash
# this script combines bedgraph files from the same group of sample
# Nhung 24 05 2023

### Merge bedgraph file 
bg_savename="$1"
bg_path="$2"
res_bedgraph_dir="$3"
shift 3
sample_list=("$@")
echo "variables in Merge_bedgraphs.sh are "
echo $bg_savename
echo $bg_path
echo $res_bedgraph_dir
echo ${sample_list[@]}
# get file path
declare -a sample_paths=()
total_sample=${#sample_list[@]}
n=0
for sample_ID in ${sample_list[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample samples---------------------- "

    sample_paths+=( "${bg_path}/${sample_ID}*.bedgraph" )

done
bedtools unionbedg -i ${sample_paths[@]} > ${res_bedgraph_dir}/${bg_savename}_union.bedgraph