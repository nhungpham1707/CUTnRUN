#!/bin/bash

# convert bam file to bigwig to make heatmap and visualize on IGV
# Nhung 11 04 2023
mkdir -p $bigwig_dir

task () {
    $bamCoverage_dir -b ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam \
        --normalizeUsing CPM \
        --effectiveGenomeSize $effectiveGenomeSize \
        --ignoreForNormalization chrX \
        -o ${bigwig_dir}/${sample_ID}.bw ;
}

#run loop

n=0
total_sample=${#sample_IDs[@]}
for sample_ID in ${sample_IDs[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample samples---------------------- "
    bw_file=$bigwig_dir/$sample_ID.bw
    if [ -s "$bw_file" ] 
    then
    echo "$bw_file exists. Skip bam2bigwig for $sample_ID!"
    else 
    task '$sample_ID' &
    fi

done

wait
echo "all done bam2bigwig"
