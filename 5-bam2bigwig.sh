#!/bin/bash

# convert bam file to bigwig to make heatmap and visualize on IGV
# Nhung 11 04 2023




##transform to bigwig
# $bamCoverage_dir -b ${rm_dup_dir}/SCC-ChIC-PMC-DRO-FH/SCC-ChIC-PMC-DRO-FH_rmdup_filt.bam \
#     --normalizeUsing CPM \
#     --effectiveGenomeSize 2913022398 \
#     --ignoreForNormalization chrX \
#     -o ${bigwig_dir}/SCC-ChIC-PMC-DRO-FH.bw

# $bamCoverage_dir -b ${rm_dup_dir}/SCC-ChIC-PMC-DRO-LH/SCC-ChIC-PMC-DRO-LH_rmdup_filt.bam \
#     --normalizeUsing CPM \
#     --effectiveGenomeSize 2913022398 \
#     --ignoreForNormalization chrX \
#     -o ${bigwig_dir}/SCC-ChIC-PMC-DRO-LH.bw

# $bamCoverage_dir -b ${rm_dup_dir}/SCC-ChIC-PMC-DRO-TH/SCC-ChIC-PMC-DRO-TH_rmdup_filt.bam \
#     --normalizeUsing CPM \
#     --effectiveGenomeSize 2913022398 \
#     --ignoreForNormalization chrX \
#     -o ${bigwig_dir}/SCC-ChIC-PMC-DRO-TH.bw

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

    task '$sample_ID' &

done

wait
echo "all done bam2bigwig"
