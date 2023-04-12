#!/bin/bash

# convert bam file to bigwig to make heatmap and visualize on IGV
# Nhung 11 04 2023




##transform to bigwig
# $bamCoverage_dir -b ${rm_dup_dir}/SCC-ChIC-PMC-DRO-FH/SCC-ChIC-PMC-DRO-FH_rmdup_filt.bam \
#     --normalizeUsing CPM \
#     --effectiveGenomeSize 2913022398 \
#     --ignoreForNormalization chrX \
#     -o ${merged_bigwig_dir}/SCC-ChIC-PMC-DRO-FH.bw

# $bamCoverage_dir -b ${rm_dup_dir}/SCC-ChIC-PMC-DRO-LH/SCC-ChIC-PMC-DRO-LH_rmdup_filt.bam \
#     --normalizeUsing CPM \
#     --effectiveGenomeSize 2913022398 \
#     --ignoreForNormalization chrX \
#     -o ${merged_bigwig_dir}/SCC-ChIC-PMC-DRO-LH.bw

# $bamCoverage_dir -b ${rm_dup_dir}/SCC-ChIC-PMC-DRO-TH/SCC-ChIC-PMC-DRO-TH_rmdup_filt.bam \
#     --normalizeUsing CPM \
#     --effectiveGenomeSize 2913022398 \
#     --ignoreForNormalization chrX \
#     -o ${merged_bigwig_dir}/SCC-ChIC-PMC-DRO-TH.bw


task () {
    bam_ID=( $(find ${rm_dup_dir}/${sample_ID} -name "*_rmdup_filt.bam") )
    echo $bam_ID
    $bamCoverage_dir -b $bam_ID \
        --normalizeUsing CPM \
        --effectiveGenomeSize 2913022398 \
        --ignoreForNormalization chrX \
        -o ${bigwig_dir}/$sample_ID.bw ;
}
# run loop

n=0
total_sample=${h3k4me3[@]}
for sample_ID in ${h3k4me3[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample samples---------------------- "

    task '$sample_ID' &

done
