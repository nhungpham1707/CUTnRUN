#!/bin/bash

# generate bedgraph files with a fix bin size for all samples to use as input for s3norm 
# nhung 25 04 2023
mkdir -p $modify_bedgraph_dir

for sample_ID in ${sample_IDs[@]}
do
bin_file=${modify_bedgraph_dir}/${sample_ID}_binsize_${bin_size}.bedgraph
if [ -s "$bin_file" ] 
    then
    echo "$bin_file exists. Skip make fix binsize bedgraph for $sample_ID!"
    else 
        $bamCoverage_dir -b ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam \
        --binSize $bin_size \
        --effectiveGenomeSize $effectiveGenomeSize \
        --numberOfProcessors 16 \
        --outFileFormat bedgraph \
        -o ${modify_bedgraph_dir}/${sample_ID}_binsize_${bin_size}.bedgraph
        fi
done