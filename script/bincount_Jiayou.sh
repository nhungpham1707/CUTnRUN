#!/bin/bash

# generate bedgraph files with a fix bin size for all samples to use as input for s3norm 
# nhung 25 04 2023

bin_size=200
sample_dir=
output_dir=
sample_list=

for sample_ID in ${sample_list[@]}
do
$bamCoverage_dir -b ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam \
        --binSize $bin_size \
        --effectiveGenomeSize $effectiveGenomeSize \
        --numberOfProcessors 16 \
        --outFileFormat bedgraph \
        -o ${bincount_dir}/${sample_ID}_binsize_${bin_size}.bedgraph
done