#!/bin/bash
# call peak from bedgraph file after data normalization 
# Nhung 02 05 2023
 task(){
#  macs2 bdgpeakcall -i ${s3norm_working_directory}/S3norm_NBP_bedgraph/${sample_ID}_nozeroes_augmented.bedgraph.NBP.s3norm.bedgraph -o ${normalize_peak_dir}/${sample_ID}_normalize.narrowPeak ;
macs2 bdgpeakcall -i ${s3norm_working_directory}/S3norm_NBP_bedgraph/${sample_ID}_*.bedgraph -o ${normalize_peak_dir}/${sample_ID}_normalize.narrowPeak ;
 } 

 n=0

for sample_ID in ${sample_IDs[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample samples---------------------- "

    task '$sample_ID' &

done

wait


echo "all done peak calling after normalization"