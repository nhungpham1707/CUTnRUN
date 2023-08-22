#!/bin/bash
# call peak from bedgraph file after data normalization 
# Nhung 02 05 2023
mkdir -p $peak_after_s3norm_dir
bedgraph_dir=$modify_bedgraph_dir/$s3norm_sample_file_name/S3norm_NBP_bedgraph
output_peak_dir=$peak_after_s3norm_dir
save_name=$s3norm_sample_file_name

echo $bedgraph_dir
echo $output_peak_dir
echo $save_name
 task(){
macs2 bdgpeakcall -i $bedgraph_dir/${sample_ID}_*.bedgraph -o ${output_peak_dir}/${sample_ID}_$save_name.narrowPeak ;

 } 

 n=0

for sample_ID in ${sample_IDs[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample samples---------------------- "
    peak_file=${output_peak_dir}/${sample_ID}_$save_name.narrowPeak
    if [ -s "$peak_file" ] 
        then
        echo "$peak_file exists. Skip peak calling for $sample_ID!"
    else 
        task '$sample_ID' &
    fi
done

wait

echo "all done peak calling after normalization"