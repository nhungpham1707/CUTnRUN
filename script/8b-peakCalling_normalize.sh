#!/bin/bash
# call peak from bedgraph file after data normalization 
# Nhung 02 05 2023

bedgraph_dir="$1"
output_peak_dir="$2"
save_name="$3"
shift 3
sample_list=("$@") 

echo $bedgraph_dir
echo $output_peak_dir
echo $save_name
echo ${sample_list[@]}
 task(){
macs2 bdgpeakcall -i $bedgraph_dir/${sample_ID}_*.bedgraph -o ${output_peak_dir}/${sample_ID}_$save_name.narrowPeak ;

 } 

 n=0

for sample_ID in ${sample_list[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample samples---------------------- "

    task '$sample_ID' &

done

wait


#  task(){
# macs2 bdgpeakcall -i $clean_normalize_dir/${sample_ID}_*.bedgraph -o ${normalize_peak_dir}/${sample_ID}_normalize.narrowPeak ;

#  } 

#  n=0

# for sample_ID in ${sample_IDs[@]}; do
#     n=$((n+1))
#     echo "-----------running $n out of $total_sample samples---------------------- "

#     task '$sample_ID' &

# done

# wait


echo "all done peak calling after normalization"