#!/bin/bash

# Merge all bam files from the same condition and then transform bam files to bigwig to prepare for heatmap generation

# Nhung 22 03 2023

# merge files from the same condition
mkdir -p ${merged_bigwig_dir}

savename="$1"
shift 1
sample_list=("$@") 

echo $savename
echo ${sample_list[@]}
echo ${#sample_list[@]}
total_sample=${#sample_list[@]}

merge_bw_file=${merged_bigwig_dir}/${savename}_merged.bw
if [ -s "$merge_bw_file" ] 
  then
    echo "$merge_bw_file exists. Skip merging replicates step for $savename!"
  else 
    declare -a bam_IDs=()
    n=0
        for sample_ID in ${sample_list[@]}; do
        n=$((n+1))
        echo "-----------running $n out of $total_sample samples---------------------- "

        bam_IDs+=( "${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam")

        done

        samtools merge -o ${merged_bigwig_dir}/${savename}_merged.bam ${bam_IDs[@]}

# samtools merge -o ${merged_bigwig_dir}/luciferase_merged.bam ${rm_dup_dir}/SCC-ChIC-PMC-DRO-L1/SCC-ChIC-PMC-DRO-L1_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-ChIC-PMC-DRO-L5/SCC-ChIC-PMC-DRO-L5_rmdup_filt.bam \
# ${rm_dup_dir}/bulkChIC-PMC-DRO-014/bulkChIC-PMC-DRO-014_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-bulkChIC-PMC-DRO-005/SCC-bulkChIC-PMC-DRO-005_rmdup_filt.bam

# samtools merge -o ${merged_bigwig_dir}/fusion_merged.bam ${rm_dup_dir}/SCC-ChIC-PMC-DRO-F1/SCC-ChIC-PMC-DRO-F1_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-ChIC-PMC-DRO-F5/SCC-ChIC-PMC-DRO-F5_rmdup_filt.bam \
# ${rm_dup_dir}/bulkChIC-PMC-DRO-015/bulkChIC-PMC-DRO-015_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-bulkChIC-PMC-DRO-002/SCC-bulkChIC-PMC-DRO-002_rmdup_filt.bam

##Index files
    # cd ${merged_bigwig_dir}
    # for INFILE in *.bam
    # do
    # samtools index $INFILE $INFILE".bai"
    # done

samtools sort ${merged_bigwig_dir}/${savename}_merged.bam -o ${merged_bigwig_dir}/${savename}_merged_sorted.bam   
samtools index ${merged_bigwig_dir}/${savename}_merged_sorted.bam ${merged_bigwig_dir}/${savename}_merged.bam.bai
##transform to bigwig
    $bamCoverage_dir -b ${merged_bigwig_dir}/${savename}_merged_sorted.bam \
    --normalizeUsing CPM \
    --effectiveGenomeSize 2913022398 \
    --ignoreForNormalization chrX \
    -o ${merged_bigwig_dir}/${savename}_merged.bw
fi

echo "finish merging replicates for $savename"
# $bamCoverage_dir -b luciferase_merged.bam \
#     --normalizeUsing CPM \
#     --effectiveGenomeSize 2913022398 \
#     --ignoreForNormalization chrX \
#     -o ${merged_bigwig_dir}/luciferase_merged.bw

# $bamCoverage_dir -b fusion_merged.bam \
#     --normalizeUsing CPM \
#     --effectiveGenomeSize 2913022398 \
#     --ignoreForNormalization chrX \
#     -o ${merged_bigwig_dir}/fusion_merged.bw
