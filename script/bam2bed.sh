#!/bin/bash
# convert bam file to bed and sorted to use for normalization
# Nhung 24 04 2023
# make task to allow parralel run
# task () {
#     bedtools bamtobed -i ${align_dir}/${sample_ID}/${sample_ID}.bam > ${align_dir}/${sample_ID}/${sample_ID}.bed

#     sortBed -chrThenSizeA -i ${align_dir}/${sample_ID}/${sample_ID}.bed > ${align_dir}/${sample_ID}/${sample_ID}_sorted.bed
# }

task () {
  # sort bam file before converting to bedgraph
  samtools sort -T $new_tmp_dir ${align_dir}/${sample_ID}/${sample_ID}.bam -o ${align_dir}/${sample_ID}/${sample_ID}_sorted.bam
  
  rm -f samtools.*.tmp*  # remove existing tmp files so samtools sort will not fail 

  # report in bedgraph format with 0 count bins
  # genomeCoverageBed -ibam ${align_dir}/${sample_ID}/${sample_ID}_sorted.bam -bga > ${align_dir}/${sample_ID}/${sample_ID}_sorted.bedgraph 

 # report at depth with 1 base coordinates
  genomeCoverageBed -ibam ${align_dir}/${sample_ID}/${sample_ID}_sorted.bam -d > ${align_dir}/${sample_ID}/${sample_ID}_1base_sorted.bedgraph ;
  }


n=0
for sample_ID in ${sample_IDs[@]}; do
  n=$((n+1))
  echo "-----------running $n out of $total_sample samples---------------------- "
  task "$sample_ID" &
done

wait

echo "bam2bed and sorted all done"