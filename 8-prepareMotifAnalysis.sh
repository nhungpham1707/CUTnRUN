#!/bin/bash

# This script prepare files for motif analysis
# Merge all bed file from the same condition and use bedops to keep only unique overlapping regions

# Nhung 23 03 2023

###### First merge all NarrowPeak bedfiles

# tfe3
#bedtools merge -o ${merged_bigwig}/tfe3_summits_merge.bed ${peakCalling}/SCC-ChIC-PMC-DRO-T1/narrow/SCC-ChIC-PMC-DRO-T1_paired_control_summits.bed \
#${peakCalling}/SCC-ChIC-PMC-DRO-T5/narrow/SCC-ChIC-PMC-DRO-T5_paired_control_summits.bed \
#${peakCalling}/bulkChIC-PMC-DRO-016/narrow/bulkChIC-PMC-DRO-016_paired_control_summits.bed \
#${peakCalling}/SCC-bulkChIC-PMC-DRO-008/narrow/SCC-bulkChIC-PMC-DRO-008_paired_control_summits.bed
echo "---start merging overlap peak in tfe3 at $(date)-------"
N=4 # overlap threshold: only keep peak that is found in at least N sample
echo $peak_dir
echo $merged_bigwig_dir

#bedops --intersect ${peak_dir}/SCC-ChIC-PMC-DRO-T1/narrow/SCC-ChIC-PMC-DRO-T1_paired_control_summits.bed \
#${peak_dir}/SCC-ChIC-PMC-DRO-T5/narrow/SCC-ChIC-PMC-DRO-T5_paired_control_summits.bed \
#${peak_dir}/bulkChIC-PMC-DRO-016/narrow/bulkChIC-PMC-DRO-016_paired_control_summits.bed \
#${peak_dir}/SCC-bulkChIC-PMC-DRO-008/narrow/SCC-bulkChIC-PMC-DRO-008_paired_control_summits.bed \
#> ${merged_bigwig_dir}/tfe3_summits_paired_control_intervals_overlapping.bed

bedops --merge ${peak_dir}/SCC-ChIC-PMC-DRO-T1/narrow/SCC-ChIC-PMC-DRO-T1_paired_control_summits.bed \
${peak_dir}/SCC-ChIC-PMC-DRO-T5/narrow/SCC-ChIC-PMC-DRO-T5_paired_control_summits.bed \
${peak_dir}/bulkChIC-PMC-DRO-016/narrow/bulkChIC-PMC-DRO-016_paired_control_summits.bed \
${peak_dir}/SCC-bulkChIC-PMC-DRO-008/narrow/SCC-bulkChIC-PMC-DRO-008_paired_control_summits.bed \
> ${merged_bigwig_dir}/tfe3_summits_paired_control_merge.bed

cut -f 1,2,3 ${merged_bigwig_dir}/tfe3_summits_paired_control_merge.bed > ${merged_bigwig_dir}/tfe3_simple.bed
bedtools getfasta -fi ${fasta_genome_dir} -bed ${merged_bigwig_dir}/tfe3_simple.bed -fo ${merged_bigwig_dir}/tfe3_simple_dreme.fasta



#bedops --everything ${peak_dir}/SCC-ChIC-PMC-DRO-T1/narrow/SCC-ChIC-PMC-DRO-T1_paired_control_summits.bed \
#${peak_dir}/SCC-ChIC-PMC-DRO-T5/narrow/SCC-ChIC-PMC-DRO-T5_paired_control_summits.bed \
#${peak_dir}/bulkChIC-PMC-DRO-016/narrow/bulkChIC-PMC-DRO-016_paired_control_summits.bed \
#${peak_dir}/SCC-bulkChIC-PMC-DRO-008/narrow/SCC-bulkChIC-PMC-DRO-008_paired_control_summits.bed \
#| bedmap --echo --echo-map-id-uniq --delim '\t' - \
#| cut -f1-3,5 - \
#| awk -vN=${N} '{ n = split($4, a, ";"); if (n==N) { print $0; } }' \
#> ${merged_bigwig_dir}/tfe3_summits_paired_control_intervals_overlapping.bed
echo "---finish merging overlap peak in tfe3 at $(date)-------"

#ref. https://www.biostars.org/p/236044/
#bedops -m ${peakCalling}/SCC-ChIC-PMC-DRO-T1/narrow/SCC-ChIC-PMC-DRO-T1_paired_control_summits.bed \
#${peakCalling}/SCC-ChIC-PMC-DRO-T5/narrow/SCC-ChIC-PMC-DRO-T5_paired_control_summits.bed \
#${peakCalling}/bulkChIC-PMC-DRO-016/narrow/bulkChIC-PMC-DRO-016_paired_control_summits.bed \
#${peakCalling}/SCC-bulkChIC-PMC-DRO-008/narrow/SCC-bulkChIC-PMC-DRO-008_paired_control_summits.bed > ${merged_bigwig}/tfe3_summits_paired_control_merge.bed

# samtools merge -o ${merged_bigwig}/luciferase_merged.bam ${rm_dup_dir}/SCC-ChIC-PMC-DRO-L1/SCC-ChIC-PMC-DRO-L1_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-ChIC-PMC-DRO-L5/SCC-ChIC-PMC-DRO-L5_rmdup_filt.bam \
# ${rm_dup_dir}/bulkChIC-PMC-DRO-014/bulkChIC-PMC-DRO-014_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-bulkChIC-PMC-DRO-005/SCC-bulkChIC-PMC-DRO-005_rmdup_filt.bam

# samtools merge -o ${merged_bigwig}/fusion_merged.bam ${rm_dup_dir}/SCC-ChIC-PMC-DRO-F1/SCC-ChIC-PMC-DRO-F1_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-ChIC-PMC-DRO-F5/SCC-ChIC-PMC-DRO-F5_rmdup_filt.bam \
# ${rm_dup_dir}/bulkChIC-PMC-DRO-015/bulkChIC-PMC-DRO-015_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-bulkChIC-PMC-DRO-002/SCC-bulkChIC-PMC-DRO-002_rmdup_filt.bam


##Index files
#cd ${merged_bigwig}
# for INFILE in *.bam
# do
#    samtools index $INFILE $INFILE".bai"
# done

