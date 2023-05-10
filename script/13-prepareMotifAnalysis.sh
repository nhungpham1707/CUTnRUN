#!/bin/bash

# This script prepare files for motif analysis:
# - merge narrow peak bed files from the same condition and extract fasta sequence that will be used as input for STREME tool meme suit

# Nhung 23 03 2023

###### First merge all NarrowPeak bedfiles

# ## tfe3
# echo "---start merging overlap peak in tfe3 at $(date)-------"
# N=4 # overlap threshold: only keep peak that is found in at least N sample
# bedops --merge ${peak_dir}/SCC-ChIC-PMC-DRO-T1/narrow/SCC-ChIC-PMC-DRO-T1_paired_control_summits.bed \
# ${peak_dir}/SCC-ChIC-PMC-DRO-T5/narrow/SCC-ChIC-PMC-DRO-T5_paired_control_summits.bed \
# ${peak_dir}/bulkChIC-PMC-DRO-016/narrow/bulkChIC-PMC-DRO-016_paired_control_summits.bed \
# ${peak_dir}/SCC-bulkChIC-PMC-DRO-008/narrow/SCC-bulkChIC-PMC-DRO-008_paired_control_summits.bed \
# > ${merged_bigwig_dir}/tfe3_summits_paired_control_merge.bed

# # get fasta seq
# cut -f 1,2,3 ${merged_bigwig_dir}/tfe3_summits_paired_control_merge.bed > ${merged_bigwig_dir}/tfe3_simple.bed
# bedtools getfasta -fi ${fasta_genome_dir} -bed ${merged_bigwig_dir}/tfe3_simple.bed -fo ${merged_bigwig_dir}/tfe3_simple_dreme.fasta

# echo "---finish merging and convert fasta from tfe3 at $(date)-------"

# # luciferase

# echo "---start merging overlap peak in luciferase at $(date)-------"
# N=4 # overlap threshold: only keep peak that is found in at least N sample
# bedops --merge ${peak_dir}/SCC-ChIC-PMC-DRO-L1/narrow/SCC-ChIC-PMC-DRO-L1_paired_control_summits.bed \
# ${peak_dir}/SCC-ChIC-PMC-DRO-L5/narrow/SCC-ChIC-PMC-DRO-L5_paired_control_summits.bed \
# ${peak_dir}/bulkChIC-PMC-DRO-014/narrow/bulkChIC-PMC-DRO-014_paired_control_summits.bed \
# ${peak_dir}/SCC-bulkChIC-PMC-DRO-005/narrow/SCC-bulkChIC-PMC-DRO-005_paired_control_summits.bed \
# > ${merged_bigwig_dir}/luciferase_summits_paired_control_merge.bed

# # get fasta seq
# cut -f 1,2,3 ${merged_bigwig_dir}/luciferase_summits_paired_control_merge.bed > ${merged_bigwig_dir}/luciferase_simple.bed
# bedtools getfasta -fi ${fasta_genome_dir} -bed ${merged_bigwig_dir}/luciferase_simple.bed -fo ${merged_bigwig_dir}/luciferase_simple_dreme.fasta

# echo "---finish merging and convert fasta from luciferase at $(date)-------"

# # fusion

# echo "---start merging overlap peak in fusion at $(date)-------"
# N=4 # overlap threshold: only keep peak that is found in at least N sample
# bedops --merge ${peak_dir}/SCC-ChIC-PMC-DRO-F1/narrow/SCC-ChIC-PMC-DRO-F1_paired_control_summits.bed \
# ${peak_dir}/SCC-ChIC-PMC-DRO-F5/narrow/SCC-ChIC-PMC-DRO-F5_paired_control_summits.bed \
# ${peak_dir}/bulkChIC-PMC-DRO-015/narrow/bulkChIC-PMC-DRO-015_paired_control_summits.bed \
# ${peak_dir}/SCC-bulkChIC-PMC-DRO-002/narrow/SCC-bulkChIC-PMC-DRO-002_paired_control_summits.bed \
# > ${merged_bigwig_dir}/fusion_summits_paired_control_merge.bed

# # get fasta seq
# cut -f 1,2,3 ${merged_bigwig_dir}/fusion_summits_paired_control_merge.bed > ${merged_bigwig_dir}/fusion_simple.bed
# bedtools getfasta -fi ${fasta_genome_dir} -bed ${merged_bigwig_dir}/fusion_simple.bed -fo ${merged_bigwig_dir}/fusion_simple_dreme.fasta

# echo "---finish merging and convert fasta from fusion at $(date)-------"

# try on normalized data
## tfe3
  save_name="$1"
  peak_dir="$2"  # ${normalize_peak_dir}
#   number_of_samples_overlap="$3" # overlap threshold: only keep peak that is found in at least N sample

  shift 2
  sample_IDs=("$@") 
  
  echo "---start merging overlap peak in $save_name at $(date)-------"

  declare -a peak_IDs=()
  total_sample=${#sample_IDs[@]}
  n=0
  for sample_ID in ${sample_IDs[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample $save_name samples---------------------- "
    peak_IDs+=( $(find ${peak_dir[@]} -maxdepth 1 -name "${sample_ID}*.narrowPeak") )

  done

bedops --merge ${peak_IDs[@]} > ${merge_dir}/savename.beds

bedops --merge ${peak_dir}/SCC-ChIC-PMC-DRO-T1/narrow/SCC-ChIC-PMC-DRO-T1_paired_control_summits.bed \
${peak_dir}/SCC-ChIC-PMC-DRO-T5/narrow/SCC-ChIC-PMC-DRO-T5_paired_control_summits.bed \
${peak_dir}/bulkChIC-PMC-DRO-016/narrow/bulkChIC-PMC-DRO-016_paired_control_summits.bed \
${peak_dir}/SCC-bulkChIC-PMC-DRO-008/narrow/SCC-bulkChIC-PMC-DRO-008_paired_control_summits.bed \
> ${merged_bigwig_dir}/tfe3_summits_paired_control_merge.bed

# get fasta seq
cut -f 1,2,3 ${merged_bigwig_dir}/tfe3_summits_paired_control_merge.bed > ${merged_bigwig_dir}/tfe3_simple.bed
bedtools getfasta -fi ${fasta_genome_dir} -bed ${merged_bigwig_dir}/tfe3_simple.bed -fo ${merged_bigwig_dir}/tfe3_simple_dreme.fasta

echo "---finish merging and convert fasta from tfe3 at $(date)-------"

# luciferase

echo "---start merging overlap peak in luciferase at $(date)-------"
N=4 # overlap threshold: only keep peak that is found in at least N sample
bedops --merge ${peak_dir}/SCC-ChIC-PMC-DRO-L1/narrow/SCC-ChIC-PMC-DRO-L1_paired_control_summits.bed \
${peak_dir}/SCC-ChIC-PMC-DRO-L5/narrow/SCC-ChIC-PMC-DRO-L5_paired_control_summits.bed \
${peak_dir}/bulkChIC-PMC-DRO-014/narrow/bulkChIC-PMC-DRO-014_paired_control_summits.bed \
${peak_dir}/SCC-bulkChIC-PMC-DRO-005/narrow/SCC-bulkChIC-PMC-DRO-005_paired_control_summits.bed \
> ${merged_bigwig_dir}/luciferase_summits_paired_control_merge.bed

# get fasta seq
cut -f 1,2,3 ${merged_bigwig_dir}/luciferase_summits_paired_control_merge.bed > ${merged_bigwig_dir}/luciferase_simple.bed
bedtools getfasta -fi ${fasta_genome_dir} -bed ${merged_bigwig_dir}/luciferase_simple.bed -fo ${merged_bigwig_dir}/luciferase_simple_dreme.fasta

echo "---finish merging and convert fasta from luciferase at $(date)-------"

# fusion

echo "---start merging overlap peak in fusion at $(date)-------"
N=4 # overlap threshold: only keep peak that is found in at least N sample
bedops --merge ${peak_dir}/SCC-ChIC-PMC-DRO-F1/narrow/SCC-ChIC-PMC-DRO-F1_paired_control_summits.bed \
${peak_dir}/SCC-ChIC-PMC-DRO-F5/narrow/SCC-ChIC-PMC-DRO-F5_paired_control_summits.bed \
${peak_dir}/bulkChIC-PMC-DRO-015/narrow/bulkChIC-PMC-DRO-015_paired_control_summits.bed \
${peak_dir}/SCC-bulkChIC-PMC-DRO-002/narrow/SCC-bulkChIC-PMC-DRO-002_paired_control_summits.bed \
> ${merged_bigwig_dir}/fusion_summits_paired_control_merge.bed

# get fasta seq
cut -f 1,2,3 ${merged_bigwig_dir}/fusion_summits_paired_control_merge.bed > ${merged_bigwig_dir}/fusion_simple.bed
bedtools getfasta -fi ${fasta_genome_dir} -bed ${merged_bigwig_dir}/fusion_simple.bed -fo ${merged_bigwig_dir}/fusion_simple_dreme.fasta

echo "---finish merging and convert fasta from fusion at $(date)-------"