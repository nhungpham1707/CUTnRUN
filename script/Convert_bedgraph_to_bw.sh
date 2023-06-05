#!/bin/bash

# ref https://www.biostars.org/p/377195/

# 2023-05-24s3norm_merge_fusion.bed is created in r (Peak_annotation_final.R) from join bedgraph files from ..sh
# $ fetchChromSizes hg38 > hg38.sizes
# /home/pmc_research/npham/PROJECTS/CUTnRUN_Maroussia/script
# cd /hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/bigwig

# # join peaks that are close to each other
# bedtools merge -d 1000 -c 4 -o sum -i 2023-05-24s3norm_merge_fusion_remove_ensembl_chr_name.bed > s3norm_merge_reps_merg_distance_remove_ensembl_chr_name_fusion.bed
# mean

# bedtools merge -d 1000 -c 4 -o mean -i 2023-05-24s3norm_merge_fusion_remove_ensembl_chr_name.bed > s3norm_merge_reps_merg_distance_remove_ensembl_chr_name_fusion_mean.bed
# # convert to bedgraph
# awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$4}' s3norm_merge_reps_merg_distance_remove_ensembl_chr_name_fusion_mean.bed > s3norm_merge_fusion_mean.bedgraph

# # LC_COLLATE=C sort -k1,1 -k2,2n s3norm_merge_fusion.bedgraph > s3norm_merge_fusion_sorted_lc_collate.bedgraph
# # sort 
# sort -k1,1 -k2,2n s3norm_merge_fusion_mean.bedgraph > s3norm_merge_fusion_sorted_mean.bedgraph

# # bedSort(s3norm_merge_fusion.bedgraph, threads = getOption("threads", 1L), sortBuffer = "1G",  sortThreads = NULL)
# bedGraphToBigWig s3norm_merge_fusion_sorted_mean.bedgraph hg38.sizes s3norm_merge_fusion_sorted_mean.bw

# head 2023-05-24s3norm_merge_fusion.bed

# sortBed -i s3norm_merge_fusion.bedgraph > s3norm_merge_fusion_sortbed.bedgraph

# # join peaks that are close to each other
input_name="$1"
result_dir="$2"
save_name="$3"
hg38_size_dir="$4"

cd $result_dir
bedtools merge -d 1000 -c 4 -o mean -i ${input_name} > ${save_name}_merge_nearby_loci_mean.bedgraph
# convert to bedgraph
# awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$4}' ${save_name}_mean.bed > ${save_name}_mean.bedgraph

# LC_COLLATE=C sort -k1,1 -k2,2n s3norm_merge_fusion.bedgraph > s3norm_merge_fusion_sorted_lc_collate.bedgraph
# sort 
sort -k1,1 -k2,2n ${save_name}_merge_nearby_loci_mean.bedgraph > ${save_name}_merge_nearby_loci_mean_sorted.bedgraph

# bedSort(s3norm_merge_fusion.bedgraph, threads = getOption("threads", 1L), sortBuffer = "1G",  sortThreads = NULL)
bedGraphToBigWig ${save_name}_merge_nearby_loci_mean_sorted.bedgraph ${hg38_size_dir}/hg38.sizes ${save_name}_merge_nearby_loci_mean_sorted.bw

