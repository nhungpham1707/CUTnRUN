#!/bin/bash

# This script is to generate heatmap visualization for peaks comparing with non differential binding and differential binding sites between groups.
# ref. https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1?step=25

# Nhung Pham 27 03 2023

# transcription unit
cores=8
# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R ${data_dir}/merged_bigwig/nonDB_peaks_Luc_TFE3.bed \
# -S ${merged_bigwig_dir}/fusion_merged.bw\
#  ${merged_bigwig_dir}/tfe3_merged.bw\
#  ${merged_bigwig_dir}/luciferase_merged.bw\
#  --skipZeros -o ${merged_bigwig_dir}/matrix_gene.mat.gz -p $cores

# plotHeatmap -m ${merged_bigwig_dir}/matrix_gene.mat.gz -out ${merged_bigwig_dir}/transcript.png --sortUsing sum


#cores=8
# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R ${data_dir}/merged_bigwig/nonDB_peaks_Luc_TFE3.bed \
# -S ${merged_bigwig_dir}/fusion_merged.bw\
#   --skipZeros -o ${merged_bigwig_dir}/matrix_gene_fusion.mat.gz -p $cores

# plotHeatmap -m ${merged_bigwig_dir}/matrix_gene_fusion.mat.gz -out ${merged_bigwig_dir}/fusion_transcript.png --sortUsing sum

# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R ${data_dir}/merged_bigwig/nonDB_peaks_Luc_TFE3.bed \
# -S ${merged_bigwig_dir}/tfe3_merged.bw\
#   --skipZeros -o ${merged_bigwig_dir}/matrix_gene_tfe3.mat.gz -p $cores

# plotHeatmap -m ${merged_bigwig_dir}/matrix_gene_tfe3.mat.gz -out ${merged_bigwig_dir}/tfe3_transcript.png --sortUsing sum

# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R ${data_dir}/merged_bigwig/nonDB_peaks_Luc_TFE3.bed \
# -S ${merged_bigwig_dir}/luciferase_merged.bw\
#   --skipZeros -o ${merged_bigwig_dir}/matrix_gene_luc.mat.gz -p $cores

# plotHeatmap -m ${merged_bigwig_dir}/matrix_gene_luc.mat.gz -out ${merged_bigwig_dir}/luc_transcript.png --sortUsing sum

# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R ${data_dir}/merged_bigwig/DBregions_fusionTFE3Luc_noHistones.bed \
# -S ${merged_bigwig_dir}/luciferase_merged.bw\
#   --skipZeros -o ${merged_bigwig_dir}/matrix_gene_luc_db.mat.gz -p $cores

# plotHeatmap -m ${merged_bigwig_dir}/matrix_gene_luc_db.mat.gz -out ${merged_bigwig_dir}/luc_db_transcript.png --sortUsing sum

# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R ${data_dir}/merged_bigwig/DBregions_fusionTFE3Luc_noHistones.bed \
# -S ${merged_bigwig_dir}/fusion_merged.bw\
#   --skipZeros -o ${merged_bigwig_dir}/matrix_gene_fusion_db.mat.gz -p $cores

# plotHeatmap -m ${merged_bigwig_dir}/matrix_gene_fusion_db.mat.gz -out ${merged_bigwig_dir}/fusion_db_transcript.png --sortUsing sum

# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R ${data_dir}/merged_bigwig/DBregions_fusionTFE3Luc_noHistones.bed \
# -S ${merged_bigwig_dir}/tfe3_merged.bw\
#   --skipZeros -o ${merged_bigwig_dir}/matrix_gene_tfe3_db.mat.gz -p $cores

# plotHeatmap -m ${merged_bigwig_dir}/matrix_gene_tfe3_db.mat.gz -out ${merged_bigwig_dir}/tfe3_db_transcript.png --sortUsing sum

# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R ${data_dir}/merged_bigwig/DBregions_fusionTFE3Luc_noHistones.bed \
# -S ${merged_bigwig_dir}/fusion_merged.bw\
#  ${merged_bigwig_dir}/tfe3_merged.bw\
#  ${merged_bigwig_dir}/luciferase_merged.bw\
#  --skipZeros -o ${merged_bigwig_dir}/matrix_gene_db.mat.gz -p $cores

# plotHeatmap -m ${merged_bigwig_dir}/matrix_gene_db.mat.gz -out ${merged_bigwig_dir}/transcript_db.png --sortUsing sum
cores=8

## plot against hg38 gene list

# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R ${hg38_dir} \
# -S ${merged_bigwig_dir}/fusion_merged.bw\
#  ${merged_bigwig_dir}/tfe3_merged.bw\
#  ${merged_bigwig_dir}/luciferase_merged.bw\
#  --skipZeros -o ${merged_bigwig_dir}/matrix_gene_hg38.mat.gz -p $cores

# plotHeatmap -m ${merged_bigwig_dir}/matrix_gene_hg38.mat.gz -out ${merged_bigwig_dir}/transcript_hg38.png --sortUsing sum

# plot for histone samples
# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R ${hg38_dir} \
# -S ${bigwig_dir}/SCC-ChIC-PMC-DRO-FH.bw\
#  ${bigwig_dir}/SCC-ChIC-PMC-DRO-LH.bw\
#  ${bigwig_dir}/SCC-ChIC-PMC-DRO-TH.bw\
#  --skipZeros -o ${bigwig_dir}/matrix_gene_histone_hg38.mat.gz -p $cores

# plotHeatmap -m ${bigwig_dir}/matrix_gene_histone_hg38.mat.gz -out ${bigwig_dir}/histone_kmean_hg38.png --sortUsing sum --zMin -3 --zMax 3 --kmeans 3

# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R ${hg38_dir} \
# -S ${bigwig_dir}/bulkChIC-PMC-DRO-011.bw\
#  ${bigwig_dir}/bulkChIC-PMC-DRO-012.bw\
#  ${bigwig_dir}/bulkChIC-PMC-DRO-013.bw\
#  --skipZeros -o ${bigwig_dir}/matrix_gene_histone_2_hg38.mat.gz -p $cores

# plotHeatmap -m ${bigwig_dir}/matrix_gene_histone_2_hg38.mat.gz -out ${bigwig_dir}/histone_2_kmean_hg38.png --sortUsing sum --zMin -3 --zMax 3 --kmeans 3

# plot heatmap for replicates in the same group

#heatmap function 
# input: 1. save name of heatmap (i.e. tfe3vshg38), 2. dir path to the reference genome and 3. array of sample IDs to plot heatmap
# ref https://askubuntu.com/questions/674333/how-to-pass-an-array-as-function-argument

# heatmap () {
#   local savename="$1"
#   local refpath="$2"
#   shift 2
#   local sample_IDs=("$@") 
  
#   declare -a bw_IDs=()
#   total_sample=${#sample_IDs[@]}
#   n=0
#   for sample_ID in ${sample_IDs[@]}; do
#     n=$((n+1))
#     echo "-----------running $n out of $total_sample $savename samples---------------------- "

#     bw_IDs+=( "${bigwig_dir}/${sample_ID}.bw" )

#   done
# echo ${bw_IDs[@]}
# # computeMatrix reference-point --referencePoint TSS \
# #   -b 1000 -a 1000 \
# #   -R $refpath \
# # -S ${bw_IDs[@]} \
# #  --skipZeros -o ${figure_dir}/matrix_${savename}.mat.gz -p 8

# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R $refpath \
# -S ${bw_IDs[@]} \
#  --skipZeros -o ${figure_dir}/matrix_${savename}.mat.gz -p 8

# plotHeatmap -m ${figure_dir}/matrix_${savename}.mat.gz -out ${figure_dir}/${savename}.png --sortUsing sum

# }

# heatmap "tfe3vshg38" "$hg38_dir" "${tfe3[@]}"

  savename="$1"
  refpath="$2"
  sample_dir="$3"
  shift 3
  sample_IDs=("$@") 
  
  declare -a bw_IDs=()
  total_sample=${#sample_IDs[@]}
  n=0
  for sample_ID in ${sample_IDs[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample $savename samples---------------------- "

    bw_IDs+=( "${sample_dir}/${sample_ID}.bw" )

  done
echo ${bw_IDs[@]}
# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R $refpath \
# -S ${bw_IDs[@]} \
#  --skipZeros -o ${figure_dir}/matrix_${savename}.mat.gz -p 8

computeMatrix reference-point --referencePoint TSS \
  -b 1000 -a 1000 \
  -R $refpath \
-S ${bw_IDs[@]} \
 --skipZeros -o ${figure_dir}/matrix_${savename}.mat.gz -p 8

# plotHeatmap -m ${figure_dir}/matrix_${savename}.mat.gz -out ${figure_dir}/${savename}.png --sortUsing sum --zMin 0 --zMax 4

plotHeatmap -m ${figure_dir}/matrix_${savename}.mat.gz -out ${figure_dir}/${savename}.png --sortUsing sum --heatmapHeight 50 --heatmapWidth	12 --samplesLabel Tubfus TubLuc --dpi 100 --plotFileFormat "pdf"

# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R ${hg38_dir} \
# -S ${bigwig_dir}/SCC-ChIC-PMC-DRO-T1.bw\
#  ${bigwig_dir}/SCC-ChIC-PMC-DRO-T5.bw\
#  ${bigwig_dir}/bulkChIC-PMC-DRO-016.bw\
#  ${bigwig_dir}/SCC-bulkChIC-PMC-DRO-008.bw\
#  --skipZeros -o ${figure_dir}/matrix_gene_tfe3_hg38.mat.gz -p $cores

#  plotHeatmap -m ${figure_dir}/matrix_gene_tfe3_hg38.mat.gz -out ${figure_dir}/tfe3_kmean_hg38.png --sortUsing sum --zMin -3 --zMax 3 --kmeans 3

# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R ${hg38_dir} \
# -S ${bigwig_dir}/SCC-ChIC-PMC-DRO-F1.bw\
#  ${bigwig_dir}/SCC-ChIC-PMC-DRO-F5.bw\
#  ${bigwig_dir}/bulkChIC-PMC-DRO-015.bw\
#  ${bigwig_dir}/SCC-bulkChIC-PMC-DRO-002.bw\
#  --skipZeros -o ${figure_dir}/matrix_gene_fusion_hg38.mat.gz -p $cores

# plotHeatmap -m ${figure_dir}/matrix_gene_fusion_hg38.mat.gz -out ${figure_dir}/fusion_kmean_hg38.png --sortUsing sum --zMin -3 --zMax 3 --kmeans 3
# luciferase
# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R ${hg38_dir} \
# -S ${bigwig_dir}/SCC-ChIC-PMC-DRO-L1.bw\
#  ${bigwig_dir}/SCC-ChIC-PMC-DRO-L5.bw\
#  ${bigwig_dir}/bulkChIC-PMC-DRO-014.bw\
#  ${bigwig_dir}/SCC-bulkChIC-PMC-DRO-005.bw\
#  --skipZeros -o ${figure_dir}/matrix_gene_luc_hg38.mat.gz -p $cores

# plotHeatmap -m ${figure_dir}/matrix_gene_luc_hg38.mat.gz -out ${figure_dir}/luc_kmean_hg38.png --sortUsing sum --zMin -3 --zMax 3 --kmeans 3

# new fusion samples

# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R ${hg38_dir} \
# -S ${bigwig_dir}/SCC-bulkChIC-PMC-DRO-023.bw\
#  ${bigwig_dir}/SCC-bulkChIC-PMC-DRO-024.bw\
#  ${bigwig_dir}/SCC-bulkChIC-PMC-DRO-025.bw\
#  --skipZeros -o ${figure_dir}/matrix_gene_new_fusion_hg38.mat.gz -p $cores

# plotHeatmap -m ${figure_dir}/matrix_gene_new_fusion_hg38.mat.gz -out ${figure_dir}/new_fusion_3_scale_hg38.png --sortUsing sum --zMin -3 -3 -3 --zMax 3 3 3
# plotHeatmap -m ${figure_dir}/matrix_gene_new_fusion_hg38.mat.gz -out ${figure_dir}/new_fusion_04_scale_hg38.png --sortUsing sum --zMin -0.4 -0.4 -0.4 --zMax 0.4 0.4 0.4
# new luc sample
# computeMatrix reference-point --referencePoint TSS \
#   -b 1000 -a 1000 \
#   -R ${hg38_dir} \
# -S ${bigwig_dir}/SCC-bulkChIC-PMC-DRO-020.bw\
#  ${bigwig_dir}/SCC-bulkChIC-PMC-DRO-021.bw\
#  ${bigwig_dir}/SCC-bulkChIC-PMC-DRO-022.bw\
#  --skipZeros -o ${figure_dir}/matrix_gene_new_luc_hg38.mat.gz -p $cores

# plotHeatmap -m ${figure_dir}/matrix_gene_new_luc_hg38.mat.gz -out ${figure_dir}/new_luc_3_scale_hg38.png --sortUsing sum --zMin -3 --zMax 3 
# plotHeatmap -m ${figure_dir}/matrix_gene_new_luc_hg38.mat.gz -out ${figure_dir}/new_luc_04_scale_hg38.png --sortUsing sum --zMin -0.4 --zMax 0.4 