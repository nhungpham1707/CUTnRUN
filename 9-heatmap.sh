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

res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test
merged_bigwig_dir=${res_dir}/merged_bigwig

computeMatrix reference-point --referencePoint TSS \
  -b 1000 -a 1000 \
  -R hg38_gene \
-S ${merged_bigwig_dir}/fusion_merged.bw\
 ${merged_bigwig_dir}/tfe3_merged.bw\
 ${merged_bigwig_dir}/luciferase_merged.bw\
 --skipZeros -o ${merged_bigwig_dir}/matrix_gene_hg38.mat.gz -p $cores

plotHeatmap -m ${merged_bigwig_dir}/matrix_gene_hg38.mat.gz -out ${merged_bigwig_dir}/transcript_hg38.png --sortUsing sum