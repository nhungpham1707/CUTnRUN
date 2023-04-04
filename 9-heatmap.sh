#!/bin/bash

# This script is to generate heatmap visualization for transcription unit and peaks.
# ref. https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1?step=25

# Nhung Pham 27 03 2023

# transcription unit
cores=8
computeMatrix reference-point --referencePoint TSS \
  -b 1000 -a 1000 \
  -R ${data_dir}/merged_bigwig/nonDB_peaks_Luc_TFE3.bed \
-S ${merged_bigwig_dir}/fusion_merged.bw\
 ${merged_bigwig_dir}/tfe3_merged.bw\
 ${merged_bigwig_dir}/luciferase_merged.bw\
 --skipZeros -o ${merged_bigwig_dir}/matrix_gene.mat.gz -p $cores

plotHeatmap -m ${merged_bigwig_dir}/matrix_gene.mat.gz -out ${merged_bigwig_dir}/transcript.png --sortUsing sum