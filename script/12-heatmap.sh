#!/bin/bash

# This script is to generate heatmap visualization for peaks comparing with non differential binding and differential binding sites between groups.
# ref. https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1?step=25

# Nhung Pham 27 03 2023

# transcription unit
cores=8

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


computeMatrix reference-point --referencePoint TSS \
  -b 1000 -a 1000 \
  -R $refpath \
-S ${bw_IDs[@]} \
 --skipZeros -o ${figure_dir}/matrix_${savename}.mat.gz -p 8

plotHeatmap -m ${figure_dir}/matrix_${savename}.mat.gz -out ${figure_dir}/${savename}.png --sortUsing sum --zMin 0 --zMax 4

# plotHeatmap -m ${figure_dir}/matrix_${savename}.mat.gz -out ${figure_dir}/${savename}.png --sortUsing sum --heatmapHeight 50 --heatmapWidth	12 --samplesLabel Tubfus TubLuc --dpi 100 --plotFileFormat "pdf"


#  plotHeatmap -m ${figure_dir}/matrix_gene_tfe3_hg38.mat.gz -out ${figure_dir}/tfe3_kmean_hg38.png --sortUsing sum --zMin -3 --zMax 3 --kmeans 3

