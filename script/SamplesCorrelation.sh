#!/bin/bash
# this script calculate the read coverages for genomic regions and generate correlation plot between samples. The idea is to check which sample is correlated and if replicates correlate with each other
# more detail here https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html and https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html

# automatically generate bam file path

declare -a bam_ID=()
for sample_ID in ${sample_IDs[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample samples---------------------- "

    # bam_ID+=( "${bigwig_dir}/${sample_ID}.bw" )
    bam_ID+=( "${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam")

done

multiBamSummary bins --bamfiles ${bam_ID[@]} -o $figure_dir/bamsummary_results.npz

plotCorrelation -in $figure_dir/bamsummary_results.npz \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o $figure_dir/heatmap_SpearmanCorr_readCounts.png   \
    --outFileCorMatrix $figure_dir/SpearmanCorr_readCounts.tab

plotCorrelation -in $figure_dir/bamsummary_results.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o $figure_dir/heatmap_PearsonCorr_readCounts.png   \
    --outFileCorMatrix $figure_dir/PearsonCorr_readCounts.tab