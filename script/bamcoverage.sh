#!/bin/bash
# nhung 24 04 2023
# generate counts for bin of different sizes
# task () {
#     cd ${res_dir} # to redirect tmp files there
#     $bamCoverage_dir -b ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam \
#         --binSize $bin_size \
#         --effectiveGenomeSize $effectiveGenomeSize \
#         --numberOfProcessors 8 \
#         -o ${bincount_dir}/${sample_ID}_binsize_${bin_size}.bedgraph ;

# }



###########
# https://www.biostars.org/p/295903/
# Generate 50bp windows using bedtools makewindows function
# Convert bedGraph to bigwig using bedGraphToBigwig
# Use bed file from step-1 and bigwig from step-2 with multiBigwigSummary from deeptools or bigWigAverageOverBed
# https://groups.google.com/g/bedtools-discuss/c/r64iIYA5tOM
# download hg38 genome
# wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

# rename file
# mv hg38.chrom.sizes hg38.genome
task () {
multiBigwigSummary BED-file -b ${bigwig_dir}/${sample_ID}.bw -o ${bincount_dir}/${sample_ID}_${bin_size}.npz --BED ${bincount_dir}/hg38_${bin_size}.bed --outRawCounts ${bincount_dir}/${sample_ID}_${bin_size}_scores_per_transcript.tab ;
}

#run loop

total_sample=${#sample_IDs[@]}
# for bin_size in 50 100 150 200 250 300 350 400 450 500 550 600 750 800 850 900 950 1000
# for bin_size in $(seq 50 50 1000)
for bin_size in $(seq 50 50 1000)
do
    echo $bin_size
    n=0
    bedtools makewindows -g /home/pmc_research/npham/PROJECTS/CUTnRUN_Maroussia/hg38.genome -w ${bin_size} > ${bincount_dir}/hg38_${bin_size}.bed
    for sample_ID in ${sample_IDs[@]}; do
        n=$((n+1))
        echo "-----------running $n out of $total_sample samples for bin size $bin_size---------------------- "

        task '$sample_ID' &

    done
    echo "all done for bin size $bin_size"
done 


wait
echo "all done bin count"

