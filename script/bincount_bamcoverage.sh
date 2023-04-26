#!/bin/bash
# nhung 24 04 2023
# generate counts for bin of different sizes
task () {
    cd ${res_dir} # to redirect tmp files there
    $bamCoverage_dir -b ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam \
        --binSize $bin_size \
        --effectiveGenomeSize $effectiveGenomeSize \
        --numberOfProcessors 32 \
        -o ${bincount_dir}/${sample_ID}_binsize_${bin_size}.bedgraph ;

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

        task '$sample_ID' 

    done
    echo "all done for bin size $bin_size"
done 


wait
echo "all done bin count"

