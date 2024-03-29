#!/bin/bash
# The fraction of reads in peaks (commonly denoted as FRiP) is a measure of how many of the total reads are found within the peaks. It is calculated as the number of reads in peaks divided by the total number of reads. This value can be used to assess the quality of the data and the specificity of the experiment. A high fraction of reads in peaks indicates that the majority of the reads are located in the regions of interest, and that the experiment has a high signal-to-noise ratio. On the other hand, a low fraction of reads in peaks may indicate that the majority of the reads are located in non-specific regions, and that the experiment has a low specificity. ref https://pluto.bio/blog/understanding-fraction-and-reads-in-peaks
# ref https://www.biostars.org/p/337872/
# Nhung 02 05 2023
# conda env cutnrun_trimgalore
mkdir -p $frip_dir
save_name=raw_data

frip_file=$frip_dir/${save_name}_read_in_peak.txt
if [ -s "$frip_file" ] 
then
    echo "$frip_file exists. Skip calculating frip step!"
else 
    n=0
    read_in_peak_list=()
    total_reads_list=()
    for sample_ID in ${sample_IDs[@]};do
    n=$((n+1))
    echo "-----------running $n out of $total_sample samples---------------"

    reads_in_peaks=$(bedtools sort -i ${peak_no_control_dir}/${sample_ID}/narrow/*.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam -b stdin -ubam | samtools view -c)
    read_in_peak_list=(${read_in_peak_list[@]} $reads_in_peaks)

    total_reads=$(samtools view -c ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam) 
    total_reads_list=(${total_reads_list[@]} $total_reads) 

    done

    echo ${total_reads_list[@]} > ${frip_dir}/${save_name}_total_reads.txt 
    echo ${read_in_peak_list[@]} > ${frip_dir}/${save_name}_read_in_peak.txt
    echo ${sample_IDs[@]} > ${frip_dir}/${save_name}_sample_IDs.txt

    # calculate and plot frip in R
    export FIGURE_DIR_VARIABLE=$figure_dir
    export FRIP_DIR_VARIABLE=$frip_dir
    export SAVE_NAME_VARIABLE=$save_name
    Rscript plot_FRiP.R

fi

echo "all done for FRiP calculation"