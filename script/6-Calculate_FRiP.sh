#!/bin/bash
# The fraction of reads in peaks (commonly denoted as FRiP) is a measure of how many of the total reads are found within the peaks. It is calculated as the number of reads in peaks divided by the total number of reads. This value can be used to assess the quality of the data and the specificity of the experiment. A high fraction of reads in peaks indicates that the majority of the reads are located in the regions of interest, and that the experiment has a high signal-to-noise ratio. On the other hand, a low fraction of reads in peaks may indicate that the majority of the reads are located in non-specific regions, and that the experiment has a low specificity. ref https://pluto.bio/blog/understanding-fraction-and-reads-in-peaks
# ref https://www.biostars.org/p/337872/
# Nhung 02 05 2023
# conda env cutnrun_trimgalore
#conda install -c bioconda subread 

# res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test
# peak_no_control_dir=${res_dir}/peakCalling_nocontrol
# frip_dir=${res_dir}/FRiP
# mkdir -p $frip_dir
# rm_dup_dir=$res_dir/rm_dup

# sample_ID="bulkChIC-PMC-DRO-011"
# sample_IDs=( "bulkChIC-PMC-DRO-011" "SCC-ChIC-PMC-DRO-T1")

# task() {
# reads_in_peaks=$(bedtools sort -i ${peak_no_control_dir}/${sample_ID}/narrow/*.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam -b stdin -ubam | samtools view -c)
# # 5407614

# total_reads=$(samtools view -c ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam) ;
# # 5407614/18278174
# }

# total_sample=${#sample_IDs[@]}
# n=0
# read_in_peak_list=()
# total_reads_list=()
# for sample_ID in ${sample_IDs[@]};do
#   n=$((n+1))
#   echo "-----------running $n out of $total_sample samples---------------"

#    task "$sample_ID" 
#     read_in_peak_list=(${read_in_peak_list[@]} $reads_in_peaks)
#     total_reads_list=(${total_reads_list[@]} $total_reads) &
# done

# wait

# echo "all done for FRiP calculation"

# echo ${total_reads[@]} > ${frip_dir}/total_reads.txt 
# echo ${read_in_peak_list[@]} > ${frip_dir}/read_in_peak.txt
# echo ${sample_IDs[@]} > ${frip_dir}/sample_IDs.txt


res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test
peak_no_control_dir=${res_dir}/peakCalling_nocontrol
frip_dir=${res_dir}/FRiP
mkdir -p $frip_dir
rm_dup_dir=$res_dir/rm_dup

task() {
# reads_in_peaks=$(bedtools sort -i ${peak_no_control_dir}/${sample_ID}/narrow/*.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam -b stdin -ubam | samtools view -c)
# 5407614

total_reads=$(samtools view -c ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam) ;
# 5407614/18278174
}
total_sample=${#sample_IDs[@]}
n=0
# read_in_peak_list=()
total_reads_list=()
for sample_ID in ${sample_IDs[@]};do
  n=$((n+1))
  echo "-----------running $n out of $total_sample samples---------------"
    # read_in_peak_list=(${read_in_peak_list[@]} $reads_in_peaks)
    total_reads=$(samtools view -c ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam) ;
    total_reads_list=(${total_reads_list[@]} $total_reads) 
done

echo ${total_reads_list[@]} > ${frip_dir}/total_reads.txt 
# echo ${read_in_peak_list[@]} > ${frip_dir}/read_in_peak.txt
# echo ${sample_IDs[@]} > ${frip_dir}/sample_IDs.txt
echo "all done for FRiP calculation"



