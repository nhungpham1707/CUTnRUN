#!/bin/bash
# The fraction of reads in peaks (commonly denoted as FRiP) is a measure of how many of the total reads are found within the peaks. It is calculated as the number of reads in peaks divided by the total number of reads. This value can be used to assess the quality of the data and the specificity of the experiment. A high fraction of reads in peaks indicates that the majority of the reads are located in the regions of interest, and that the experiment has a high signal-to-noise ratio. On the other hand, a low fraction of reads in peaks may indicate that the majority of the reads are located in non-specific regions, and that the experiment has a low specificity. ref https://pluto.bio/blog/understanding-fraction-and-reads-in-peaks
# ref https://www.biostars.org/p/337872/
# Nhung 02 05 2023
# conda env cutnrun_trimgalore

# total_sample=${#sample_IDs[@]}
# n=0
# read_in_peak_list=()
# total_reads_list=()
# for sample_ID in ${sample_IDs[@]};do
#   n=$((n+1))
#   echo "-----------running $n out of $total_sample samples---------------"

#   reads_in_peaks=$(bedtools sort -i ${peak_no_control_dir}/${sample_ID}/narrow/*.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam -b stdin -ubam | samtools view -c)
#   read_in_peak_list=(${read_in_peak_list[@]} $reads_in_peaks)

#   total_reads=$(samtools view -c ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam) 
#   total_reads_list=(${total_reads_list[@]} $total_reads) 

# done

# echo ${total_reads_list[@]} > ${frip_dir}/total_reads.txt 
# echo ${read_in_peak_list[@]} > ${frip_dir}/read_in_peak.txt
# echo ${sample_IDs[@]} > ${frip_dir}/sample_IDs.txt



## after normalization

# normalize_sample_count_NBP=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/modify_bedgraph/S3norm_NBP_bedgraph/
# normalize_peak_file=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/peak_s3norm_antiTfe3_new_control

# s3norm_sample_file=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/modify_bedgraph/s3norm_antiTFE3_Samples/anti_TFE3_new_control.csv
# sample_list=($(awk -F "\"*\\t\"*" '{print $1}' $s3norm_sample_file)) 
# sample_list=`awk -F "\"*\\t\"*" '{print $1}' /hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/modify_bedgraph/s3norm_antiTFE3_Samples/anti_TFE3_new_control.csv`
# sample_count_dir="$1"
# peak_dir="$2"
# save_name="$3"
# shift 3
# sample_list=("$@")
# echo $peak_dir
# echo $sample_count_dir
# echo $save_name
# echo ${sample_list[1]}
# echo ${sample_list[2]}
# sample_count_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/modify_bedgraph/s3norm_antiTFE3_Samples/with_new_control/S3norm_rc_bedgraph
# peak_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/peak_s3norm_antiTfe3_new_control
# sample_list=( "bulkChIC-PMC-DRO-014"\
#             "bulkChIC-PMC-DRO-015"\
#             "bulkChIC-PMC-DRO-016"\
#             "SCC-bulkChIC-PMC-DRO-002"\
#             "SCC-bulkChIC-PMC-DRO-005"\
#             "SCC-bulkChIC-PMC-DRO-008"\
#             "SCC-ChIC-PMC-DRO-L5"\
#             "SCC-ChIC-PMC-DRO-F1"\
#             "SCC-ChIC-PMC-DRO-F5"\
#             "SCC-ChIC-PMC-DRO-T1"\
#             "SCC-ChIC-PMC-DRO-T5"\
#             "SCC-ChIC-PMC-DRO-L1" )
# save_name=normalize_anti_TFE3_new_control

# total_sample=${#sample_list[@]}
# n=0
# read_in_peak_list=()
# total_reads_list=()
# for sample_ID in ${sample_list[@]};do
#   n=$((n+1))
#   echo "-----------running $n out of $total_sample samples---------------"

# #   reads_in_peaks=$(bedtools sort -i $normalize_peak_file/${sample_ID}*.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a $normalize_sample_count/${sample_ID}*.bedgraph -b stdin -ubam | samtools view -c)
#   reads_in_peaks=$(awk '{ sum += $5; } END { print sum; }' "$@" ${peak_dir}/${sample_ID}*.narrowPeak )
#   read_in_peak_list=(${read_in_peak_list[@]} $reads_in_peaks)

#   total_reads=$(awk '{ sum += $4; } END { print sum; }' "$@" ${sample_count_dir}/${sample_ID}*.bedgraph)
#   total_reads_list=(${total_reads_list[@]} $total_reads) 

# #   total_reads_nbp= $(awk '{ sum += $4; } END { print sum; }' "$@" $normalize_sample_count/${sample_ID}*.bedgraph) 
# done
# echo ${total_reads_list[@]} > ${frip_dir}/${save_name}_total_reads.txt 
# echo ${read_in_peak_list[@]} > ${frip_dir}/${save_name}_read_in_peak.txt
# echo ${sample_list[@]} > ${frip_dir}/${save_name}_sample_IDs.txt


# # calculate and plot frip in R
# export FIGURE_DIR_VARIABLE=$figure_dir
# export FRIP_DIR_VARIABLE=$frip_dir
# export SAVE_NAME_VARIABLE=$save_name
# Rscript plot_FRiP.R
# echo "all done for FRiP calculation"



#!/bin/bash
# The fraction of reads in peaks (commonly denoted as FRiP) is a measure of how many of the total reads are found within the peaks. It is calculated as the number of reads in peaks divided by the total number of reads. This value can be used to assess the quality of the data and the specificity of the experiment. A high fraction of reads in peaks indicates that the majority of the reads are located in the regions of interest, and that the experiment has a high signal-to-noise ratio. On the other hand, a low fraction of reads in peaks may indicate that the majority of the reads are located in non-specific regions, and that the experiment has a low specificity. ref https://pluto.bio/blog/understanding-fraction-and-reads-in-peaks
# ref https://www.biostars.org/p/337872/
# Nhung 02 05 2023
# conda env cutnrun_trimgalore

# total_sample=${#sample_IDs[@]}
# n=0
# read_in_peak_list=()
# total_reads_list=()
# for sample_ID in ${sample_IDs[@]};do
#   n=$((n+1))
#   echo "-----------running $n out of $total_sample samples---------------"

#   reads_in_peaks=$(bedtools sort -i ${peak_no_control_dir}/${sample_ID}/narrow/*.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam -b stdin -ubam | samtools view -c)
#   read_in_peak_list=(${read_in_peak_list[@]} $reads_in_peaks)

#   total_reads=$(samtools view -c ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam) 
#   total_reads_list=(${total_reads_list[@]} $total_reads) 

# done

# echo ${total_reads_list[@]} > ${frip_dir}/total_reads.txt 
# echo ${read_in_peak_list[@]} > ${frip_dir}/read_in_peak.txt
# echo ${sample_IDs[@]} > ${frip_dir}/sample_IDs.txt



## after normalization

# sample_count_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/modify_bedgraph/s3norm_antiTFE3_Samples/with_new_control/S3norm_rc_bedgraph
# peak_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/peak_s3norm_antiTfe3_new_control

# save_name=normalize_anti_TFE3_new_control

normalize_sample_count=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/modify_bedgraph/s3norm_antiTFE3_Samples/with_new_control/S3norm_rc_bedgraph
normalize_sample_count=${res_dir}/modify_bedgraph/s3norm_antiTFE3_Samples
# normalize_sample_count_NBP=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/modify_bedgraph/S3norm_NBP_bedgraph/
# normalize_peak_file=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/peak_s3norm_antiTfe3_new_control
# normalize_peak_calling=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/peak_s3norm_antiTfe3_new_control
normalize_peak_calling=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/peakCalling
# s3norm_sample_file=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/modify_bedgraph/s3norm_antiTFE3_Samples/anti_TFE3_new_control.csv
# sample_list=($(awk -F "\"*\\t\"*" '{print $1}' $s3norm_sample_file)) 

sample_list=( "bulkChIC-PMC-DRO-014"\
            "bulkChIC-PMC-DRO-015"\
            "bulkChIC-PMC-DRO-016"\
            "SCC-bulkChIC-PMC-DRO-002"\
            "SCC-bulkChIC-PMC-DRO-005"\
            "SCC-bulkChIC-PMC-DRO-008"\
            "SCC-ChIC-PMC-DRO-L5"\
            "SCC-ChIC-PMC-DRO-F1"\
            "SCC-ChIC-PMC-DRO-F5"\
            "SCC-ChIC-PMC-DRO-T1"\
            "SCC-ChIC-PMC-DRO-T5"\
            "SCC-ChIC-PMC-DRO-L1" )
# sample_list=`awk -F "\"*\\t\"*" '{print $1}' /hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/modify_bedgraph/s3norm_antiTFE3_Samples/anti_TFE3_new_control.csv`
total_sample=${#sample_list[@]}
n=0
read_in_peak_list=()
total_reads_list=()
for sample_ID in ${sample_list[@]};do
  n=$((n+1))
  echo "-----------running $n out of $total_sample samples---------------"

#   reads_in_peaks=$(bedtools sort -i $normalize_peak_file/${sample_ID}*.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -nonamecheck -a $normalize_sample_count/${sample_ID}*.bedgraph -b stdin -ubam | samtools view -c)
# reads_in_peaks=$(awk '{ sum += $5; } END { print sum; }' "$@" $normalize_peak_calling/${sample_ID}*.narrowPeak )
#   read_in_peak_list=(${read_in_peak_list[@]} $reads_in_peaks)
reads_in_peaks=$(awk '{ sum += $5; } END { print sum; }' "$@" $normalize_peak_calling/${sample_ID}/narrow/${sample_ID}*.narrowPeak )

  read_in_peak_list=(${read_in_peak_list[@]} $reads_in_peaks)
  total_reads=$(awk '{ sum += $4; } END { print sum; }' "$@" $normalize_sample_count/${sample_ID}*.bedgraph)
  total_reads_list=(${total_reads_list[@]} $total_reads) 

#   total_reads_nbp= $(awk '{ sum += $4; } END { print sum; }' "$@" $normalize_sample_count/${sample_ID}*.bedgraph) 
done
save_name=raw_count_peak_file
echo ${total_reads_list[@]} > ${frip_dir}/${save_name}_total_reads.txt 
echo ${read_in_peak_list[@]} > ${frip_dir}/${save_name}_read_in_peak.txt
echo ${sample_list[@]} > ${frip_dir}/${save_name}_sample_IDs.txt


# calculate and plot frip in R
export FIGURE_DIR_VARIABLE=$figure_dir
export FRIP_DIR_VARIABLE=$frip_dir
export SAVE_NAME_VARIABLE=$save_name
Rscript plot_FRiP.R
echo "all done for FRiP calculation"