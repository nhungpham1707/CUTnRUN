#!/bin/bash
# this script extract peaks from each condition after peakcalling. Either by intersect, merge or multiinter
# Nhung 11 04 2023

########## find overlap peaks among replicates in the same condition #################

# bedtools multiinter -i /hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/peakCalling/SCC-ChIC-PMC-DRO-F1/narrow/SCC-ChIC-PMC-DRO-F1_paired_control_peaks.narrowPeak /hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/peakCalling/SCC-ChIC-PMC-DRO-F5/narrow/SCC-ChIC-PMC-DRO-F5_paired_control_peaks.narrowPeak | awk '{if ($4 == 2) print}' > test.txt


# # tfe3
# bedtools multiinter -i ${peak_dir}/SCC-ChIC-PMC-DRO-T1/narrow/SCC-ChIC-PMC-DRO-T1_paired_control_peaks.narrowPeak ${peak_dir}/bulkChIC-PMC-DRO-016/narrow/bulkChIC-PMC-DRO-016_paired_control_peaks.narrowPeak ${peak_dir}/SCC-bulkChIC-PMC-DRO-008/narrow/SCC-bulkChIC-PMC-DRO-008_paired_control_peaks.narrowPeak ${peak_dir}/SCC-ChIC-PMC-DRO-T5/narrow/SCC-ChIC-PMC-DRO-T5_paired_control_peaks.narrowPeak -header > ${peak_analysis_dir}/tfe3_multiinter.txt

# # luc 
# bedtools multiinter -i ${peak_dir}/SCC-ChIC-PMC-DRO-L1/narrow/SCC-ChIC-PMC-DRO-L1_paired_control_peaks.narrowPeak ${peak_dir}/SCC-ChIC-PMC-DRO-L5/narrow/SCC-ChIC-PMC-DRO-L5_paired_control_peaks.narrowPeak ${peak_dir}/bulkChIC-PMC-DRO-014/narrow/bulkChIC-PMC-DRO-014_paired_control_peaks.narrowPeak ${peak_dir}/SCC-bulkChIC-PMC-DRO-005/narrow/SCC-bulkChIC-PMC-DRO-005_paired_control_peaks.narrowPeak -header > ${peak_analysis_dir}/luc_multiinter.txt

# fusion 
# bedtools multiinter -i ${peak_dir}/SCC-ChIC-PMC-DRO-F1/narrow/SCC-ChIC-PMC-DRO-F1_paired_control_peaks.narrowPeak ${peak_dir}/SCC-ChIC-PMC-DRO-F5/narrow/SCC-ChIC-PMC-DRO-F5_paired_control_peaks.narrowPeak ${peak_dir}/bulkChIC-PMC-DRO-015/narrow/bulkChIC-PMC-DRO-015_paired_control_peaks.narrowPeak ${peak_dir}/SCC-bulkChIC-PMC-DRO-002/narrow/SCC-bulkChIC-PMC-DRO-002_paired_control_peaks.narrowPeak -header > ${peak_analysis_dir}/fusion_multiinter.txt

## histone samples
# bedtools multiinter -i ${peak_no_control_dir}/SCC-ChIC-PMC-DRO-FH/narrow/SCC-ChIC-PMC-DRO-FH_paired_peaks.narrowPeak ${peak_no_control_dir}/SCC-ChIC-PMC-DRO-LH/narrow/SCC-ChIC-PMC-DRO-LH_paired_peaks.narrowPeak ${peak_no_control_dir}/SCC-ChIC-PMC-DRO-TH/narrow/SCC-ChIC-PMC-DRO-TH_paired_peaks.narrowPeak ${peak_no_control_dir}/bulkChIC-PMC-DRO-011/narrow/bulkChIC-PMC-DRO-011_paired_peaks.narrowPeak ${peak_no_control_dir}/bulkChIC-PMC-DRO-012/narrow/bulkChIC-PMC-DRO-012_paired_peaks.narrowPeak ${peak_no_control_dir}/bulkChIC-PMC-DRO-013/narrow/bulkChIC-PMC-DRO-013_paired_peaks.narrowPeak -header > ${peak_analysis_dir}/h3k4me3_multiinter.txt

###### merge peak from the same condition
# samtools merge -o ${peak_analysis_dir}/tfe3_merged.bam ${rm_dup_dir}/SCC-ChIC-PMC-DRO-T1/SCC-ChIC-PMC-DRO-T1_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-ChIC-PMC-DRO-T5/SCC-ChIC-PMC-DRO-T5_rmdup_filt.bam \
# ${rm_dup_dir}/bulkChIC-PMC-DRO-016/bulkChIC-PMC-DRO-016_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-bulkChIC-PMC-DRO-008/SCC-bulkChIC-PMC-DRO-008_rmdup_filt.bam

# samtools merge -o ${peak_analysis_dir}/luciferase_merged.bam ${rm_dup_dir}/SCC-ChIC-PMC-DRO-L1/SCC-ChIC-PMC-DRO-L1_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-ChIC-PMC-DRO-L5/SCC-ChIC-PMC-DRO-L5_rmdup_filt.bam \
# ${rm_dup_dir}/bulkChIC-PMC-DRO-014/bulkChIC-PMC-DRO-014_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-bulkChIC-PMC-DRO-005/SCC-bulkChIC-PMC-DRO-005_rmdup_filt.bam

# samtools merge -o ${peak_analysis_dir}/fusion_merged.bam ${rm_dup_dir}/SCC-ChIC-PMC-DRO-F1/SCC-ChIC-PMC-DRO-F1_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-ChIC-PMC-DRO-F5/SCC-ChIC-PMC-DRO-F5_rmdup_filt.bam \
# ${rm_dup_dir}/bulkChIC-PMC-DRO-015/bulkChIC-PMC-DRO-015_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-bulkChIC-PMC-DRO-002/SCC-bulkChIC-PMC-DRO-002_rmdup_filt.bam

# samtools merge -o ${peak_analysis_dir}/h3k4me3_merged.bam ${rm_dup_dir}/SCC-ChIC-PMC-DRO-FH/SCC-ChIC-PMC-DRO-FH_rmdup_filt.bam ${rm_dup_dir}/SCC-ChIC-PMC-DRO-LH/SCC-ChIC-PMC-DRO-LH_rmdup_filt.bam ${rm_dup_dir}/SCC-ChIC-PMC-DRO-TH/SCC-ChIC-PMC-DRO-TH_rmdup_filt.bam ${rm_dup_dir}/bulkChIC-PMC-DRO-011/bulkChIC-PMC-DRO-011_rmdup_filt.bam ${rm_dup_dir}/bulkChIC-PMC-DRO-012/bulkChIC-PMC-DRO-012_rmdup_filt.bam ${rm_dup_dir}/bulkChIC-PMC-DRO-013/bulkChIC-PMC-DRO-013_rmdup_filt.bam

##Index files
# cd ${peak_analysis_dir}
# for INFILE in *.bam
# do
#    samtools index $INFILE $INFILE".bai"
# done

## find overlap. More details on flags and input can be found here https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html 
# # tfe3
# bedtools intersect -a ${peak_dir}/SCC-ChIC-PMC-DRO-T1/narrow/SCC-ChIC-PMC-DRO-T1_paired_control_peaks.narrowPeak -b ${peak_dir}/bulkChIC-PMC-DRO-016/narrow/bulkChIC-PMC-DRO-016_paired_control_peaks.narrowPeak ${peak_dir}/SCC-bulkChIC-PMC-DRO-008/narrow/SCC-bulkChIC-PMC-DRO-008_paired_control_peaks.narrowPeak ${peak_dir}/SCC-ChIC-PMC-DRO-T5/narrow/SCC-ChIC-PMC-DRO-T5_paired_control_peaks.narrowPeak -names bulkChIC-PMC-DRO-016 SCC-bulkChIC-PMC-DRO-008 SCC-ChIC-PMC-DRO-T5 -sorted  > ${peak_analysis_dir}/tfe3_intersect.txt

# # luc 
# bedtools intersect -a ${peak_dir}/SCC-ChIC-PMC-DRO-L1/narrow/SCC-ChIC-PMC-DRO-L1_paired_control_peaks.narrowPeak -b ${peak_dir}/SCC-ChIC-PMC-DRO-L5/narrow/SCC-ChIC-PMC-DRO-L5_paired_control_peaks.narrowPeak ${peak_dir}/bulkChIC-PMC-DRO-014/narrow/bulkChIC-PMC-DRO-014_paired_control_peaks.narrowPeak ${peak_dir}/SCC-bulkChIC-PMC-DRO-005/narrow/SCC-bulkChIC-PMC-DRO-005_paired_control_peaks.narrowPeak -names SCC-ChIC-PMC-DRO-L5 bulkChIC-PMC-DRO-014 SCC-bulkChIC-PMC-DRO-005 -sorted > ${peak_analysis_dir}/luc_intersect.txt

# fusion 

# bedtools intersect -i ${peak_dir}/SCC-ChIC-PMC-DRO-F1/narrow/SCC-ChIC-PMC-DRO-F1_paired_control_peaks.narrowPeak ${peak_dir}/SCC-ChIC-PMC-DRO-F5/narrow/SCC-ChIC-PMC-DRO-F5_paired_control_peaks.narrowPeak ${peak_dir}/bulkChIC-PMC-DRO-015/narrow/bulkChIC-PMC-DRO-015_paired_control_peaks.narrowPeak ${peak_dir}/SCC-bulkChIC-PMC-DRO-002/narrow/SCC-bulkChIC-PMC-DRO-002_paired_control_peaks.narrowPeak -names SCC-ChIC-PMC-DRO-F5 bulkChIC-PMC-DRO-015 SCC-bulkChIC-PMC-DRO-002 -sorted > ${peak_analysis_dir}/fusion_intersect.txt
# get fusion sample

declare -a fusion_IDs=()
total_sample=${#fusion[@]}
n=0
for sample_ID in ${fusion[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample fusion samples---------------------- "

    fusion_IDs+=( "${peak_dir}/${sample_ID}/narrow/${sample_ID}_paired_control_peaks.narrowPeak" )

done

bedtools intersect -a ${fusion_IDs[0]} -b ${fusion_IDs[@]:1} -names ${fusion[@]} -sorted > ${peak_analysis_dir}/luc_intersect.txt
## histone samples
# bedtools intersect -a ${peak_no_control_dir}/SCC-ChIC-PMC-DRO-FH/narrow/SCC-ChIC-PMC-DRO-FH_paired_peaks.narrowPeak -b ${peak_no_control_dir}/SCC-ChIC-PMC-DRO-LH/narrow/SCC-ChIC-PMC-DRO-LH_paired_peaks.narrowPeak ${peak_no_control_dir}/SCC-ChIC-PMC-DRO-TH/narrow/SCC-ChIC-PMC-DRO-TH_paired_peaks.narrowPeak ${peak_no_control_dir}/bulkChIC-PMC-DRO-011/narrow/bulkChIC-PMC-DRO-011_paired_peaks.narrowPeak ${peak_no_control_dir}/bulkChIC-PMC-DRO-012/narrow/bulkChIC-PMC-DRO-012_paired_peaks.narrowPeak ${peak_no_control_dir}/bulkChIC-PMC-DRO-013/narrow/bulkChIC-PMC-DRO-013_paired_peaks.narrowPeak -names SCC-ChIC-PMC-DRO-LH SCC-ChIC-PMC-DRO-TH bulkChIC-PMC-DRO-012 bulkChIC-PMC-DRO-013 bulkChIC-PMC-DRO-011 -sorted > ${peak_analysis_dir}/h3k4me3_multiinter.txt

# get histone file names
declare -a histone_IDs=()
total_sample=${#h3k4me3_all[@]}
n=0
for sample_ID in ${h3k4me3_all[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample samples---------------------- "

    histone_IDs+=( "${peak_no_control_dir}/${sample_ID}/narrow/${sample_ID}_paired_peaks.narrowPeak" )

done

 bedtools intersect -a ${histone_IDs[0]} -b ${histone_IDs[@]:1} -names ${h3k4me3_all[@]} -sorted > ${peak_analysis_dir}/h3k4me3_intersect.txt

bedtools multiinter -i ${histone_IDs[@]} -header > ${peak_analysis_dir}/h3k4me3_multiinter.txt