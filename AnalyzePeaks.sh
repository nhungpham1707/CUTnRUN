#!/bin/bash
# this script extract peaks from each condition after peakcalling. 
# Nhung 11 04 2023

########## find overlap peaks among replicates in the same condition #################

# bedtools multiinter -i /hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/peakCalling/SCC-ChIC-PMC-DRO-F1/narrow/SCC-ChIC-PMC-DRO-F1_paired_control_peaks.narrowPeak /hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/peakCalling/SCC-ChIC-PMC-DRO-F5/narrow/SCC-ChIC-PMC-DRO-F5_paired_control_peaks.narrowPeak | awk '{if ($4 == 2) print}' > test.txt


# tfe3
bedtools multiinter -i ${peak_dir}/SCC-ChIC-PMC-DRO-T1/narrow/SCC-ChIC-PMC-DRO-T1_paired_control_peaks.narrowPeak ${peak_dir}/bulkChIC-PMC-DRO-016/narrow/bulkChIC-PMC-DRO-016_paired_control_peaks.narrowPeak ${peak_dir}/SCC-bulkChIC-PMC-DRO-008/narrow/SCC-bulkChIC-PMC-DRO-008_paired_control_peaks.narrowPeak ${peak_dir}/SCC-ChIC-PMC-DRO-T5/narrow/SCC-ChIC-PMC-DRO-T5_paired_control_peaks.narrowPeak -header > ${peak_analysis_dir}/tfe3_multiinter.txt

# luc 
bedtools multiinter -i ${peak_dir}/SCC-ChIC-PMC-DRO-L1/narrow/SCC-ChIC-PMC-DRO-L1_paired_control_peaks.narrowPeak ${peak_dir}/SCC-ChIC-PMC-DRO-L5/narrow/SCC-ChIC-PMC-DRO-L5_paired_control_peaks.narrowPeak ${peak_dir}/bulkChIC-PMC-DRO-014/narrow/bulkChIC-PMC-DRO-014_paired_control_peaks.narrowPeak ${peak_dir}/SCC-bulkChIC-PMC-DRO-005/narrow/SCC-bulkChIC-PMC-DRO-005_paired_control_peaks.narrowPeak -header > ${peak_analysis_dir}/luc_multiinter.txt

# fusion 
bedtools multiinter -i ${peak_dir}/SCC-ChIC-PMC-DRO-F1/narrow/SCC-ChIC-PMC-DRO-F1_paired_control_peaks.narrowPeak ${peak_dir}/SCC-ChIC-PMC-DRO-F5/narrow/SCC-ChIC-PMC-DRO-F5_paired_control_peaks.narrowPeak ${peak_dir}/bulkChIC-PMC-DRO-015/narrow/bulkChIC-PMC-DRO-015_paired_control_peaks.narrowPeak ${peak_dir}/SCC-bulkChIC-PMC-DRO-002/narrow/SCC-bulkChIC-PMC-DRO-002_paired_control_peaks.narrowPeak -header > ${peak_analysis_dir}/fusion_multiinter.txt

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

##Index files
# cd ${peak_analysis_dir}
# for INFILE in *.bam
# do
#    samtools index $INFILE $INFILE".bai"
# done