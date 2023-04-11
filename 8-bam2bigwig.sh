#!/bin/bash

# Merge all bam files from the same condition and then transform bam files to bigwig to prepare for heatmap generation

# Nhung 22 03 2023

# merge files from the same condition

# samtools merge -o ${merged_bigwig_dir}/tfe3_merged.bam ${rm_dup_dir}/SCC-ChIC-PMC-DRO-T1/SCC-ChIC-PMC-DRO-T1_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-ChIC-PMC-DRO-T5/SCC-ChIC-PMC-DRO-T5_rmdup_filt.bam \
# ${rm_dup_dir}/bulkChIC-PMC-DRO-016/bulkChIC-PMC-DRO-016_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-bulkChIC-PMC-DRO-008/SCC-bulkChIC-PMC-DRO-008_rmdup_filt.bam


# samtools merge -o ${merged_bigwig_dir}/luciferase_merged.bam ${rm_dup_dir}/SCC-ChIC-PMC-DRO-L1/SCC-ChIC-PMC-DRO-L1_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-ChIC-PMC-DRO-L5/SCC-ChIC-PMC-DRO-L5_rmdup_filt.bam \
# ${rm_dup_dir}/bulkChIC-PMC-DRO-014/bulkChIC-PMC-DRO-014_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-bulkChIC-PMC-DRO-005/SCC-bulkChIC-PMC-DRO-005_rmdup_filt.bam

# samtools merge -o ${merged_bigwig_dir}/fusion_merged.bam ${rm_dup_dir}/SCC-ChIC-PMC-DRO-F1/SCC-ChIC-PMC-DRO-F1_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-ChIC-PMC-DRO-F5/SCC-ChIC-PMC-DRO-F5_rmdup_filt.bam \
# ${rm_dup_dir}/bulkChIC-PMC-DRO-015/bulkChIC-PMC-DRO-015_rmdup_filt.bam \
# ${rm_dup_dir}/SCC-bulkChIC-PMC-DRO-002/SCC-bulkChIC-PMC-DRO-002_rmdup_filt.bam

##Index files
cd ${merged_bigwig_dir}
for INFILE in *.bam
do
   samtools index $INFILE $INFILE".bai"
done

##transform to bigwig
$bamCoverage_dir -b tfe3_merged.bam \
    --normalizeUsing CPM \
    --effectiveGenomeSize 2913022398 \
    --ignoreForNormalization chrX \
    -o ${merged_bigwig_dir}/tfe3_merged.bw

$bamCoverage_dir -b luciferase_merged.bam \
    --normalizeUsing CPM \
    --effectiveGenomeSize 2913022398 \
    --ignoreForNormalization chrX \
    -o ${merged_bigwig_dir}/luciferase_merged.bw

$bamCoverage_dir -b fusion_merged.bam \
    --normalizeUsing CPM \
    --effectiveGenomeSize 2913022398 \
    --ignoreForNormalization chrX \
    -o ${merged_bigwig_dir}/fusion_merged.bw
