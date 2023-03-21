#!/bin/bash
# This script defines variables that are needed to run all scripts to analyze cut and run data. Modify these variables when running for new samples or using differnt tools or directory 
# Nhung, 21 03 2023

data_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/SCC_ChIC-PMC-DRO_plates_20210520_run1

# result dir
res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test

qc_dir=${res_dir}/fastqc
mkdir -p $qc_dir

trim_dir=${res_dir}/trimmed
mkdir -p $trim_dir

align_dir=${res_dir}/alignment
mkdir -p $align_dir

rm_dup_dir=$res_dir/rm_dup
mkdir -p $rm_dup_dir

peak_dir=${res_dir}/peakCalling/${sample_ID}
mkdir -p ${peak_dir}

peak_no_control_dir=${res_dir}/peakCalling_nocontrol
mkdir -p ${peak_no_control_dir}

# tool dir
bowtie2Index=/hpc/pmc_drost/SOURCES/Genomes/human/bowtie2/human_gencode37_hg38

picardTool=/hpc/pmc_drost/PROJECTS/swang/software/picard.jar
new_tmp_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/tmp # to solve the out of space with the temporary output from picard

# sample_Ids was generated as text file from ls filename.txt from the data_dir

sample_IDs=( "bulkChIC-PMC-DRO-011" \
            "bulkChIC-PMC-DRO-012"\
            "bulkChIC-PMC-DRO-013"\
            "bulkChIC-PMC-DRO-014"\
            "bulkChIC-PMC-DRO-015"\
            "bulkChIC-PMC-DRO-016"\
            "SCC-bulkChIC-PMC-DRO-002"\
            "SCC-bulkChIC-PMC-DRO-005"\
            "SCC-bulkChIC-PMC-DRO-008"\
            "SCC-ChIC-PMC-DRO-L5"\
            "SCC-ChIC-PMC-DRO-LH"\
            "SCC-ChIC-PMC-DRO-F1"\
            "SCC-ChIC-PMC-DRO-F5"\
            "SCC-ChIC-PMC-DRO-FH"\
            "SCC-ChIC-PMC-DRO-T1"\
            "SCC-ChIC-PMC-DRO-T5"\
            "SCC-ChIC-PMC-DRO-TH"\ 
            "SCC-ChIC-PMC-DRO-L1")

total_sample=${#sample_IDs[@]}

# classify sample for peakcalling
tfe3="SCC-ChIC-PMC-DRO-T1 SCC-ChIC-PMC-DRO-T5 bulkChIC-PMC-DRO-016 SCC-bulkChIC-PMC-DRO-008"
luciferase="SCC-ChIC-PMC-DRO-L1 SCC-ChIC-PMC-DRO-L5 bulkChIC-PMC-DRO-014 SCC-bulkChIC-PMC-DRO-005"
fusion="SCC-ChIC-PMC-DRO-F1 SCC-ChIC-PMC-DRO-F5 bulkChIC-PMC-DRO-015 SCC-bulkChIC-PMC-DRO-002"
allT="$tfe3 $luciferase $fusion"
tfe3C=$res_dir/rm_dup/bulkChIC-PMC-DRO-013/bulkChIC-PMC-DRO-013_rmdup_filt.bam
luciferaseC=$res_dir/rm_dup/bulkChIC-PMC-DRO-011/bulkChIC-PMC-DRO-011_rmdup_filt.bam
fusionC=$res_dir/rm_dup/bulkChIC-PMC-DRO-012/bulkChIC-PMC-DRO-012_rmdup_filt.bam

# source variables to other scripts in the pipeline 
. ./qualityCheck.sh
. ./2-trimming.sh
. ./3-alignment.sh
. ./4-filtering.sh
. ./5-peakCalling.sh






