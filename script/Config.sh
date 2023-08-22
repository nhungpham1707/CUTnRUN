#!/bin/bash

## Define variables for all modules

# directory to data and result

data_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/SCC_ChIC-PMC-DRO_plates_20210520_run1
# result dir
res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test

qc_dir=${res_dir}/qualityCheck

trim_dir=${res_dir}/trimmed

align_dir=${res_dir}/alignment

rm_dup_dir=$res_dir/rm_dup

frip_dir=${res_dir}/FRiP

modify_bedgraph_dir=${res_dir}/modify_bedgraph

peak_dir=${res_dir}/peakCalling

peak_no_control_dir=${res_dir}/peakCalling_nocontrol

normalize_peak_dir=${res_dir}/normalize_peakCalling
mkdir -p ${normalize_peak_dir}

motif_dir=${res_dir}/motif
mkdir -p ${motif_dir}

merged_bigwig_dir=${res_dir}/merged_bigwig

bigwig_dir=${res_dir}/bigwig

peak_analysis_dir=${res_dir}/peak_analysis
mkdir -p $peak_analysis_dir

figure_dir=${res_dir}/figures

diffBind_res_dir=${res_dir}/diffBind_analysis

# Define tool dir 
bowtie2Index=/hpc/pmc_drost/SOURCES/Genomes/human/bowtie2/human_gencode37_hg38

picardTool=/hpc/pmc_drost/PROJECTS/swang/software/picard.jar
new_tmp_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/tmp # to solve the out of space with the temporary output from picard

homer_dir=/hpc/pmc_drost/nhung/anaconda3/envs/cutnrun_trimgalore/bin/homer
findMotif_dir=/hpc/pmc_drost/nhung/anaconda3/envs/cutnrun_trimgalore/bin/findMotifsGenome.pl

bamCoverage_dir=/hpc/pmc_drost/nhung/anaconda3/envs/cutnrun_trimgalore/bin/bamCoverage # to merge 
effectiveGenomeSize=2913022398 
# fasta genome directory for motif finding. It has to be the same with the reference genome that was used for alignment
fasta_genome_dir=/hpc/pmc_drost/SOURCES/Genomes/human/gencode37_GRCh38_primary_assembly_genome.fa
 
# Define samples
# sample_IDs=( "bulkChIC-PMC-DRO-011" \
#             "bulkChIC-PMC-DRO-012"\
#             "bulkChIC-PMC-DRO-013"\
#             "bulkChIC-PMC-DRO-014"\
#             "bulkChIC-PMC-DRO-015"\
#             "bulkChIC-PMC-DRO-016"\
#             "SCC-bulkChIC-PMC-DRO-002"\
#             "SCC-bulkChIC-PMC-DRO-005"\
#             "SCC-bulkChIC-PMC-DRO-008"\
#             "SCC-ChIC-PMC-DRO-L5"\
#             "SCC-ChIC-PMC-DRO-LH"\
#             "SCC-ChIC-PMC-DRO-F1"\
#             "SCC-ChIC-PMC-DRO-F5"\
#             "SCC-ChIC-PMC-DRO-FH"\
#             "SCC-ChIC-PMC-DRO-T1"\
#             "SCC-ChIC-PMC-DRO-T5"\
#             "SCC-ChIC-PMC-DRO-TH"\ 
#             "SCC-ChIC-PMC-DRO-L1"\ 
#             "SCC-bulkChIC-PMC-DRO-020"\
#             "SCC-bulkChIC-PMC-DRO-021"\
#             "SCC-bulkChIC-PMC-DRO-022"\
#             "SCC-bulkChIC-PMC-DRO-023"\
#             "SCC-bulkChIC-PMC-DRO-024"\
#             "SCC-bulkChIC-PMC-DRO-025")

# sample_IDs=( "bulkChIC-PMC-DRO-011" \
#             "bulkChIC-PMC-DRO-012")
sample_IDs=('SCC-bulkChIC-PMC-DRO-020' 'SCC-bulkChIC-PMC-DRO-023')

total_sample=${#sample_IDs[@]}

# Define sample groups

tfe3=( "SCC-ChIC-PMC-DRO-T1" "SCC-ChIC-PMC-DRO-T5" "bulkChIC-PMC-DRO-016" "SCC-bulkChIC-PMC-DRO-008")
luciferase=( "SCC-ChIC-PMC-DRO-L1" "SCC-ChIC-PMC-DRO-L5" "bulkChIC-PMC-DRO-014" "SCC-bulkChIC-PMC-DRO-005")
fusion=( "SCC-ChIC-PMC-DRO-F1" "SCC-ChIC-PMC-DRO-F5" "bulkChIC-PMC-DRO-015" "SCC-bulkChIC-PMC-DRO-002")

tfe3C=$res_dir/rm_dup/bulkChIC-PMC-DRO-013/bulkChIC-PMC-DRO-013_rmdup_filt.bam
luciferaseC=$res_dir/rm_dup/bulkChIC-PMC-DRO-011/bulkChIC-PMC-DRO-011_rmdup_filt.bam
fusionC=$res_dir/rm_dup/bulkChIC-PMC-DRO-012/bulkChIC-PMC-DRO-012_rmdup_filt.bam
# general set up 
 turn of core file generation to not overload home dir
 ulimit -c 0 # to disable coredump  https://community.hpe.com/t5/system-administration/how-to-disable-or-restrict-core-dumps/td-p/4276097#.ZDzoIpNByDU

 # set up for module 1
frip_sample_count=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/modify_bedgraph/s3norm_antiTFE3_Samples/with_new_control/S3norm_rc_bedgraph 
frip_peak_count=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/peak_s3norm_antiTfe3_new_control
frip_save_name=normalize_anti_TFE3_new_control
frip_samples=${anti_tfe3_sample_IDs[@]}


###checks if file extension is .fq or .fastq
input_format=$(ls -R $data_dir/${sample_IDs[1]}/*.gz | head -1)
if [[ $input_format == *"fastq"* ]]
then
input_extension="fastq.gz"
elif [[ $input_format == *"fq"* ]]
then
input_extension="fq.gz"
fi

echo $input_extension

# Module 2: variables for data normalization
bin_size=200

export s3norm_sample_list=${sample_IDs[0]}
for i in "${sample_IDs[@]:1}"; do
   s3norm_sample_list+=,$i
done

## where s3norm was downloaded
s3norm_script_directory='/hpc/pmc_drost/nhung/S3norm'
# s3norm_yichao allows without control samples
# s3norm_script_directory='/hpc/pmc_drost/nhung/s3norm_yichao/S3norm'

s3norm_working_directory=$modify_bedgraph_dir 
# sample list csv file and bedgraph files need to be in this working dir. result will be stored in the working dir. 
s3norm_sample_file_name=anti_SFPQ_w_control 

# directory to save peak from normalize data
peak_after_s3norm_dir=$normalize_peak_dir/$s3norm_sample_file_name

# module 3: peak analysis
## variables for generating sample sheet for diffBind
diffbind_save_name='diffbind_normalize_samples' # the same name for the sample sheet file
export SAMPLE_SHEET_DIR=$diffBind_res_dir
export BAM_DIR=$rm_dup_dir
export PEAK_DIR=$peak_after_s3norm_dir 
#clean_normalize_peak_dir
export SAVE_NAME=$diffbind_save_name