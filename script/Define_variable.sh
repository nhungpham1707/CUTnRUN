#!/bin/bash

# Define fix variable 
anti_tfe3_sample_IDs=( "bulkChIC-PMC-DRO-014"\
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
            tfe3=( "SCC-ChIC-PMC-DRO-T1" "SCC-ChIC-PMC-DRO-T5" "bulkChIC-PMC-DRO-016" "SCC-bulkChIC-PMC-DRO-008")
luciferase=( "SCC-ChIC-PMC-DRO-L1" "SCC-ChIC-PMC-DRO-L5" "bulkChIC-PMC-DRO-014" "SCC-bulkChIC-PMC-DRO-005")
fusion=( "SCC-ChIC-PMC-DRO-F1" "SCC-ChIC-PMC-DRO-F5" "bulkChIC-PMC-DRO-015" "SCC-bulkChIC-PMC-DRO-002")

res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test
modify_bedgraph_dir=${res_dir}/modify_bedgraph
rm_dup_dir=$res_dir/rm_dup
diffBind_res_dir=${res_dir}/diffBind_analysis
merged_bigwig_dir=${res_dir}/merged_bigwig
figure_dir=${res_dir}/figures
IGV_input_dir=${res_dir}/IGV_input
mkdir -p $IGV_input_dir

save_name=antiTFE3_Samples
s3norm_script_directory='/hpc/pmc_drost/nhung/S3norm'
s3norm_working_directory=${modify_bedgraph_dir}/s3norm_$save_name

s3norm_new_output_dir=$s3norm_working_directory/with_new_control
# mkdir -p $s3norm_new_output_dir
sample_list=(${anti_tfe3_sample_IDs[@]})
heatmap_samples_list=("fusion_merged" "tfe3_merged" "luciferase_merged" )
heatmap_sample_dir=${merged_bigwig_dir[@]}
s3norm_sample_file_name=anti_H3k4me3.csv

bedgraph_dir=$s3norm_new_output_dir/S3norm_NBP_bedgraph
output_peak_dir=${res_dir}/peak_$save_name
mkdir -p ${output_peak_dir}

sample_count=$s3norm_new_output_dir/S3norm_rc_bedgraph
peak_calling=$output_peak_dir

refpath=$diffBind_res_dir/*$save_name*.bed
heatmap_savename=DE_$save_name