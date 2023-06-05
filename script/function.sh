#!/bin/bash
#SBATCH --job-name=module2
#SBATCH --output=after_peak_module.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=90G
#SBATCH --cpus-per-task=8
#SBATCH --gres=tmpspace:30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl


# module2: normalization, peak calling, frip plot, and diffBind

# prepare sample file for s3norm manually

# Define flexible variables (change everytime run for new set up)
save_name=antiTFE3_Samples
s3norm_new_output_dir=$s3norm_working_directory/with_new_control
mkdir -p $s3norm_new_output_dir
sample_list=(${anti_tfe3_sample_IDs[@]})
heatmap_samples_list=("fusion_merged" "tfe3_merged" "luciferase_merged" )
heatmap_sample_dir=${merged_bigwig_dir[@]}
s3norm_sample_file_name=anti_H3k4me3.csv

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
res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test
modify_bedgraph_dir=${res_dir}/modify_bedgraph
rm_dup_dir=$res_dir/rm_dup
diffBind_res_dir=${res_dir}/diffBind_analysis
merged_bigwig_dir=${res_dir}/merged_bigwig
figure_dir=${res_dir}/figures

s3norm_script_directory='/hpc/pmc_drost/nhung/S3norm'
s3norm_working_directory=${modify_bedgraph_dir}/s3norm_$save_name

bedgraph_dir=$s3norm_new_output_dir/S3norm_NBP_bedgraph
output_peak_dir=${res_dir}/peak_$save_name
mkdir -p ${output_peak_dir}

sample_count=$s3norm_new_output_dir/S3norm_rc_bedgraph
peak_calling=$output_peak_dir


refpath=$(find $diffBind_res_dir -name *$save_name*.bed)
heatmap_savename=DE_$save_name

# Normalize data
# . ./7-run_s3norm.sh 

# move s3norm output to a separate folder if needed
# mv $s3norm_working_directory/average_ref_bedgraph $s3norm_new_output_dir
# mv $s3norm_working_directory/NBP_bedgraph $s3norm_new_output_dir
# mv $s3norm_working_directory/S3norm_NBP_bedgraph $s3norm_new_output_dir
# mv $s3norm_working_directory/S3norm_rc_bedgraph $s3norm_new_output_dir

# call peak from normalize data
# echo "-------------------step 8. running peak calling after normalization----------------------"

# . ./8b-peakCalling_normalize.sh $bedgraph_dir $output_peak_dir $save_name ${sample_list[@]}  # worked! tested!

# plot frip 

# . ./7-Calculate_FRiP.sh $sample_count $output_peak_dir $save_name ${sample_list[@]}

## Identify DE sites
# generate sample sheet file (if there is no)
# echo "---------------step 11. generate sample sheet for peak differential analysis---------------"
# export SAMPLE_SHEET_DIR=$diffBind_res_dir
# export BAM_DIR=$rm_dup_dir
# export PEAK_DIR=$output_peak_dir 
# #clean_normalize_peak_dir
# export SAVE_NAME=$save_name
# Rscript Generate_diffbind_sample_sheet.R # need to change sample_IDs manually inside the script # worked! #tested!

# run diffBind
echo "---------------step 11. running peak differential analysis---------------"
export SAMPLE_SHEET_DIR_VARIABLE=$diffBind_res_dir/$save_name.csv
export SAVE_NAME_VARIABLE=$save_name
diffBind_res_sub_dir=$diffBind_res_dir/$SAVE_NAME_VARIABLE
mkdir -p $diffBind_res_sub_dir
export DIFFBIND_RESULT_DIR_VARIABLE=$diffBind_res_sub_dir
Rscript DiffBind_analysis.R

# Annotate DE site -> need to install packages
# path=($(find $diffBind_res_sub_dir -name *$save_name*.bed))
# export SAMPLE_BED_FILE=$path
# Rscript Annotate_peaks.R

# plot heatmap for DE site
echo "---------------step 12. running heatmap generation---------------"
. ./12-heatmap.sh "$heatmap_savename" "$refpath" "$heatmap_sample_dir" "${heatmap_samples_list[@]}"

echo "finished!"
echo "normalize data are in $s3norm_new_output_dir"
echo " heatmap figure is in $figure_dir"
echo " different peak sites and annotation are in $diffBind_res_sub_dir"
echo "frip figure is in $figure_dir"
echo "frip text is in $frip_dir "


