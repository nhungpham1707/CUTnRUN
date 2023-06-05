#!/bin/bash
#SBATCH --job-name=module3
#SBATCH --output=generate_igv_track.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=90G
#SBATCH --cpus-per-task=8
#SBATCH --gres=tmpspace:30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl

# Generate signal track for IGV visualization
# Nhung 01 06 2023 
# . ./Define_variable.sh 
# merge NBP bedgraph files from the same group
# signal_track_path=$s3norm_new_output_dir/NBP_bedgraph
# IGV_bedgraph_out_dir=$IGV_input_dir

# tfe3_bg_savename='new_control_tfe3'
# fusion_bg_savename='new_control_fusion'
# luc_bg_savename='new_control_luciferase'

# echo "start bedgraph union from samples in the same group"
# . ./Merge_bedgraphs.sh $tfe3_bg_savename $signal_track_path $IGV_bedgraph_out_dir ${tfe3[@]}

# . ./Merge_bedgraphs.sh $fusion_bg_savename $signal_track_path $IGV_bedgraph_out_dir ${fusion[@]}

# . ./Merge_bedgraphs.sh $luc_bg_savename $signal_track_path $IGV_bedgraph_out_dir ${luciferase[@]}

# # calculate mean from replicates from merg bedgraph files
# echo "start merging replicates in the bedgraph union files"
# export RESULT_DIR_VARIABLE=$IGV_bedgraph_out_dir

# export FILE_PATH_VARIABLE=$IGV_bedgraph_out_dir/${tfe3_bg_savename}_union.bedgraph
# export SAVE_NAME_VARIABLE=$tfe3_bg_savename
# Rscript Calculate_mean_bedgraph.R 

# export FILE_PATH_VARIABLE=$IGV_bedgraph_out_dir/${fusion_bg_savename}_union.bedgraph
# export SAVE_NAME_VARIABLE=$fusion_bg_savename
# Rscript Calculate_mean_bedgraph.R 

# export FILE_PATH_VARIABLE=$IGV_bedgraph_out_dir/${luc_bg_savename}_union.bedgraph
# export SAVE_NAME_VARIABLE=$luc_bg_savename
# Rscript Calculate_mean_bedgraph.R 

# convert bedgraph to bigwig
# echo "convert bedgraph to bigwig"

# res_dir=$IGV_bedgraph_out_dir
# hg38_size_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/bigwig
# save_name=$fusion_bg_savename
# input_name=(find $IGV_bedgraph_out_dir -name *${save_name}_mean_rep.bedgraph)
# . ./Convert_bedgraph_to_bw.sh $input_name $res_dir $save_name $hg38_size_dir

# save_name=$luc_bg_savename
# input_name=(find $IGV_bedgraph_out_dir -name *${save_name}_mean_rep.bedgraph)
# . ./Convert_bedgraph_to_bw.sh $input_name $res_dir $save_name $hg38_size_dir

# save_name=$tfe3_bg_savename
# input_name=(find $IGV_bedgraph_out_dir -name *${save_name}_mean_rep.bedgraph)
# . ./Convert_bedgraph_to_bw.sh $input_name $res_dir $save_name $hg38_size_dir

#### histone antibodies

res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test
IGV_input_dir=${res_dir}/IGV_input
IGV_bedgraph_out_dir=$IGV_input_dir

h3k4me3=( "bulkChIC-PMC-DRO-011" "bulkChIC-PMC-DRO-012" "bulkChIC-PMC-DRO-013")
h3k4me1=("SCC-bulkChIC-PMC-DRO-021" "SCC-bulkChIC-PMC-DRO-024")
h3k27ac=("SCC-bulkChIC-PMC-DRO-022" "SCC-bulkChIC-PMC-DRO-025")
sfpq=("SCC-bulkChIC-PMC-DRO-020" "SCC-bulkChIC-PMC-DRO-023")

me1_signal_track_path=${res_dir}/modify_bedgraph/s3norm_anti_H3k4me1_samples/NBP_bedgraph
h3k4me1_bg_savename='antiH3k4me1'

me3_signal_track_path=${res_dir}/modify_bedgraph/s3norm_anti_H3k4me3_samples/NBP_bedgraph
h3k4me3_bg_savename='antiH3k4me4'

h3k27ac_signal_track_path=${res_dir}/modify_bedgraph/s3norm_anti_H3k27ac_samples/NBP_bedgraph
h3k27ac_bg_savename='antiH3k27ac'

sfpq_signal_track_path=${res_dir}/modify_bedgraph/s3norm_anti_SFPQ_samples/NBP_bedgraph
sfpq_bg_savename='antiSFPQ'

echo "start bedgraph union from samples in the same group"

. ./Merge_bedgraphs.sh $h3k4me1_bg_savename $me1_signal_track_path $IGV_bedgraph_out_dir ${h3k4me1[@]}

. ./Merge_bedgraphs.sh $h3k4me3_bg_savename $me3_signal_track_path $IGV_bedgraph_out_dir ${h3k4me3[@]}

. ./Merge_bedgraphs.sh $h3k27ac_bg_savename $h3k27ac_signal_track_path $IGV_bedgraph_out_dir ${h3k27ac[@]}

. ./Merge_bedgraphs.sh $sfpq_bg_savename $sfpq_signal_track_path $IGV_bedgraph_out_dir ${sfpq[@]}
# calculate mean from replicates from merg bedgraph files
echo "start merging replicates in the bedgraph union files"
export RESULT_DIR_VARIABLE=$IGV_bedgraph_out_dir

export FILE_PATH_VARIABLE=$IGV_bedgraph_out_dir/${sfpq_bg_savename}_union.bedgraph
export SAVE_NAME_VARIABLE=$sfpq_bg_savename
Rscript Calculate_mean_bedgraph.R 

export FILE_PATH_VARIABLE=$IGV_bedgraph_out_dir/${h3k4me3_bg_savename}_union.bedgraph
export SAVE_NAME_VARIABLE=$h3k4me3_bg_savename
Rscript Calculate_mean_bedgraph.R 

export FILE_PATH_VARIABLE=$IGV_bedgraph_out_dir/${h3k4me1_bg_savename}_union.bedgraph
export SAVE_NAME_VARIABLE=$h3k4me1_bg_savename
Rscript Calculate_mean_bedgraph.R 

export FILE_PATH_VARIABLE=$IGV_bedgraph_out_dir/${h3k27ac_bg_savename}_union.bedgraph
export SAVE_NAME_VARIABLE=$h3k27ac_bg_savename
Rscript Calculate_mean_bedgraph.R 


echo "finished!"
echo "files that can be visualized in IGV are sorted bedgraph files (mean, merge loci and sort based on chr) in $IGV_bedgraph_out_dir" 
echo "files for plotting heatmap are bigwig files in $IGV_bedgraph_out_dir"