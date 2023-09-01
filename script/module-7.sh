# modify these variables below before running

# define which sample to plot and where are the bw files 
samples_list=("fusion_merged" "luciferase_merged" )
sample_dir=${merged_bigwig_dir[@]}

# define genome region to plot
refpath=$diffBind_res_dir/Rerun_14062023/2023-06-14common_peaks_wo_sig_peaks_fusion_luc_s3norm_all_except_low_dephistone.bed 

# define plot save name
savename=common_peak_wo_tfe3_increase_size.svg
. ./12-heatmap.sh "$savename" "$refpath" "$sample_dir" "${samples_list[@]}"

