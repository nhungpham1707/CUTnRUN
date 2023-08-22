# module3: peak analysis: differential analysis, occupancy analysis and annotation
mkdir -p $diffBind_res_dir

# Identify DE sites
# generate sample sheet file (if there is no)
echo "---------------step 11. generate sample sheet for peak differential analysis---------------"

Rscript Generate_diffbind_sample_sheet.R # need to change sample_IDs manually inside the script # worked! #tested!

# run diffBind
echo "---------------step 11. running peak differential analysis---------------"
export SAMPLE_SHEET_DIR_VARIABLE=$diffBind_res_dir/${diffbind_save_name}.csv
export SAVE_NAME_VARIABLE=$diffbind_save_name
diffBind_res_sub_dir=$diffBind_res_dir/$SAVE_NAME_VARIABLE
mkdir -p $diffBind_res_sub_dir
export DIFFBIND_RESULT_DIR_VARIABLE=$diffBind_res_sub_dir
Rscript DiffBind_analysis.R

# occupancyn analysis and peak annotation are conducted on a different R script on a local machine