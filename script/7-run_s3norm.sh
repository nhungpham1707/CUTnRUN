#!/bin/bash

# normalize data using s3norm package
# nhung 25 04 2023

## set up 
# conda env s3norm
# clone s3norm package
# mkdir s3norm
# cd s3norm 
# git clone https://github.com/guanjue/S3norm.git

# run s3norm
## Entering working directory
cd $s3norm_working_directory
### Run S3norm
### old # time python $s3norm_script_directory'/src/s3norm_pipeline.py' -s $s3norm_script_directory'/src/' -t ${s3norm_working_directory}/nozeros_all_samples.txt
# s3norm_script_directory='/hpc/pmc_drost/nhung/S3norm'

time python $s3norm_script_directory'/src/s3norm_pipeline.py' -s $s3norm_script_directory'/src/' -t ${s3norm_working_directory}/${s3norm_sample_file_name}

## remove line with 0 reads in each sample to prepare for diffBind. data with 0 cause problem in deseq2  # Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
#   every gene contains at least one zero, cannot compute log geometric means # no need to run this. it can be removed in r
# export DATA_DIR_VARIABLE=$s3norm_working_directory/S3norm_NBP_bedgraph
# export AFTER_NORMALIZE_DIR_VARIABLE=$clean_normalize_dir

# find $s3norm_working_directory/S3norm_NBP_bedgraph  -name *.bedgraph.NBP.s3norm.bedgraph > $s3norm_working_directory/S3norm_NBP_bedgraph/sample_IDs.csv

# python clean_normalized_data.py