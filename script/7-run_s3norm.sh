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
# cd $s3norm_working_directory
### Run S3norm
# time python $s3norm_script_directory'/src/s3norm_pipeline.py' -s $s3norm_script_directory'/src/' -t ${s3norm_working_directory}/nozeros_all_samples.txt

# time python $s3norm_script_directory'/src/s3norm_pipeline.py' -s $s3norm_script_directory'/src/' -t ${s3norm_working_directory}/no_zeros_testsamples.txt
# s3norm_script_directory='/hpc/pmc_drost/nhung/S3norm'

s3norm_working_directory_test=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/test_FH/

cd $s3norm_working_directory_test

time python $s3norm_script_directory'/src/s3norm_pipeline.py' -s $s3norm_script_directory'/src/' -t ${s3norm_working_directory_test}/test_FH_sample.txt