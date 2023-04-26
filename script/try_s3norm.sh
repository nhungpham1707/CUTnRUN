#!/bin/bash
#SBATCH --job-name=tabtestS3
#SBATCH --output=test_s3norm_tab_file.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --gres=tmpspace:30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl



# normalize data using s3norm package. run with their file_list.txt and bedgraph files to check if the code work
# nhung 25 04 2023
# conda env s3norm
# clone s3norm package
# mkdir s3norm
# cd s3norm 
# git clone https://github.com/guanjue/S3norm.git

# run s3norm
### Setting script directory
# script_directory='/home/pmc_research/npham/PROJECTS/CUTnRUN_Maroussia/script/s3norm/S3norm'

script_directory='/hpc/pmc_drost/nhung/S3norm'
### Setting working directory

working_directory=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test

## Entering working directory
cd $working_directory
### Run S3norm
# time python $script_directory'/src/s3norm_pipeline.py' -s $script_directory'/src/' -t s3norm_1sample.csv

# time python $script_directory'/src/s3norm_pipeline.py' -s $script_directory'/src/' -t ${working_directory}/s3norm_file_list.txt

time python $script_directory'/src/s3norm_pipeline.py' -s $script_directory'/src/' -t ${working_directory}/s3norm_tab_file_list.txt