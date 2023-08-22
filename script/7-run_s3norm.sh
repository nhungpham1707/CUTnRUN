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


DIR=${modify_bedgraph_dir}/$s3norm_sample_file_name/S3norm_NBP_bedgraph
# init
# look for empty dir
if [ -d "$DIR" ]
then
	if [ "$(ls -A $DIR)" ]; then
     echo "$DIR is not Empty. Skipped s3norm normalization for these samples"
	else
    time python $s3norm_script_directory'/src/s3norm_pipeline.py' -s $s3norm_script_directory'/src/' -t ${s3norm_working_directory}/${s3norm_sample_file_name}.csv -r mean -m mean -i 2.0 -f 0.05 -l 0.001 -a 100000 -b 0 -p z -k 0 -g 0

    mkdir -p $s3norm_sample_file_name
    mv $modify_bedgraph_dir/average_ref_bedgraph $modify_bedgraph_dir/$s3norm_sample_file_name/
    mv $modify_bedgraph_dir/NBP_bedgraph $modify_bedgraph_dir/$s3norm_sample_file_name/
    mv $modify_bedgraph_dir/S3norm_NBP_bedgraph $modify_bedgraph_dir/$s3norm_sample_file_name/
    mv $modify_bedgraph_dir/S3norm_rc_bedgraph $modify_bedgraph_dir/$s3norm_sample_file_name/
	fi
else
	echo "Directory $DIR not found."
fi

