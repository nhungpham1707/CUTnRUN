# Module 2: Data normalization and peak calling after normalization 

# generate bedgraph files with a fix bin size for all samples to use as input for s3norm 
. ./make_fix_bin_size_bedgraph.sh

# modify bedgraph file to remove rows with 0 count in all samples. If 0 values are more than 10% in the sample, s3norm will fail (log0 is inf), hence add 1 to these 0 values in all samples

python modifyBedgraphForS3norm.py 

# start normalization on the modified bedgraph files

. ./7-run_s3norm.sh 

