# Nhung 09 05 2023 
# The s3norm requires all files to have exactly the same first 3 columns (chr start end).
# This results in many chr with 0 reads in some samples. This creates problem for deseq2
# This script remove those lines before calling peaks and run diffbind
import pandas as pd
import os
# export DATA_DIR_VARIABLE=$s3norm_working_directory/S3norm_NBP_bedgraph
# export AFTER_NORMALIZE_DIR_VARIABLE=$clean_normalize_dir
#find $s3norm_working_directory/S3norm_NBP_bedgraph  -name *.bedgraph.NBP.s3norm.bedgraph > $s3norm_working_directory/S3norm_NBP_bedgraph/sample_IDs.csv

# sample_IDs = os.environ["SAMPLE_VARIABLE"]
data_dir = os.environ["DATA_DIR_VARIABLE"]
res_dir = os.environ["AFTER_NORMALIZE_DIR_VARIABLE"]
sample_IDs_df = pd.read_csv(data_dir + '/sample_IDs.csv',  header = None)
sample_IDs = sample_IDs_df[0].tolist()
samples_IDs_without_path = [sample.split('/')[-1] for sample in sample_IDs]
# data_dir=$s3norm_working_directory/S3norm_NBP_bedgraph
for sample, sample_without_path in zip(sample_IDs, samples_IDs_without_path):
    print (sample, sample_without_path)
# sample_IDs=$(find $s3norm_working_directory/S3norm_NBP_bedgraph  -name *.bedgraph.NBP.s3norm.bedgraph) 
# data_dir = '/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/modify_bedgraph/S3norm_NBP_bedgraph'
# sample_IDs = ['SCC-bulkChIC-PMC-DRO-024_nozeroes_augmented.bedgraph.NBP.s3norm.bedgraph', 'SCC-bulkChIC-PMC-DRO-025_nozeroes_augmented.bedgraph.NBP.s3norm.bedgraph']
# read bedgraph files as dataframes
# dfs = [pd.read_csv(data_dir + '/' + path, sep='\t', comment='t', header=None) for path in sample_IDs]

# for i in range(sample_IDs.size[0]):
#     dfs.append(pd.read_csv(sample_IDs[i][0], sep='\t', comment='t', header=None))


dfs = []
for sample in sample_IDs:
    dfs.append(pd.read_csv(sample, sep='\t', comment='t', header=None))

cleaned_dfs = [df[df[3]>0] for df in dfs]
print (cleaned_dfs[0])
for df, sample_without_path in zip(cleaned_dfs, samples_IDs_without_path):
   df.to_csv(res_dir + '/' + sample_without_path[:-9]+'_process_after_normalize.bedgraph', index=False, sep='\t', header=None)
