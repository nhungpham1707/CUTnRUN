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
sample_IDs = pd.read_csv(data_dir + '/sample_IDs.csv',  header = None)
# data_dir=$s3norm_working_directory/S3norm_NBP_bedgraph
for i in range(len(sample_IDs)):
    print (sample_IDs.iloc[i][0])
# sample_IDs=$(find $s3norm_working_directory/S3norm_NBP_bedgraph  -name *.bedgraph.NBP.s3norm.bedgraph) 
# data_dir = '/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/modify_bedgraph/S3norm_NBP_bedgraph'
# sample_IDs = ['SCC-bulkChIC-PMC-DRO-024_nozeroes_augmented.bedgraph.NBP.s3norm.bedgraph', 'SCC-bulkChIC-PMC-DRO-025_nozeroes_augmented.bedgraph.NBP.s3norm.bedgraph']
# read bedgraph files as dataframes
# dfs = [pd.read_csv(data_dir + '/' + path, sep='\t', comment='t', header=None) for path in sample_IDs]

# for i in range(sample_IDs.size[0]):
#     dfs.append(pd.read_csv(sample_IDs[i][0], sep='\t', comment='t', header=None))


dfs = []
for i in range(len(sample_IDs)):
    dfs.append(pd.read_csv(sample_IDs.iloc[i][0], sep='\t', comment='t', header=None))

sets = []
for df in dfs:
    #retrain only the rows where count==0. The fourth colmun ([3]) is the counts. From this result take the indexes and put them in a list. Make this list a set.
    sets.append(set(list(df[df[3]==0].index)))

sets[0]

cleaned_dfs = []
# for i in range(len(sample_IDs)):
#     cleaned_dfs.append(dfs[i].drop(index=sets[i])])

for df, my_set in zip(dfs, sets):
    cleaned_dfs.append(df.drop(index=my_set))

cleaned_dfs[0]

for df, sample in zip(cleaned_dfs, sample_IDs):
    df.to_csv(res_dir + sample[:-9]+'_process_after_normalize.bedgraph', index=False, sep='\t', header=None)
