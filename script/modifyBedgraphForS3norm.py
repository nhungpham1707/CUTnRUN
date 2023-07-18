# prepare bedgraph files for s3norm normalization
# Input: bedgraph files with equal bin sizes for all samples. The 1st 3 columns in all samples are identical. 
# These input files are generated with bamcoverage (bincount_bamcoverage.sh)
# Steps:
# - read bedgraph files as dataframe
# - remove rows that are 0 in all samples
# - modify 0 value in samples. Replace them with 1 to prevent inf log transform
# - write to bedgraph files 
# output: bedgraph files that are processed and are ready to use as input for s3norm

# nhung 01 05 2023
# conda env s3norm
import pandas as pd
import os
from pathlib import Path

def replace_zeroes_df(df, value_to_replace_zeroes):
    df[3] = df[3].apply(lambda v: value_to_replace_zeroes if v==0 else v)
    return df

def add_one(df):
    df[3] = df[3]+1
    return df

# get variables from config.sh
directory = os.environ["modify_bedgraph_dir"]
bin_size = os.environ["bin_size"]

s3norm_sample_list = os.getenv('s3norm_sample_list').split(',')


paths = [Path(directory + '/' + sample + 'bin_size' + bin_size +  '_nozeroes_augmented.bedgraph') for sample in s3norm_sample_list]

new_sample_list = []
for path, sample in zip(paths, s3norm_sample_list):
    if path.is_file() == False:
        new_sample_list.append(sample)

#if list is empty
if not new_sample_list:
    print ('all files exist. Skip modify bedgraph !')
else:

    #read bedgraph files as dataframes
    dfs = [pd.read_csv(directory + '/' + sample + '_binsize_' + bin_size + '.bedgraph', sep='\t', comment='t', header=None) for sample in s3norm_sample_list]

    sets = []
    for df in dfs:
        #retrain only the rows where count==0. The fourth colmun ([3]) is the counts. From this result take the indexes and put them in a list. Make this list a set.
        sets.append(set(list(df[df[3]==0].index)))

    #take the intersection of the sets in order to find the rows where all dataframes have count==0
    intersection_set = set.intersection(*sets)  # check if all intersect are 0 or a nonzeros value
    #convert the intersection set to a list
    rows_for_removal  = list(intersection_set)
    #drop the corresponding rows from all the dataframes
    cleaned_dfs = [df.drop(index=rows_for_removal) for df in dfs]

    cleaned_augmented_dfs= [replace_zeroes_df(cleaned_df,1) for cleaned_df in cleaned_dfs]

    #save the resulting dataframes as bedgraph files
    for df, sample in zip(cleaned_augmented_dfs, s3norm_sample_list):
        df.to_csv(directory + '/' + sample + 'bin_size' + bin_size +  '_nozeroes_augmented.bedgraph', index=False, sep='\t', header=None)





