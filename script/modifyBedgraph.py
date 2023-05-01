# nhung 29 04 2023
import pandas as pd

if False:
    # path_to_bedgraph_file1 = '/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/bincount_bamcoverage_fixbinsize/bulkChIC-PMC-DRO-011_binsize_200.bedgraph'
    path_to_bedgraph_file1 = '/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/bincount_bamcoverage_fixbinsize/SCC-ChIC-PMC-DRO-T5_binsize_200.bedgraph'
    # path_to_bedgraph_file2 = '/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/bincount_bamcoverage_fixbinsize/bulkChIC-PMC-DRO-012_binsize_200.bedgraph'
    path_to_bedgraph_file2 = '/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/bincount_bamcoverage_fixbinsize/bulkChIC-PMC-DRO-013_binsize_200.bedgraph'
    path_to_bedgraph_file3 = '/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/bincount_bamcoverage_fixbinsize/SCC-ChIC-PMC-DRO-T1_binsize_200.bedgraph'

    df1 = pd.read_csv(path_to_bedgraph_file1, sep='\t', comment='t', header=None)
    df2 = pd.read_csv(path_to_bedgraph_file2, sep='\t', comment='t', header=None)
    df3 = pd.read_csv(path_to_bedgraph_file3, sep='\t', comment='t', header=None)

    df1['df2counts'] = df2[3]
    df1['df3counts'] = df3[3]
    df2['df1counts'] = df1[3]
    df2['df3counts'] = df3[3]
    df3['df1counts'] = df1[3]
    df3['df2counts'] = df2[3]

    df1 = df1[(df1[3]>0) & (df1['df2counts']>0) & (df1['df3counts']>0)].drop(columns=['df2counts', 'df3counts'])
    df2 = df2[(df2[3]>0) & (df2['df1counts']>0) & (df2['df3counts']>0)].drop(columns=['df1counts', 'df3counts'])
    df3 = df3[(df3[3]>0) & (df3['df1counts']>0) & (df3['df2counts']>0)].drop(columns=['df1counts', 'df2counts'])

    df1.to_csv(path_to_bedgraph_file1[:-9]+'_nozeroes.bedgraph', index=False, sep='\t', header=None)
    df2.to_csv(path_to_bedgraph_file2[:-9]+'_nozeroes.bedgraph', index=False, sep='\t', header=None)
    df3.to_csv(path_to_bedgraph_file3[:-9]+'_nozeroes.bedgraph', index=False, sep='\t', header=None)

directory = '/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/bincount_bamcoverage_fixbinsize/'
sample_list = ['SCC-ChIC-PMC-DRO-T5_binsize_200.bedgraph',
                'bulkChIC-PMC-DRO-013_binsize_200.bedgraph',
                'SCC-ChIC-PMC-DRO-T1_binsize_200.bedgraph']

#read bedgraph files as dataframes
dfs = [pd.read_csv(directory + path, sep='\t', comment='t', header=None) for path in sample_list]

sets = []
for df in dfs:
    #retrain only the rows where count==0. The fourth colmun ([3]) is the counts. From this result take the indexes and put them in a list. Make this list a set.
    sets.append(set(list(df[df[3]==0].index)))

#take the intersection of the sets in order to find the rows where all dataframes have count==0
intersection_set = set.intersection(*sets)
#convert the intersection set to a list
rows_for_removal  = list(intersection_set)
#drop the corresponding rows from all the dataframes
cleaned_dfs = [df.drop(index=rows_for_removal) for df in dfs]

def replace_zeroes_df(df, value_to_replace_zeroes):
    df[3] = df[3].apply(lambda v: value_to_replace_zeroes if v==0 else v)
    return df

def add_one(df):
    df[3] = df[3]+1
    return df

cleaned_augmented_dfs = [add_one(cleaned_df) for cleaned_df in cleaned_dfs]


#for df, cleaned_df, cleaned_augmented_df in zip(dfs, cleaned_dfs, cleaned_augmented_dfs):
#    print(df[[3]].describe())
#    print(cleaned_df[[3]].describe())
#    print(cleaned_augmented_df[[3]].describe())

#save the resulting dataframes as bedgraph files
for df, sample in zip(cleaned_augmented_dfs, sample_list):
    df.to_csv(directory + sample[:-9]+'_nozeroes_augmented_addone.bedgraph', index=False, sep='\t', header=None)

# s3norm samples
#df1 = pd.read_csv('/hpc/pmc_drost/nhung/S3norm/example_file/sig1.ctrl.sorted.bedgraph', sep='\t', comment='t', header=None)
#df1[[3]].describe()
#df2 = pd.read_csv('/hpc/pmc_drost/nhung/S3norm/example_file/sig2.ctrl.sorted.bedgraph', sep='\t', comment='t', header=None)
#df2[[3]].describe()



