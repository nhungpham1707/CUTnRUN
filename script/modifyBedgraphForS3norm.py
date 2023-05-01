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

directory = '/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/bincount_bamcoverage_fixbinsize/'
out_directory= '/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/modify_bedgraph/'

sample_list= ["bulkChIC-PMC-DRO-011",
            "bulkChIC-PMC-DRO-012",
            "bulkChIC-PMC-DRO-013",
            "bulkChIC-PMC-DRO-014",
            "bulkChIC-PMC-DRO-015",
            "bulkChIC-PMC-DRO-016",
            "SCC-bulkChIC-PMC-DRO-002",
            "SCC-bulkChIC-PMC-DRO-005",
            "SCC-bulkChIC-PMC-DRO-008",
            "SCC-ChIC-PMC-DRO-L5",
            "SCC-ChIC-PMC-DRO-LH",
            "SCC-ChIC-PMC-DRO-F1",
            "SCC-ChIC-PMC-DRO-F5",
            "SCC-ChIC-PMC-DRO-FH",
            "SCC-ChIC-PMC-DRO-T1",
            "SCC-ChIC-PMC-DRO-T5",
            "SCC-ChIC-PMC-DRO-TH", 
            "SCC-ChIC-PMC-DRO-L1",
            "SCC-bulkChIC-PMC-DRO-020",
            "SCC-bulkChIC-PMC-DRO-021",
            "SCC-bulkChIC-PMC-DRO-022",
            "SCC-bulkChIC-PMC-DRO-023",
            "SCC-bulkChIC-PMC-DRO-024",
            "SCC-bulkChIC-PMC-DRO-025"]

#read bedgraph files as dataframes
dfs = [pd.read_csv(directory + path+ '_binsize_200.bedgraph', sep='\t', comment='t', header=None) for path in sample_list]

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

def replace_zeroes_df(df, value_to_replace_zeroes):
    df[3] = df[3].apply(lambda v: value_to_replace_zeroes if v==0 else v)
    return df

def add_one(df):
    df[3] = df[3]+1
    return df

# cleaned_augmented_dfs = [add_one(cleaned_df) for cleaned_df in cleaned_dfs]

cleaned_augmented_dfs= [replace_zeroes_df(cleaned_df,1) for cleaned_df in cleaned_dfs]
#for df, cleaned_df, cleaned_augmented_df in zip(dfs, cleaned_dfs, cleaned_augmented_dfs):
#    print(df[[3]].describe())
#    print(cleaned_df[[3]].describe())
#    print(cleaned_augmented_df[[3]].describe())

#save the resulting dataframes as bedgraph files
for df, sample in zip(cleaned_augmented_dfs, sample_list):
    df.to_csv(out_directory + sample + '_nozeroes_augmented.bedgraph', index=False, sep='\t', header=None)

# s3norm samples
#df1 = pd.read_csv('/hpc/pmc_drost/nhung/S3norm/example_file/sig1.ctrl.sorted.bedgraph', sep='\t', comment='t', header=None)
#df1[[3]].describe()
#df2 = pd.read_csv('/hpc/pmc_drost/nhung/S3norm/example_file/sig2.ctrl.sorted.bedgraph', sep='\t', comment='t', header=None)
#df2[[3]].describe()



