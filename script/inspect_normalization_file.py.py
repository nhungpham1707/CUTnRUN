import pandas as pd
df2 = pd.read_csv('/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/bincount_bamcoverage_fixbinsize/S3norm_NBP_bedgraph/SCC-ChIC-PMC-DRO-T5_binsize_200_nozeroes_augmented.bedgraph.NBP.s3norm.bedgraph', sep='\t', comment='t', header=None)
df1 = pd.read_csv('/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/bincount_bamcoverage_fixbinsize/SCC-ChIC-PMC-DRO-T5_binsize_200_nozeroes_augmented.bedgraph', sep='\t', comment='t', header=None)
df1['normalized']=df2[3]
df1['diffs']=df1[3]-df1['normalized']
df1.sort_values('diffs')
# df.to_csv(directory + sample[:-9]+'_nozeroes_augmented.bedgraphe)