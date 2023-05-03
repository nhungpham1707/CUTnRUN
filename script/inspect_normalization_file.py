import pandas as pd
import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

#our samples
#sample 1
df1 = pd.read_csv('/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/bincount_bamcoverage_fixbinsize/SCC-ChIC-PMC-DRO-T5_binsize_200_nozeroes_augmented_addone.bedgraph', sep='\t', comment='t', header=None)
df2 = pd.read_csv('/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/bincount_bamcoverage_fixbinsize/S3norm_NBP_bedgraph/SCC-ChIC-PMC-DRO-T5_binsize_200_nozeroes_augmented_addone.bedgraph.NBP.s3norm.bedgraph', sep='\t', comment='t', header=None)
#sample 2
df3 = pd.read_csv('/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/bincount_bamcoverage_fixbinsize/SCC-ChIC-PMC-DRO-T1_binsize_200_nozeroes_augmented_addone.bedgraph', sep='\t', comment='t', header=None)
df4 = pd.read_csv('/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/bincount_bamcoverage_fixbinsize/S3norm_NBP_bedgraph/SCC-ChIC-PMC-DRO-T1_binsize_200_nozeroes_augmented_addone.bedgraph.NBP.s3norm.bedgraph', sep='\t', comment='t', header=None)


#package samples
#sample 3
df5 = pd.read_csv('/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/s3norm_example_file/sig3.sorted.bedgraph', sep='\t', comment='t', header=None)
df6 = pd.read_csv('/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/s3norm_example_file/S3norm_NBP_bedgraph/sig3.sorted.bedgraph.NBP.s3norm.bedgraph', sep='\t', comment='t', header=None)
# sample 4
df7 = pd.read_csv('/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/s3norm_example_file/sig2.sorted.bedgraph', sep='\t', comment='t', header=None)
df8 = pd.read_csv('/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/s3norm_example_file/S3norm_NBP_bedgraph/sig2.sorted.bedgraph.NBP.s3norm.bedgraph', sep='\t', comment='t', header=None)

fig, axes = plt.subplots(4, 2, figsize=(14,8))
plt.tight_layout()
binwidth = 20

min_ = min(df1[3]) if min(df1[3])<min(df3[3]) else min(df3[3])
max_ = max(df1[3]) if max(df1[3])>max(df3[3]) else max(df3[3])
#sample 1
df1[3].plot.hist(ax=axes[0,0], bins=np.arange(min_, max_ + binwidth, binwidth), log=True, ylim=(0,10**7)).set_title('original_counts')
df2[3].plot.hist(ax=axes[0,1], bins=np.arange(min_, max_ + binwidth, binwidth), log=True, ylim=(0,10**7)).set_title('normalized')
#sample 2
df3[3].plot.hist(ax=axes[1,0], bins=np.arange(min_, max_ + binwidth, binwidth), log=True, ylim=(0,10**7)).set_title('original_counts')
df4[3].plot.hist(ax=axes[1,1], bins=np.arange(min_, max_ + binwidth, binwidth), log=True, ylim=(0,10**7)).set_title('normalized')

min_ = min(df5[3]) if min(df5[3])<min(df7[3]) else min(df7[3])
max_ = max(df5[3]) if max(df5[3])>max(df7[3]) else max(df7[3])
binwidth = 1000
#sample 3
df5[3].plot.hist(ax=axes[2,0], bins=np.arange(min_, max_ + binwidth, binwidth), log=True, ylim=(0,10**4)).set_title('original_counts')
df6[3].plot.hist(ax=axes[2,1], bins=np.arange(min_, max_ + binwidth, binwidth), log=True, ylim=(0,10**4)).set_title('normalized')
#sample 4
df7[3].plot.hist(ax=axes[3,0], bins=np.arange(min_, max_ + binwidth, binwidth), log=True, ylim=(0,10**4)).set_title('original_counts')
df8[3].plot.hist(ax=axes[3,1], bins=np.arange(min_, max_ + binwidth, binwidth), log=True, ylim=(0,10**4)).set_title('normalized')
plt.savefig('/home/pmc_research/npham/PROJECTS/CUTnRUN_Maroussia/Figures/histograms.png') 

df1 = pd.read_csv('/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/bincount_bamcoverage_fixbinsize/SCC-ChIC-PMC-DRO-T5_binsize_200_nozeroes_augmented_addone.bedgraph', sep='\t', comment='t', header=None)
df1[3].plot.hist(ax=axes[1,0], bins=np.arange(min_, max_ + binwidth, binwidth), log=True, ylim=(0,10**7)).set_title('original_counts')
plt.savefig('/home/pmc_research/npham/PROJECTS/CUTnRUN_Maroussia/Figures/histograms.png')