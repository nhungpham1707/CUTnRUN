# Nhung 09 05 2023
# modify sample sheet file for diffBind (change peak files)
import pandas as pd
import glob

sample_sheet = pd.read_csv('/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/R/diffbind_normalize_samples.csv', sep = ';')

path_to_new_peak = '/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/normalize_peakCalling/'

for i in range(len(sample_sheet["SampleID"])):
    for name in glob.glob(path_to_new_peak + sample_sheet["SampleID"][i]+ '*'):
        print name
        sample_sheet["Peaks"][i] = name


