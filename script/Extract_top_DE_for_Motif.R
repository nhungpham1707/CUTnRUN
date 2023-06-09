

cat 2023-06-01antiTFE3_Samples_fold_change2_FDR0.05.bed | awk 'NR >= 0  && NR <= 2001 { print }'  > top_2000_DE_Sites.bed

setwd('/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/diffBind_analysis/antiTFE3_Samples/')
all_peak <- read.csv('2023-06-01antiTFE3_Samples_fold_change2_FDR0.05.bed', 
            sep = '\t', header = FALSE)
head(all_peak)
length(all_peak[,1])
top_sites <- all_peak[1:2000,]
write.table(top_sites, 'top_2000_DE_sites.bed', sep = '\t', row.names = FALSE, col.names = FALSE)
