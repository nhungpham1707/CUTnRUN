# Nhung 03 05 2023
# ref https://bioinformatics-core-shared-training.github.io/Quantitative-ChIPseq-Workshop/articles/Quantitative-ChIPseq-Workshop.html
# ref https://bioconductor.org/packages/devel/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("DiffBind")

library(DiffBind)
library(dplyr)

source('/home/pmc_research/npham/PROJECTS/CUTnRUN_Maroussia/script/Identify_DE_peaks_functions.R')
# to test on terminal
# sample_sheet_name_dir <- '/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/R/diffBind_sample_sheet_s3norm_data_no_control.csv'
# save_name <- 'with_no_control_hpc6'
# res_dir <- '/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/diffBind_analysis/with_no_control_hpc6'
# get paths from environment
res_dir <- Sys.getenv("DIFFBIND_RESULT_DIR_VARIABLE")
save_name <- Sys.getenv("SAVE_NAME_VARIABLE")
message(paste('result will be saved at ', res_dir))
sample_sheet_name_dir <- Sys.getenv("SAMPLE_SHEET_DIR_VARIABLE")
samples <- read.csv(sample_sheet_name_dir)
head(samples)

# create dba object for further analysis
print ("start peak")
peaks <- dba(sampleSheet=samples)
# remove peaks in the problematic regions
print ("start applying black list")
peaks  <- dba.blacklist(peaks, greylist = FALSE)    
## Consensus peaks sets and counting reads
olap.rate <- dba.overlap(peaks, mode=DBA_OLAP_RATE)
reso <- 600
png(paste0(res_dir,'/', Sys.Date(), '-', save_name ,'-overlapping_peaks.png'), width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
plot(olap.rate, xlab="Overlapping samples", ylab="Overlapping peaks", type="l", ylim = c(0, max(olap.rate)))
text(x = seq(1,length(olap.rate), by=1) , y = olap.rate - 500, labels=as.character(olap.rate))
dev.off()
# consensus.peaks <- dba.peakset(peaks_narrow, bRetrieve=TRUE) # extract peaks that are in at least 2 samples
# consensus.peaks[,0] # 47 sequences 
# counting overlapping peaks
peak_counts <- dba.count(peaks)
save.image(paste0(res_dir, '/', Sys.Date(),'-', save_name,"-peak_count_DiffBind.RData"))   
#
# message('load Rdata')
# load('/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/diffBind_analysis/remove_low_depth_histone_samples/2023-05-21-remove_low_depth_histone_samples-peak_count_DiffBind.RData')   
# plot correlation 
FDR_thredshold <- 0.05
fold_change_thredshold <- 2
contrast_number <- 3
column_to_keep = c('Conc_luciferase','Conc_fusion')
peak_model_test <-  conduct_diffBind_analysis(peaks, peak_counts, res_dir, save_name,column_to_keep = column_to_keep, contrast_number = 3)
# write to bed file without using the function
res_deseq_contrast <- dba.report(peak_model_test, method=DBA_DESEQ2, contrast = 3, th=1)
res_deseq_contrast_df <- as.data.frame(res_deseq_contrast)

bed_contrast <- res_deseq_contrast_df %>% 
  dplyr::filter(FDR < FDR_thredshold & abs(Fold) > fold_change_thredshold ) %>% 
  dplyr::select(seqnames, start, end, column_to_keep, Fold, p.value, FDR)
message(paste ('the number of DE sites is:', length(bed_contrast$seqnames) ))
bed_constrast <- bed_contrast[order(-abs(bed_contrast$Fold)),]
# plot(bed_contrast$Fold)
write.table(bed_contrast, file=paste0(res_dir, '/', Sys.Date() ,save_name,
            '_fold_change', fold_change_thredshold, 
            '_FDR',FDR_thredshold ,'.bed'),
            sep="\t", quote=F, row.names=F, col.names=F)

reso <- 1200
png(filename = paste0(res_dir,'/', Sys.Date() ,'-', save_name,'-diffbind_correlation.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
plot(peak_counts, mode=DBA_ID)
dev.off()

message('make and write report')
report <- dba.show(peak_model_test, bContrasts = TRUE) # the last column show the number of differential peaks in each contrast  using the default threshold of FDR <= 0.05
write.csv(report,(paste0(res_dir,'/', Sys.Date(),'-', save_name,'-dba-show-report.csv')))


# png(filename = paste0(res_dir, '/', Sys.Date() ,'-', save_name,'-diffbind_Volcano_contrast3.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
# dba.plotVolcano(peak_model, contrast = 3, dotSize = 2, bFlip = TRUE)
# dev.off()

# write table result
message('make and write result table')
result <- data.frame(total_peak = peaks$totalMerged,DE_contrast_1 = as.numeric(report$DB.DESeq2[1]), DE_contrast_2 = as.numeric(report$DB.DESeq2[2]), DE_contrast_3 = as.numeric(report$DB.DESeq2[3]))

write.csv(result, (paste0(res_dir,'/', Sys.Date(),'-', save_name,'-Number_of_differential_peak_per_contrast.csv')))

save.image(paste0(res_dir, '/', Sys.Date(),'-', save_name,"-DiffBind-peak-model.RData"))     

png(filename = paste0(res_dir, '/', Sys.Date(), "_", save_name, "_volcano_plot.png") , width = 1200 * reso/72, 
    height = 700 * reso/72, units ="px", res = reso)
dba.plotVolcano(peak_model_test, contrast = contrast_number, dotSize = 2)  # bFlip = TRUE
dev.off()

message('finished!')