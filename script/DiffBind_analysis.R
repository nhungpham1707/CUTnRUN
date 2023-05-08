# Nhung 03 05 2023
# ref https://bioinformatics-core-shared-training.github.io/Quantitative-ChIPseq-Workshop/articles/Quantitative-ChIPseq-Workshop.html
# ref https://bioconductor.org/packages/devel/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("DiffBind")

library(DiffBind)
# library(diffbind)
# library(Gviz)
# library(profileplyr)
library(dplyr)
# library(ChIPpeakAnno)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(EnsDb.Hsapiens.v86)

# get paths from environment
res_dir <- Sys.getenv("DIFFBIND_RESULT_DIR_VARIABLE")
sample_sheet_name_dir <- Sys.getenv("SAMPLE_SHEET_DIR_VARIABLE")
samples <- read.csv(sample_sheet_name_dir, sep = ';')
head(samples)

try({

# create dba object for further analysis
print ("start peak")
peaks <- dba(sampleSheet=samples)
# remove peaks in the problematic regions
print ("start applying black list")
peaks  <- dba.blacklist(peaks, greylist = FALSE)    
## Consensus peaks sets and counting reads
olap.rate <- dba.overlap(peaks, mode=DBA_OLAP_RATE)
reso <- 600
png(paste0(res_dir,'/', Sys.Date() ,'-overlapping_peaks.png'), width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
plot(olap.rate, xlab="Overlapping samples", ylab="Overlapping peaks", type="l", ylim = c(0, max(olap.rate)))
text(x = seq(1,length(olap.rate), by=1) , y = olap.rate - 500, labels=as.character(olap.rate))
dev.off()
#consensus.peaks <- dba.peakset(peaks_narrow, bRetrieve=TRUE) # extract peaks that are in at least 2 samples
#consensus.peaks[,0] # 47 sequences 
# counting overlapping peaks
peak_counts <- dba.count(peaks)
save.image(paste0(Sys.Date(),"-peak_count_DiffBind.RData"))      
# plot correlation 

reso <- 1200
png(filename = paste0(res_dir,'/', Sys.Date() ,'-diffbind_correlation.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
plot(peak_counts, mode=DBA_ID)
dev.off()

# plot pca
png(filename = paste0(res_dir,'/', Sys.Date() ,'-diffbind_PCA.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
dba.plotPCA(peak_counts, label = DBA_TREATMENT)
dev.off()

# plot contrast
# png(filename = paste0(res_dir, "/diffbind_contrast.png") , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
#               dba.plotMA(counts_narrow, bNormalized=FALSE,
#               contrast=list(fusion=counts_narrow$masks$fusion,
#                                 tfe3=counts_narrow$masks$tfe3))
# dev.off()

# generate differential analysis model 
print ("start peak model")       
peak_model <- dba.contrast(peak_counts, reorderMeta=list(Condition="tfe3"), design = "~Condition")
peak_model <- dba.analyze(peak_model,bParallel=FALSE, bGreylist = FALSE)
report <- dba.show(peak_model, bContrasts = TRUE) # the last column show the number of differential peaks in each contrast  using the default threshold of FDR <= 0.05

# plotMA
png(filename = paste0(res_dir, '/', Sys.Date() ,'-diffbind_plotMA_contrast1.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
dba.plotMA(peak_model,contrast=1)
dev.off()

png(filename = paste0(res_dir, '/', Sys.Date() ,'-diffbind_plotMA_contrast2.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
dba.plotMA(peak_model,contrast=2)
dev.off()
       
png(filename = paste0(res_dir, '/', Sys.Date() ,'-diffbind_plotMA_contrast3.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
dba.plotMA(peak_model,contrast=3)
dev.off()

# plot vocanol      
png(filename = paste0(res_dir, '/', Sys.Date() ,'-diffbind_Volcano_contrast1.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
dba.plotVolcano(peak_model, contrast = 1)
dev.off()

png(filename = paste0(res_dir, '/', Sys.Date() ,'-diffbind_Volcano_contrast2.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
dba.plotVolcano(peak_model, contrast = 2)
dev.off()

png(filename = paste0(res_dir, '/', Sys.Date() ,'-diffbind_Volcano_contrast3.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
dba.plotVolcano(peak_model, contrast = 3)
dev.off()

#dba.show(counts_narrow,attributes=c(DBA_ID,DBA_FRIP))

# Venn diagram of Gain vs Loss differentially bound sites. Gain" sites (those that
# increase binding enrichment in the Resistant condition) and the "Loss" sites (those with lower enrichment)
png(filename = paste0(res_dir, '/', Sys.Date() ,'-diffbind_Venn_contrast1.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
dba.plotVenn(peak_model, contrast = 1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)
dev.off()
       
png(filename = paste0(res_dir, '/', Sys.Date() ,'-diffbind_Venn_contrast2.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
dba.plotVenn(peak_model, contrast = 2, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)
dev.off()
       
png(filename = paste0(res_dir, '/', Sys.Date() ,'-diffbind_Venn_contrast3.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
dba.plotVenn(peak_model, contrast = 3, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)
dev.off()
       
# A PCA plot using only the 246 differentially bound sites (corresponding to Figure 3), using
# an FDR threshold of 0.05
png(filename = paste0(res_dir, '/', Sys.Date() ,'-diffbind_PCA_DE_contrast1.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
dba.plotPCA(peak_model, contrast = 1)
dev.off()

png(filename = paste0(res_dir, '/', Sys.Date() ,'-diffbind_PCA_DE_contrast2.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
dba.plotPCA(peak_model, contrast = 2)
dev.off()

png(filename = paste0(res_dir, '/', Sys.Date() ,'-diffbind_PCA_DE_contrast3.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
dba.plotPCA(peak_model, contrast = 3)
dev.off()       
       
# write table result 
result <- data.frame(total_peak = counts_narrow$totalMerged,
                     DE_contrast_1 = as.numeric(report$DB.DESeq2[1]) ,
                     DE_contrast_2 = as.numeric(report$DB.DESeq2[2]) ,
                     DE_contrast_3 = as.numeric(report$DB.DESeq2[3]))

write.csv(result, (paste0(res_dir,'/', Sys.Date(),'-Number_of_differential_peak_per_contrast.csv'))
# Plot heatmap profile
profiles_all <- dba.plotProfile(peak_model, merge=c(DBA_FACTOR, DBA_REPLICATE, DBA_TREATMENT))
png(filename = paste0(res_dir, '/', Sys.Date() ,'-diffbind_profile.png') , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
dba.plotProfile(profiles_all)
dev.off()

## write DE binding site result
# Create bed files for each keeping only significant peaks (p < 0.05)
res_deseq_contrast1 <- dba.report(peak_model, method=DBA_DESEQ2, contrast = 1, th=1)
res_deseq_contrast1_df <- as.data.frame(res_deseq_contrast1)
bed_contrast1 <- res_deseq_contrast1_df %>% 
        dplyr::filter(FDR > 0.05 & abs(Fold)  > 2) %>% 
        dplyr::select(seqnames, start, end)
write.table(bed_contrast1, file=paste0(res_dir,'/', Sys.Date() ,'-diffBind_contrast1.bed'), sep="\t", quote=F, row.names=F, col.names=F)

res_deseq_contrast2 <- dba.report(peak_model, method=DBA_DESEQ2, contrast = 2, th=1)
res_deseq_contrast2_df <- as.data.frame(res_deseq_contrast2)
bed_contrast2 <- res_deseq_contrast2_df %>% 
        dplyr::filter(FDR > 0.05 & abs(Fold)  > 2) %>% 
        dplyr::select(seqnames, start, end)
write.table(bed_contrast2, file=paste0(res_dir,'/', Sys.Date() ,'-diffBind_contrast2.bed'), sep="\t", quote=F, row.names=F, col.names=F)

res_deseq_contrast3 <- dba.report(peak_model, method=DBA_DESEQ2, contrast = 3, th=1)
res_deseq_contrast3_df <- as.data.frame(res_deseq_contrast3)
bed_contrast3 <- res_deseq_contrast3_df %>% 
        dplyr::filter(FDR > 0.05 & abs(Fold)  > 2) %>% 
        dplyr::select(seqnames, start, end)
write.table(bed_contrast3, file=paste0(res_dir,'/', Sys.Date() ,'-diffBind_contrast3.bed'), sep="\t", quote=F, row.names=F, col.names=F)

})
save.image(paste0(Sys.Date(),"-DiffBind.RData"))      