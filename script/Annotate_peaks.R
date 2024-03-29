# This script annotate location, find enrich GO term,
# and annotate genes with ensembl and symbol names for DE sites 

# Nhung 25 05 2023


# source('/home/pmc_research/npham/PROJECTS/CUTnRUN_Maroussia/script/Identify_DE_peaks_functions.R')

# library(DiffBind)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(ChIPseeker)
# library(clusterProfiler)
library(AnnotationDbi)
# library(org.Hs.eg.db)
library(ggplot2)
library(viridis)
library(ggpubr)

res_dir <- Sys.getenv("DIFFBIND_RESULT_DIR_VARIABLE")
save_name <- Sys.getenv("SAVE_NAME_VARIABLE")
samplefiles <- Sys.getenv("SAMPLE_BED_FILE")
################## Annotate significant peaks 
# location annotation
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
edb <- EnsDb.Hsapiens.v86

samplefiles <- as.list(samplefiles)
names(samplefiles) <- "luc_vs_fusion"
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE) # this step takes a while

save.image(paste0(res_dir,'/','Peak_AnnoList.RData'))

# plot the distribution
annotate_df_luc_vs_fusion <- c()
for (i in 1:length(peakAnnoList[[1]]@anno$annotation)){
  annotate_df_luc_vs_fusion = data.frame(annotate = strsplit(peakAnnoList[[1]]@anno$annotation[i]," +")[[1]][1]) %>% rbind(annotate_df_luc_vs_fusion, .)
}

annotate_table_luc_vs_fusion<- as.data.frame(table(annotate_df_luc_vs_fusion));annotate_table_luc_vs_fusion


##### plot chromosome coverage
peakFusion <- readPeakFile(samplefiles[[1]])

#### Gene annotation 
fusionLuc_annot <- data.frame(peakAnnoList[["luc_vs_fusion"]]@anno)
fusionLuc_annot_entrez <- fusionLuc_annot$geneId
# Return the gene symbol for the set of Entrez IDs
fusionluc_annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                                   keys = fusionLuc_annot_entrez,
                                                   columns = c("GENENAME"),
                                                   keytype = "ENTREZID")
# conduct GO enrichment analysis
ego_fuLuc <- enrichGO(gene = fusionluc_annotations_edb$ENTREZID, 
                      keyType = "ENTREZID", 
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = TRUE)

save.image(paste0(res_dir, '/', "annotation_enrichment.RData"))

cluster_summary_luc <- data.frame(ego_fuLuc)
write.csv(cluster_summary_luc, paste0(res_dir, '/', Sys.Date(),"GO_enrichment.csv"))

## ensembl annotation
## get the txdb and peak files
files <- samplefiles
peak <- ChIPseeker::readPeakFile(files[[1]])

## change the seqlevels of peak files
oldseqlevels <- data.frame(old=GenomeInfoDb::seqlevels(peak),
                           new=gsub("chr","",GenomeInfoDb::seqlevels(peak)))
newseqlevels <- data.frame(new=c(1:22,"X","Y"))
combine <- merge(oldseqlevels,newseqlevels,by="new")
newnames <- combine$new
names(newnames) <- combine$old
peak <- GenomeInfoDb::renameSeqlevels(peak,newnames)
annotatePeak_result <- ChIPseeker::annotatePeak(peak = peak,TxDb = edb) # take time 
LucFusion_annotEDB <- data.frame(annotatePeak_result@anno)
LucFusion_annotEDB
colnames(LucFusion_annotEDB)[6:10] = c('Conc_luciferase','Conc_fusion', 'Fold', 'p.value', 'FDR')

# Get gene symbols
dataset <- LucFusion_annotEDB
ensembl <- dataset$geneId

# Return the gene symbol for the set of ensemble ids
annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                         keys = ensembl,
                                         columns = c("GENENAME"),
                                         keytype = "GENEID")
annotations_edb$GENEID <- as.character(annotations_edb$GENEID)
# Write to file
peak_annotation <- dataset %>% 
  left_join(annotations_edb, by=c("geneId"="GENEID"))

length(peak_annotation$seqnames) #4181
write.table(peak_annotation, file=paste0(res_dir, '/', Sys.Date(),"functional_analysis_EDB.txt"), sep="\t", quote=F, row.names=F)

save.image(paste0(res_dir, '/',"s3norm_all_sig_sites_GO_enrich_ensembl_annotation.RData"))

# plot
png(filename = paste0(res_dir, '/', Sys.Date(), "-chromosome_coverage.png") , width = 1400 * reso/72, height = 700 * reso/72, units ="px", res = reso)
covplot(peakFusion, weightCol="V5", chrs=c(paste0("chr", c(1:22, "X"))))
dev.off

png(filename = paste0(res_dir, '/', Sys.Date(), "-GO_enrichment.png") , width = 1400 * reso/72, height = 700 * reso/72, units ="px", res = reso)

plot = dotplot(ego_fuLuc, showCategory=15, font.size = 30) 
annotate_figure(plot, top = text_grob("GO enrichment fusion vs luciferase", 
                                      color = "black", face = "bold", size = 30))
dev.off()

png(filename = paste0(res_dir, '/', Sys.Date(), "-chromosome_distribution.png") , width = 1400 * reso/72, height = 700 * reso/72, units ="px", res = reso)
plot(LucFusion_annotEDB$seqnames)
dev.off()

reso <- 600
png(filename = paste0(res_dir, '/', Sys.Date(), "-annotation.png") , width = 1400 * reso/72, height = 700 * reso/72, units ="px", res = reso)
annotate_table_luc_vs_fusion %>% ggplot(aes(x = annotate, y = Freq, fill = annotate)) +
  geom_bar(stat="identity") + theme_bw(base_size = 24) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +   scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +   xlab("") + ylab("Frequency") + ggtitle( "A. Annotation of DE binding sites (fusion vs luciferase)") + theme(plot.title = element_text(size = 20, face = "bold"))

dev.off()