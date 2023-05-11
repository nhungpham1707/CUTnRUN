# Annotate differential binding sites resulted from DiffBind_analysis.R.
# Output: peak annotation with plot on elelment distribution (promoters, intron, exon, ect), enriched GO terms and pathways 
# Ref. https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/12_annotation_functional_analysis.html
# Nhung 09 05 2023
# library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(ChIPseeker)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
##Functional analysis
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
edb <- EnsDb.Hsapiens.v86
res_dir <- Sys.getenv("DIFFBIND_RESULT_DIR_VARIABLE")

setwd("~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN")

#check that we are working with the same sequence levels
# samplefiles <- c(
#   "fusion_overlap.bed", "tfe3_overlap.bed",
#   "luciferase_overlap.bed",
#   "fusionLuc_batchCorrectFC2.bed",
#   "fusionTfe3_batchCorrectFC2.bed")
samplefiles <- c("test_contrast1.bed", "test_contrast3.bed")
samplefiles <- c("test_contrast1_remove_control.bed", "test_contrast3_remove_control.bed")
samplefiles <- c("test_contrast1_remove_control_fold1.bed", "test_contrast3_remove_control_fold1.bed")

setwd("~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/FROM_HPC")
setwd('~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/FROM_HPC/test_diffBind_HPC')
load("2023-10-05-diffBind-s3norm-peak-annotation.RData")
samplefiles <- c("2023-05-10-diffBind_contrast1_s3norm.bed", "2023-05-10-diffBind_contrast3_s3norm.bed")
samplefiles <- c("2023-05-10-diffBind_contrast1_s3norm_fold1.bed", "2023-05-10-diffBind_contrast3_s3norm_fold1.bed")
samplefiles <- c("2023-05-10-diffBind_contrast1_s3norm_fold1_without_control.bed", "2023-05-10-diffBind_contrast3_s3norm_fold1_without_control.bed")
samplefiles <- c("2023-05-10-diffBind_contrast1_s3norm_fold1_without_control_fix_fdr.bed", "2023-05-10-diffBind_contrast3_s3norm_fold1_without_control_fix_FDR.bed")
samplefiles <- as.list(samplefiles)
# names(samplefiles) <- c("Fusion_vs_tfe3", "Fusion_vs_luc")
names(samplefiles) <- c("tfe3_vs_fusion", "luc_vs_fusion")
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE) # this step takes a while
reso = 600
png(filename = "DE_bind_feature_distribution_s3norm.png" , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)

plotAnnoBar(peakAnnoList[1:2])
dev.off()

png(filename = "Distribution_of_TF_s3norm.png" , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
plotDistToTSS(peakAnnoList[1:2], title="Distribution of transcription factor-binding loci \n relative to TSS")
dev.off()

peakFusion <- readPeakFile(samplefiles[[1]])
peakTFE3<- readPeakFile(samplefiles[[2]])
peakLuc <- readPeakFile(samplefiles[[3]])
covplot(peakFusion, weightCol="V5", chrs=c(paste0("chr", c(1:22, "X"))))
covplot(peakTFE3, weightCol="V5", chrs=c(paste0("chr", c(1:22, "X"))))

##Profile of ChIP peaks binding to TSS -- this step takes time and memory
promoter <- getPromoters(TxDb=txdb, upstream=100000, downstream=100000)
tagMatrixList <- lapply(samplefiles, getTagMatrix, windows=promoter) # this step is time-consuming and memory-consuming ~30 mins 
plotAvgProf(tagMatrixList, xlim=c(-100000, 100000), conf=0.95,resample=500, facet="row")
plotPeakProf2(samplefiles, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body",
              TxDb = txdb, facet = "row", nbin = 800)

##Functional comparison pathways
fusionLuc_annot <- data.frame(peakAnnoList[["luc_vs_fusion"]]@anno)
fusionTfe3_annot <- data.frame(peakAnnoList[["tfe3_vs_fusion"]]@anno)

fusiontfe3_annot_entrez <- fusionTfe3_annot$geneId
fusionLuc_annot_entrez <- fusionLuc_annot$geneId
# Return the gene symbol for the set of Entrez IDs
fusiontfe3_annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                         keys = fusiontfe3_annot_entrez,
                                         columns = c("GENENAME"),
                                         keytype = "ENTREZID")
fusiontfe3_annotations_edb

ego <- enrichGO(gene = fusiontfe3_annot_entrez, 
                keyType = "ENTREZID", 
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
ego
cluster_summary <- data.frame(ego)
cluster_summary
 write.csv(cluster_summary, paste0(Sys.Date(),"clusterProfiler_fusiontfe3_fold2_s3norm_w_control.csv"))
dotplot(ego, showCategory=30)

# fusion vs luc
fusionluc_annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                                    keys = fusionLuc_annot_entrez,
                                                    columns = c("GENENAME"),
                                                    keytype = "ENTREZID")
fusionluc_annotations_edb

ego_fuLuc <- enrichGO(gene = fusionluc_annotations_edb$ENTREZID, 
                keyType = "ENTREZID", 
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
ego_fuLuc
dotplot(ego_fuLuc, showCategory=30)
cluster_summary_luc <- data.frame(ego_fuLuc)
write.csv(cluster_summary_luc, paste0(Sys.Date(),"clusterProfiler_fusionLuc_fold2_s3norm_w_control.csv"))

# genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
# names(genes) = sub("_", "\n", names(genes))
# compKEGG <- compareCluster(geneCluster   = genes,
# fun           = "enrichKEGG",
# pvalueCutoff  = 0.05,
# pAdjustMethod = "BH")
# dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")
genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)


##Using ensembl annotation
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
tfe3Fusion_annotEDB <- data.frame(annotatePeak_result@anno)
tfe3Fusion_annotEDB
# Get the entrez IDs
dataset <- tfe3Fusion_annotEDB
ensembl <- dataset$geneId

# Return the gene symbol for the set of Entrez IDs
annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                         keys = ensembl,
                                         columns = c("GENENAME"),
                                         keytype = "GENEID")
annotations_edb$GENEID <- as.character(annotations_edb$GENEID)
# Write to file
dataset %>% 
  left_join(annotations_edb, by=c("geneId"="GENEID")) %>% 
  write.table(file="functional_analysis_tfe3vsfusion_EDB.txt", sep="\t", quote=F, row.names=F)


# ## get the txdb and peak files
peak_luc_fu <- ChIPseeker::readPeakFile(files[[2]])

## change the seqlevels of peak files
oldseqlevels <- data.frame(old=GenomeInfoDb::seqlevels(peak_luc_fu),
                           new=gsub("chr","",GenomeInfoDb::seqlevels(peak_luc_fu)))
newseqlevels <- data.frame(new=c(1:22,"X","Y"))
combine <- merge(oldseqlevels,newseqlevels,by="new")
newnames <- combine$new
names(newnames) <- combine$old
peak_luc_fu <- GenomeInfoDb::renameSeqlevels(peak_luc_fu,newnames)
annotatePeak_result <- ChIPseeker::annotatePeak(peak = peak_luc_fu,TxDb = edb) # take time 
lucFusion_annotEDB <- data.frame(annotatePeak_result@anno)
lucFusion_annotEDB
# Get the entrez IDs
dataset <- lucFusion_annotEDB
ensembl <- dataset$geneId

# Return the gene symbol for the set of Entrez IDs
annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                         keys = ensembl,
                                         columns = c("GENENAME"),
                                         keytype = "GENEID")
annotations_edb$GENEID <- as.character(annotations_edb$GENEID)
# Write to file
dataset %>% 
  left_join(annotations_edb, by=c("geneId"="GENEID")) %>% 
  write.table(file="functional_analysis_lucvsfusion_EDB.txt", sep="\t", quote=F, row.names=F)
# ####### compare my and terezinha enrich pathways

# from terezinha
to_compare <- read.csv("~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/Terezinha_clusterProfiler_fusiontfe3.csv")
# from my analysis
fusiontfe <- read.csv("clusterProfiler_fusiontfe3_fold1.csv")

similar <- intersect(to_compare$ID, fusiontfe$ID)
length(similar) # 590 pathways

in_my_analysis <- setdiff(to_compare$ID, fusiontfe$ID)
length(in_my_analysis) # 235 pathways

in_other <- setdiff(fusiontfe$ID, to_compare$ID)
length(in_other) # 319 pathways

# compare top pathways 

head(fusiontfe)
sort_fusiontfe <- fusiontfe
sort_fusiontfe <- sort_fusiontfe[order(sort_fusiontfe$p.adjust),]

sort_to_compare <- to_compare[order(to_compare$p.adjust),]

top_similar <- intersect(sort_fusiontfe[1:100,3], sort_to_compare[1:100,3])
length(top_similar) # 53 pathways the same in the top 100
sort_to_compare[1:10,1:8]
sort_fusiontfe[1:10,1:8]
