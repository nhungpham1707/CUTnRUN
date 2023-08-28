# This script is for steps in module 3,5 and 6. Adapt the directories wherever needed before running. 
# The script is a continutiation after peaks and peakcounts are extracted from module-3.sh. 
# Nhung 14 06 2023

library(DiffBind)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(ChIPseeker)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(viridis)
library(ggpubr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(uuutils)
library(rtracklayer)
library(VennDiagram)
library(DESeq2)
############################ Set directory to store output#########################

main_dir <- "~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/Rerun_14062023" 
DE_dir <- 'DE_output'
Annotation_dir <- 'Annotation_output'
Occupancy_dir <- 'Occupancy_output'
compare_dir <- 'Compare_output'
motif_dir <-'Motif'
dir.create(file.path(main_dir, DE_dir))
dir.create(file.path(main_dir, Annotation_dir))
dir.create(file.path(main_dir, Occupancy_dir))
dir.create(file.path(main_dir, compare_dir))

######################### load data and function ########################

setwd(main_dir)
load("~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/Rerun_14062023/2023-05-21-remove_low_depth_histone_samples-peak_count_DiffBind.RData")
source("~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/Rdata_Rscript/Identify_DE_peaks_function.R")

# universal set for hypergeometric 
# download gtf from here 'https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz'

# gtf_file <- import("~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/Rerun_14062023/Homo_sapiens.GRCh38.109.gtf")
# universal_set <- unique(gtf_file$gene_id)
# hypergeo(length(sc_in_data$seqnames), length(sc_symbol2ensemble$ENSEMBL), length(peak_annotation$geneId), length(universal_set), test.depletion = FALSE) # 7.05987e-07
# write.csv(universal_set, '~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/Rerun_14062023/ensembl_gene_set_hg38.csv')
universal_set <- read.csv('~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/Rerun_14062023/ensembl_gene_set_hg38.csv')
universal_set <- universal_set[,2] ; head(universal_set)
################################### Get differential binding sites #################################
setwd(file.path(main_dir,DE_dir))

contrast_number <- 3
save_name <- 'rerun_wo_lowdepth_histone'
fold_change_thredshold <- 2

test_peak_counts <- replace_zero_value(peaks, peak_counts, value_to_replace = 1)

message(paste('the number of 0 value: ', 
              length(which(test_peak_counts$peaks[[1]]$Reads == 0))/ length(test_peak_counts$peaks[[1]]$Reads))) 

message('Building contrast ...')        
test_contrast <- dba.contrast(test_peak_counts, design = '~Factor')

message('Analyzing contrast ....')
peak_model <- dba.analyze(test_contrast,bParallel=FALSE, bGreylist = FALSE)

message('saving volcano plot ...')
reso <- 600
png(filename = paste(Sys.time(), "_", save_name, "_volcano_plot.png") , width = 1200 * reso/72, 
    height = 700 * reso/72, units ="px", res = reso)
dba.plotVolcano(peak_model, contrast = contrast_number, dotSize = 2, bFlip = TRUE)  # bFlip = TRUE
dev.off()

png(filename = paste(Sys.time(), "_", save_name, "Fold2_volcano_plot.png") , width = 1200 * reso/72, 
    height = 700 * reso/72, units ="px", res = reso)
dba.plotVolcano(peak_model, contrast = contrast_number, fold = 1, dotSize = 2, bFlip = TRUE)  # bFlip = TRUE
dev.off()

res_deseq_contrast <- dba.report(peak_model, method=DBA_DESEQ2, contrast = 3, th=1)
res_deseq_contrast_df <- as.data.frame(res_deseq_contrast)
bed_contrast <- res_deseq_contrast_df %>% 
  dplyr::filter(FDR < 0.05 & Fold <  -fold_change_thredshold ) %>% 
  dplyr::select(seqnames, start, end, 'Conc_luciferase', 'Conc_fusion', Fold, p.value, FDR) # change the column names to yours
message(paste ('the number of DE sites is:', length(bed_contrast$seqnames) ))
bed_contrast <- bed_contrast[order(-abs(bed_contrast$Fold)),]
plot(bed_contrast$Fold)
bed_contrast$combineChrStart <- paste0(bed_contrast$seqnames, bed_contrast$start); head(bed_contrast)
write.table(bed_contrast[,c(1:3,9,4:8)], file=paste0(Sys.Date() ,save_name,
                                                     '_fold_change', fold_change_thredshold, 
                                                     '_FDR',0.05 ,'.bed'),
            sep="\t", quote=F, row.names=F, col.names=F) # this bed file is input for the next steps

luc_sites <- res_deseq_contrast_df %>% 
  dplyr::filter(FDR < 0.05) %>% 
  dplyr::select(seqnames, start, end, 'Conc_luciferase', 'Conc_fusion', Fold, p.value, FDR) # change the column names to yours
luc_sites[luc_sites$Fold > 0,]


fusion_specific_sites <- read.csv('~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/Rerun_14062023/DE_output/2023-06-15rerun_wo_lowdepth_histone_fold_change2_FDR0.05.bed', sep = '\t', header = F); head(fusion_specific_sites)
plot(fusion_specific_sites$V7)
write.table(fusion_specific_sites[1:1000,], 'top1000_fusion_specific_sites.bed', sep="\t", quote=F, row.names=F, col.names=F) # this file is for motif finding

######################################### Annotate peak #####################

setwd(file.path(main_dir,Annotation_dir))
samplefiles <- c(file.path(main_dir, DE_dir,"2023-06-15rerun_wo_lowdepth_histone_fold_change2_FDR0.05.bed"))
name <- 'lucvsfusion'

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
edb <- EnsDb.Hsapiens.v86

samplefiles <- as.list(samplefiles)
names(samplefiles) <- name
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE) # this step takes a while


# plot the distribution

plotAnnoBar(peakAnnoList[[name]])

#### Gene annotation 
PeakAnno_df <- data.frame(peakAnnoList[[name]]@anno)
PeakAnno_entrez <- PeakAnno_df$geneId
# Return the gene symbol for the set of Entrez IDs
PeakAnno_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                      keys = PeakAnno_entrez,
                                      columns = c("GENENAME"),
                                      keytype = "ENTREZID"); head(PeakAnno_edb)

go_enrich <- enrichGO(gene = PeakAnno_edb$ENTREZID, 
                      keyType = "ENTREZID", 
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = TRUE)

go_enrich # 51 enrich terms

plot = dotplot(go_enrich, showCategory=15, font.size = 30) 
annotate_figure(plot, top = text_grob("GO enrichment fusion vs luciferase", 
                                      color = "black", face = "bold", size = 30))
go_enrich_df <- data.frame(go_enrich); head(go_enrich_df)
length(unique(go_enrich_df$ID))
write.csv(go_enrich_df, paste0(Sys.time(),"GO_enrichment.csv"))

## ensembl annotation
## get the txdb and peak files
peak <- ChIPseeker::readPeakFile(samplefiles[[1]])

## change the seqlevels of peak files
oldseqlevels <- data.frame(old=GenomeInfoDb::seqlevels(peak),
                           new=gsub("chr","",GenomeInfoDb::seqlevels(peak)))
newseqlevels <- data.frame(new=c(1:22,"X","Y"))
combine <- merge(oldseqlevels,newseqlevels,by="new")
newnames <- combine$new
names(newnames) <- combine$old
peak <- GenomeInfoDb::renameSeqlevels(peak,newnames)
peak
annotatePeak_result <- ChIPseeker::annotatePeak(peak = peak,TxDb = edb) # take time 
annotEDB_df <- data.frame(annotatePeak_result@anno); head(annotEDB_df)

go_enrich_ensem <- enrichGO(gene = annotEDB_df$geneId,
                            keyType = "ENSEMBL", 
                            OrgDb = org.Hs.eg.db, 
                            ont = "BP", 
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05, 
                            readable = TRUE)

dotplot(go_enrich_ensem) # --> different result than the other within 1000TSS
plotAnnoBar(annotatePeak_result)
# colname for DE peaks
colnames(annotEDB_df)[6:11] = c('combineChrStart', 'Conc_luciferase','Conc_fusion', 'Fold', 'p.value', 'FDR')
#colname for common peaks
# colnames(annotEDB_df)[6:9] = c('combineChrStart', 'score', 'RPKM', 'reads')

head(annotEDB_df)
# Get gene symbols
dataset <- annotEDB_df
ensembl <- dataset$geneId

# Return the gene symbol for the set of ensemble ids
annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                         keys = ensembl,
                                         columns = c("GENENAME"),
                                         keytype = "GENEID")
annotations_edb$GENEID <- as.character(annotations_edb$GENEID)


annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                         keys = annotEDB_df$geneId,
                                         columns = c("GENENAME"),
                                         keytype = "GENEID")


# Write to file
peak_annotation <- dataset %>% 
  left_join(annotations_edb, by=c("geneId"="GENEID"))
head(peak_annotation)
length(peak_annotation$seqnames) #2653
length(unique(peak_annotation$GENENAME))
write.table(peak_annotation, file=paste0(Sys.Date(),"functional_analysis_in_", name, ".txt"), sep="\t", quote=F, row.names=F)

save.image(file.path(main_dir,paste0(name,"_annotation.RData")))

############################# Occupancy analysis ###################################
setwd(file.path(main_dir, Occupancy_dir))

tfe3_peak <- peak_model[[1]][1:4] # peak_model is output from diffbind contrast analyze
fusion_peak <- peak_model[[1]][5:8]
luc_peak <- peak_model[[1]][9:12]

thredshold <- 10
tfe3_all_rep <- write_top_peak(tfe3_peak, thredshold, 'tfe3')
tfe3_all_rep
fusion_all_rep <- write_top_peak(fusion_peak,thredshold , 'fusion')
luc_all_rep <- write_top_peak(luc_peak, thredshold, 'luc')

common_peak <- intersect(fusion_all_rep$combineChrStart, luc_all_rep$combineChrStart)
length(common_peak)  # 3291 in all except low depth histone
venn.diagram(x = list(fusion_all_rep$combineChrStart, luc_all_rep$combineChrStart),
             filename = paste(Sys.time(),'s3norm_all_except_low_depth_histone_fusion_luc_rpmk', thredshold,'.png'), 
             category.names = c("fusion", "luc") )

# remove insignificant luc sites 
luc_sites <- setdiff(luc_all_rep$combineChrStart, fusion_all_rep$combineChrStart)
venn.diagram(x = list(fusion_all_rep$combineChrStart, common_peak),
             filename = paste(Sys.time(),'s3norm_all_except_low_depth_histone_fusion_luc_rpmk', thredshold,'_remove_insig_luc_sites.png'), 
             category.names = c("fusion", "luc"), alpha = rep(0.5,2), fill = c('light blue', 'pink') ) 

# remove insignificant fusion sites
sig_fusion_sites <- paste0('chr',bed_contrast$seqnames, bed_contrast$start)
fusion_sites_combine <- c(sig_fusion_sites, common_peak)
venn.diagram(x = list(fusion_sites_combine, common_peak),
             filename = paste(Sys.time(),'s3norm_all_except_low_depth_histone_fusion_luc_rpmk', 
                              thredshold,'_remove_insig_fus_sites.png'), 
             category.names = c("fusion", "luc"), alpha = rep(0.5,2), fill = c('light blue', 'pink') ) 

# check if these top 2000 common peaks have FDR < 0.05 in DE site
all_DE <- res_deseq_contrast_df %>% 
  dplyr::filter(FDR < 0.05) %>% 
  dplyr::select(seqnames, start, end, 'Conc_luciferase', 'Conc_fusion', Fold, p.value, FDR)
all_DE$combineChrStart <- paste0(all_DE$seqnames, all_DE$start)
length(intersect(all_DE$combineChrStart, common_peak_res_df$combineChrStart[1:2000])) # 683 sites are there 

# remove sites with FDR < 0.05 from common peak 
common_peak_without_sig_FDR <- setdiff(common_peak, all_DE$combineChrStart); head(common_peak_without_sig_FDR); length(common_peak_without_sig_FDR)

# write to file for motif analysis 

# extract top 200 sites from common peak for motif analysis 

common_peak_in_luc <- luc_all_rep[luc_all_rep$combineChrStart %in% common_peak_without_sig_FDR,]
length(common_peak_in_luc$combineChrStart)

common_peak_in_fusion <- fusion_all_rep[fusion_all_rep$combineChrStart %in% common_peak_without_sig_FDR,]
length(common_peak_in_fusion$combineChrStart)

data.frame(common_peak_in_fusion$combineChrStart, common_peak_in_luc$combineChrStart)
plot(common_peak_in_fusion$RPKM, common_peak_in_luc$RPKM)

common_peak_df <- common_peak_in_luc[order(-common_peak_in_luc$RPKM),]
plot(common_peak_df$RPKM)

common_peak_res_df <- data.frame(Chr = common_peak_df$Chr,
                                 Start = common_peak_df$Start,
                                 End = common_peak_df$End,
                                 combineChrStart = common_peak_df$combineChrStart,
                                 Score = common_peak_df$Score,
                                 RPKM = common_peak_df$RPKM,
                                 Reads = common_peak_df$Reads); head(common_peak_res_df)

write.table(common_peak_res_df, file=paste0(Sys.Date(),'common_peaks_wo_sig_peaks_fusion_luc_s3norm_all_except_low_dephistone.bed'),
            sep="\t", quote=F, row.names=F, col.names=F)

################################ Compare with Others ##################################
peak_annotation <- read.csv("~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/Rerun_14062023/Annotation_output/2023-06-15functional_analysis_in_DE.txt", sep = '\t'); head(peak_annotation)
setwd(file.path(main_dir, compare_dir))
## With Sc RNA seq ----
sc_rna <- read.csv('~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/FROM_PHILIP/e16_DE_genes-byCMO.tsv', sep = '\t') # file with genes, log2fc, padj
head(sc_rna)

# by name
same_sc <- intersect(sc_rna$gene, peak_annotation$GENENAME)
length(same_sc) # 20 gene similar 

# by ensembl
sc_symbol2ensemble <- clusterProfiler::bitr(sc_rna$gene, fromType = "SYMBOL",toType = c("ENTREZID","SYMBOL","ENSEMBL"), OrgDb = org.Hs.eg.db)

same_peak_sc <- intersect(sc_symbol2ensemble$ENSEMBL, peak_annotation$geneId); head(same_peak_sc)
length(same_peak_sc)

sc_in_data <- peak_annotation[peak_annotation$geneId %in% same_peak_sc,]; head(sc_in_data); length(sc_in_data$seqnames)
for (i in 1:length(sc_in_data$GENENAME)){
  index_in_sc <- which(sc_symbol2ensemble$ENSEMBL == sc_in_data$geneId[i])
  sc_in_data$sc_avg_log2FC [i] <- sc_rna[sc_rna$gene == sc_symbol2ensemble$SYMBOL[index_in_sc],3]
  sc_in_data$sc_cluster[i] <- sc_rna[sc_rna$gene == sc_symbol2ensemble$SYMBOL[index_in_sc],7]
}
head(sc_in_data)
sc_in_data[,c(1:3,11,19,21:23)]
length(sc_in_data$seqnames) #29 genes 
write.csv(sc_in_data, "shared_DE_Sites_with_SC_RNAseq.csv")

unique_gene_type <- data.frame(type = unique(sc_in_data$transcriptBiotype))
for (i in 1 : length(unique_gene_type$type)){
  unique_gene_type$count [i] <- length(sc_in_data[sc_in_data$transcriptBiotype == unique_gene_type$type[i],1])
}
unique_gene_type
reso <- 600

png(filename = paste0(Sys.time(), "-compare_DE_w_scRNA_genetypes.png") , width = 1400 * reso/72, height = 700 * reso/72, units ="px", res = reso)
unique_gene_type %>% ggplot(aes(x = type, y = count, fill = type)) +
  geom_bar(stat="identity") + theme_bw(base_size = 24) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                               panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  xlab("") + ylab("") + ggtitle( "Transcript types of overlapping genes between cutnrun and ScRNA)") +
  theme(plot.title = element_text(size = 20, face = "bold"))

dev.off()


unique_cluster <- data.frame(type=unique(sc_in_data$sc_cluster))
for (i in 1 : length(unique_cluster$type)){
  unique_cluster$count [i] <- length(sc_in_data[sc_in_data$sc_cluster == unique_cluster$type[i],1])
}
unique_cluster

png(filename = paste0(Sys.time(), "-compare_DE_w_scRNA_cluster.png") , width = 1400 * reso/72, height = 700 * reso/72, units ="px", res = reso)
unique_cluster %>% ggplot(aes(x = type, y = count, fill = type)) +
  geom_bar(stat="identity") + theme_bw(base_size = 24) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                               panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  xlab("") + ylab("") + ggtitle( "Cluster of overlapping genes between cutnrun and ScRNA)") +
  theme(plot.title = element_text(size = 20, face = "bold"))

dev.off()

########## compare with literature ----

zb_etal <- read.csv("~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/FROM_LITERATURE/ZB_etal_mmc5.csv", sep = ";") 
head(zb_etal)
length(zb_etal$Genes) # 139 genes

# compare by name
same_lit <- intersect(zb_etal$Genes, peak_annotation$GENENAME); length(same_lit) 
in_lit <- zb_etal[zb_etal$Genes%in%same_lit,]
in_data <- peak_annotation[peak_annotation$GENENAME %in% same_lit,]
for (i in 1:length(in_data$GENENAME)){
  in_data$ZB_etal [i] <- in_lit[in_lit$Genes == in_data$GENENAME[i],2]
}
in_data[,c(1:3,11,19,21,22)]
in_lit
write.csv(in_data, "compare_DE_Sites_with_ZB_et_al.csv")

# compare by ensembl
zb_symbol2ensemble <- clusterProfiler::bitr(zb_etal$Genes, fromType = "SYMBOL",toType = c("ENTREZID","SYMBOL","ENSEMBL"), OrgDb = org.Hs.eg.db)
same_ensem <- intersect(zb_symbol2ensemble$ENSEMBL, peak_annotation$geneId)
length(same_ensem) # 14
gene_lit <- zb_symbol2ensemble[zb_symbol2ensemble$ENSEMBL %in% same_ensem,]
in_lit_ensem <- zb_etal[zb_etal$Genes %in% gene_lit$SYMBOL,]
in_lit_ensem # similar genes with by name

# hypergeometric
hypergeo(length(in_data$seqnames), length(zb_etal$Genes), length(peak_annotation$seqnames), length(universal_set), test.depletion = TRUE)

## compare with bulk RNAseq ----
load("~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/FROM_BULK_RNAseq/CR_RNA_analysis_Aug2022.RData")
# 
de_genes_df <- data.frame(results(dds3, contrast = c("condition", "Fusion", "Luciferase"), 
                                  pAdjustMethod = "BH",cooksCutoff=FALSE)) %>%
  dplyr::filter(padj < 0.05)
#
length(de_genes_df$baseMean)
deseq_res <- results(dds3, contrast = c("condition", "Fusion", "Luciferase"), 
                     pAdjustMethod = "BH",cooksCutoff=FALSE)
saveRDS(deseq_res, '~/Dropbox/POSTDOC_MAXIMA/PROJECTS/tRCC_CnR/data/processed/bulk_RNA_result.RDS')

de_genes_from_file <- read.csv("~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/FROM_BULK_RNAseq/DEGs_FusionVLuc_308G_BH0.05.csv", sep = ';', header = FALSE); head(de_genes_from_file)
colnames(de_genes_from_file) <- c('GeneName', 'baseMean', 'log2FC', 'lfcSE', 'stat', 'pvalue', 'padj')
de_genes_from_file <- de_genes_from_file[2:nrow(de_genes_from_file),c(1:7)]; head(de_genes_from_file)
length(de_genes_from_file$GeneName)

de_genes_from_file$log2FC <- sub(",",".",de_genes_from_file$log2FC)
de_genes_from_file$baseMean <- sub(",",".",de_genes_from_file$baseMean)
de_genes_from_file$lfcSE <- sub(",",".",de_genes_from_file$lfcSE)
de_genes_from_file$stat <- sub(",",".",de_genes_from_file$stat)
de_genes_from_file$pvalue <- sub(",",".",de_genes_from_file$pvalue)
de_genes_from_file$padj <- sub(",",".",de_genes_from_file$padj); head(de_genes_from_file)

# same_bulk <- intersect(bulk_rna, peak_annotation$GENENAME)
same_bulk <- intersect(de_genes_from_file$GeneName, peak_annotation$GENENAME)
length(same_bulk) #168 genes # 61 genes 
which(same_bulk == 'GPNMB')
same_peak_bulk <- peak_annotation[peak_annotation$GENENAME %in% same_bulk,]
which(same_peak_bulk$GENENAME == 'GPNMB')

which(peak_annotation$GENENAME == 'GPNMB')
which(de_genes_from_file$GeneName == 'GPNMB')
hypergeo(length(same_peak_bulk$seqnames), length(bulk_rna), length(peak_annotation$seqnames), length(universal_set), test.depletion = F)

same_peak_bulk
same_peak_bulk_in_peak <- peak_annotation[peak_annotation$GENENAME %in% same_peak_bulk$GENENAME,]
same_peak_bulk_in_peak <- same_peak_bulk_in_peak[order(same_peak_bulk_in_peak$Fold),]
which(same_peak_bulk_in_peak$GENENAME == 'GPNMB')

same_peak_bulk_in_bulk <- de_genes_from_file[de_genes_from_file$GeneName %in% same_peak_bulk$GENENAME,]
head(same_peak_bulk_in_bulk)
merge_cr_rna <- merge(same_peak_bulk_in_peak, same_peak_bulk_in_bulk, by.x = 'GENENAME', by.y = 'GeneName')
head(merge_cr_rna)

# write.csv(merge_cr_rna, 'shared_DE_Sites_w_bulkrna_from_xls_file.csv')

same_bulk_sc <- intersect(same_bulk, same_sc)
same_cutnrun_bulk_sc <- peak_annotation[peak_annotation$GENENAME %in% same_bulk_sc,c(1:3,11,19,21)]
write.csv(same_cutnrun_bulk_sc, "overlap_DE_in_cutnrun_scRNA_bulkRNA.csv")
# compare with rnaseq in ZB paper
zb_rna <- c("SNCB", 'TRIM67', 'IRX6', 'SQSTM1', 'TMEM64', 'SLC39A1', 'RAB7A', 'RHEB', 'RRAGC','ATP6V1C1') # SQSTM1 is in 
intersect(zb_rna, zb_etal$Genes)
intersect(peak_annotation$GENENAME, zb_rna)
peak_annotation[peak_annotation$GENENAME=='SQSTM1',]

# identify genes in C&R and bulk RNAseq with motif ----
de_mitf <- read.csv("~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/Rerun_14062023/Motif/14062013_DE_sites_len10/MITF_motif/peak_with_MITF_motif.csv"); head(de_mitf)
cr_rna_genes_w_mitf_motif <- intersect(de_mitf$GENENAME, merge_cr_rna$GENENAME) #12 genes 
cr_rna_mitf_motif <- merge_cr_rna[merge_cr_rna$GENENAME %in% cr_rna_genes_w_mitf_motif,]
length(unique(cr_rna_mitf_motif$GENENAME))
write.csv(cr_rna_mitf_motif, 'overlap_DE_in_cutnrun_bulk_RNA_w_MITF_motif.csv')

de_tfe3 <- read.csv("~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/Rerun_14062023/Motif/14062013_DE_sites_len81012/TFE3_in_DE_sites/all_peaks_w_motif_TFE3_in_DE_sites.csv"); head(de_tfe3)
cr_rna_genes_w_tfe3_motif <- intersect(de_tfe3$GENENAME, merge_cr_rna$GENENAME); head(cr_rna_genes_w_tfe3_motif)
length(cr_rna_genes_w_tfe3_motif) #12 genes 
cr_rna_tfe3_motif <- merge_cr_rna[merge_cr_rna$GENENAME %in% cr_rna_genes_w_tfe3_motif,]
write.csv(cr_rna_tfe3_motif, 'overlap_DE_in_cutnrun_bulk_RNA_w_tfe3_motif.csv')
which(de_tfe3$GENENAME == 'GPNMB') # not in DE site probably because they are in the common site
length(de_tfe3$combineChrStart) # 682 sites 


common_tfe3 <- read.csv("~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/Rerun_14062023/Motif/14062013_common_sites_DEfold2_rpkm10_len10/all_peaks_w_motif_ TFE3_in_common_sites .csv"); head(common_tfe3)
common_tfe3[which(common_tfe3$GENENAME == 'GPNMB'),] # yes because they are in the common peaks

# inspect tfe3 sites in de sites ----
de_tfe3$p.value[de_tfe3$GENENAME %in% cr_rna_genes_w_tfe3_motif]
de_tfe3[order(de_tfe3$p.value),]

