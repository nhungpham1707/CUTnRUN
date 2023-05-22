# this script contains function to contrast and identify differential binding sites from peak data obtained from diffBind

# Nhung 22 May 2023

# if data contain 0 values, Deseq will not be able to run differential analysis. 1st option is to add
# a pseudo 1 to replace 0 in all the data 
# Input: 
# - peaks : output from  peaks <- dba(sampleSheet=samples)
# - peak_counts : output from peak_counts <- dba.count(peaks)
# - value_to_replace : value to replace for 0, by default 1
replace_zero_value <- function(peaks, peak_counts, value_to_replace=1) {
  for (i in 1: length(peaks$samples[[1]])) {
      for (x in 1: length(peak_counts$peaks[[i]]$Reads)){
      if (peak_counts$peaks[[i]]$Reads[x] == 0) { 
        peak_counts$peaks[[i]]$Reads[x] = value_to_replace
      } 
    }
  }
  message(paste("replace_value is", value_to_replace))
  return (peak_counts)
  }


## write bed file for DE sites
# Input: 
#   - peak_model output from peak_model <- dba.analyze(test_contrast,bParallel=FALSE, bGreylist = FALSE)
#   - contrast_number: specify which contrast to identify DE. Can be checked from peak_model
#   - fold_change_thredshold: above this value is considered significantly differen. Default: 2
#   - FDR_thredshold: thredshold for false discovery rate
#   - column_to_keep: beside seqnames, start, end, fold, pvalue and FDR. check dba.report(peak_model)
# to identify colum names
    # - save_name: name to save file
# Output: bed file save in the current dir


write_DE_bed <- function(peak_model, contrast_number = 3, fold_change_thredshold = 2, 
                         FDR_thredshold = 0.05, 
                         column_to_keep = c('Conc_luciferase', 'Conc_fusion'), 
                         save_name ) {
res_deseq_contrast <- dba.report(peak_model, method=DBA_DESEQ2, contrast = contrast_number, th=1)
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
}

## conduct desqe2 analysis and save plot and bed file
## conduct the diffBind analysis
# Input:
# - peaks : output from  peaks <- dba(sampleSheet=samples)
# - peak_counts : output from peak_counts <- dba.count(peaks)
# - value_to_replace : value to replace for 0, by default 1
#   - peak_model output from peak_model <- dba.analyze(test_contrast,bParallel=FALSE, bGreylist = FALSE)
#   - contrast_number: specify which contrast to identify DE. Can be checked from peak_model
#   - fold_change_thredshold: above this value is considered significantly differen. Default: 2
#   - FDR_thredshold: thredshold for false discovery rate
#   - column_to_keep: beside seqnames, start, end, fold, pvalue and FDR. check dba.report(peak_model)
# to identify colum names
# - save_name: name to save file
# Output: 
# - bed file save in the current dir
# - volcano plot in png
# - peak_model 

conduct_diffBind_analysis <- function(peaks, peak_counts, save_name, 
                                      column_to_keep = c('Conc_luciferase', 'Conc_fusion'),
                                      reso = 600, contrast_number = 3, 
                                      fold_change_thredshold = 2, 
                                      FDR_thredshold = 0.05, 
                                      value_to_replace = 1){
  
message('Replacing 0 value before conducting Deseq ...')
test_peak_counts <- replace_zero_value(peaks, peak_counts, value_to_replace)

message(paste('the number of 0 value: ', 
              length(which(test_peak_counts$peaks[[1]]$Reads == 0))/ length(test_peak_counts$peaks[[1]]$Reads))) 
        
message('Building contrast ...')        
test_contrast <- dba.contrast(test_peak_counts, design = '~Factor')

message('Analyzing contrast ....')
peak_model <- dba.analyze(test_contrast,bParallel=FALSE, bGreylist = FALSE)
return (peak_model)

message('saving volcano plot ...')
png(filename = paste(res_dir, '/', Sys.time(), "_", save_name, "_volcano_plot.png") , width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
dba.plotVolcano(peak_model, contrast = contrast_number, dotSize = 2)  # bFlip = TRUE
dev.off()

message('write DE to bed file ...')
write_DE_bed(peak_model, save_name, contrast_number, fold_change_thredshold, FDR_thredshold,column_to_keep)

message('finished!')
}
