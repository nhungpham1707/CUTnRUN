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


write_DE_bed <- function(peak_model, res_dir, contrast_number = 3, fold_change_thredshold = 2, 
                         FDR_thredshold = 0.05, 
                         column_to_keep = c('Conc_luciferase', 'Conc_fusion'), 
                         save_name ) {
message('start write to bed file')
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

conduct_diffBind_analysis <- function(peaks, peak_counts, save_name, res_dir,  
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

message('write DE to bed file ...')
write_DE_bed(peak_model, res_dir, save_name, contrast_number, fold_change_thredshold, FDR_thredshold,column_to_keep)

message('saving volcano plot ...')
png(filename = paste(res_dir, '/', Sys.time(), "_", save_name, "_volcano_plot.png") , width = 1200 * reso/72, 
    height = 700 * reso/72, units ="px", res = reso)
dba.plotVolcano(peak_model, contrast = contrast_number, dotSize = 2)  # bFlip = TRUE
dev.off()

message('finished!')
}

# make function for finding overlap peaks between samples

# # function with bigger start in sample
# 
# overlay_peak_all_cases <- function(sample_peak, validate_peak){
#   de_chr_list <- unique(sample_peak$seqnames)
#   result <- c()
#   for (i in 1: length(de_chr_list)){
#     # extract peak on one chromosome in DE sites
#     sample_peak_chr <- sample_peak[sample_peak$seqnames == de_chr_list[i] ,]
#     head(sample_peak_chr)
#     length(sample_peak_chr$seqnames) # 224 chr1 158 chr10
#     sample_peak_chr_sort <- sample_peak_chr[order(sample_peak_chr$start),]
#     head(sample_peak_chr_sort)
#     print(paste('the sig peak on chr',de_chr_list[i],'start at', min(sample_peak_chr_sort$start)))
#     print(paste('the last sig peak on chr', de_chr_list[i], 'end with', max(sample_peak_chr_sort$end)))
#     # extract peak on one chromosome in fusion H3k27ac
#     validate_peak_chr <- validate_peak[validate_peak$chr == paste0('chr',de_chr_list[i]),]
#     length(validate_peak_chr$chr) # 127
#     head(validate_peak_chr)
#     print(paste('peak in h3k27ac on chr', de_chr_list[i], 'start at', min(validate_peak_chr$start)))
#     print(paste('the last sig peak on chr', de_chr_list[i], 'end with', max(validate_peak_chr$end)))
#     start_in_validate_peak <- min(validate_peak_chr$start)
#     
#     if (max(validate_peak_chr$end) < min(sample_peak_chr_sort$start)) {
#       result <- data.frame(chr = de_chr_list[i], 
#                            start_in_sample = 'peaks in validate end before peaks in sample start', 
#                            end_in_sample = 'peaks in validate end before peaks in sample start', 
#                            start_in_validate = 'peaks in validate end before peaks in sample start',
#                            end_in_validate = 'peaks in validate end before peaks in sample start') %>% rbind(result,.)} 
#     else if (max(sample_peak_chr_sort$end) < min(validate_peak_chr$start)) {
#       result <- data.frame(chr = de_chr_list[i], 
#                            start_in_sample = 'peaks in sample end before peaks in validate start', 
#                            end_in_sample = 'peaks in sample end before peaks in validate start', 
#                            start_in_validate = 'peaks in sample end before peaks in validate start',
#                            end_in_validate = 'peaks in sample end before peaks in validate start') %>% rbind(result,.)} 
#     else if ((min(sample_peak_chr_sort$start) > start_in_validate_peak) & (max(sample_peak_chr_sort$end) > min(validate_peak_chr$start))){
#       
#       result <- data.frame(chr = de_chr_list[i], 
#                            start_in_sample = 'start in sample is bigger than validate sample but end is bigger', 
#                            end_in_sample = 'start in sample is bigger than validate sample', 
#                            start_in_validate = 'start in sample is bigger than validate sample',
#                            end_in_validate = 'start in sample is bigger than validate sample') %>% rbind(result,.) }
#     
#     else if ((min(sample_peak_chr_sort$start) < min(validate_peak_chr$start)) & (min(sample_peak_chr_sort$start) < max(validate_peak_chr$end))) {
#       tail(sample_peak_chr_sort[sample_peak_chr_sort$start < start_in_validate_peak,], n=1)
#       find_validate_start_in_sample <- max(sample_peak_chr_sort$start[sample_peak_chr_sort$start < start_in_validate_peak])
#       validate_start_in_sample_index <- which(sample_peak_chr_sort$start == find_validate_start_in_sample)
#       min(abs(as.numeric(sample_peak_chr_sort$start[validate_start_in_sample_index]) - as.numeric(validate_peak_chr$start))) 
#       for (j in validate_start_in_sample_index : length(sample_peak_chr_sort$seqnames)) {
#         distance_sample_to_validate <- as.numeric(sample_peak_chr_sort$start[j]) - as.numeric(validate_peak_chr$start)
#         index_in_validate <- which(abs(distance_sample_to_validate) == min(abs(distance_sample_to_validate)))
#         
#         result <- data.frame(chr = de_chr_list[i],
#                              start_in_sample = sample_peak_chr_sort$start[j],
#                              end_in_sample = sample_peak_chr_sort$end[j],
#                              start_in_validate = validate_peak_chr$start[index_in_validate],
#                              end_in_validate = validate_peak_chr$end[index_in_validate]) %>% rbind(result,.)
#       }
#     }
#   }
#   return (result) 
# }

## new approach only 2 cases with the start sites 
overlay_peak_start_site <- function(sample_peak, validate_peak){
  de_chr_list <- unique(sample_peak$seqnames)
  result <- c()
  count_big_start_sample <- 0
  count_small_start_sample <- 0
  for (i in 1: length(de_chr_list)){
    # extract peak on one chromosome in DE sites
    sample_peak_chr <- sample_peak[sample_peak$seqnames == de_chr_list[i] ,]
    head(sample_peak_chr)
    length(sample_peak_chr$seqnames) # 224 chr1 158 chr10
    sample_peak_chr_sort <- sample_peak_chr[order(sample_peak_chr$start),]
    head(sample_peak_chr_sort)
    print(paste('the sig peak on chr',de_chr_list[i],'start at', min(sample_peak_chr_sort$start)))
    
    # extract peak on one chromosome in fusion H3k27ac
    validate_peak_chr <- validate_peak[validate_peak$chr == paste0('chr',de_chr_list[i]),]
    length(validate_peak_chr$chr) # 127
    head(validate_peak_chr)
    print(paste('peak in h3k27ac on chr', de_chr_list[i], 'start at', min(validate_peak_chr$start)))
    start_in_validate_peak <- min(validate_peak_chr$start)
    if (min(sample_peak_chr_sort$start) > start_in_validate_peak){
      count_big_start_sample <- count_big_start_sample + 1
      close_start_site_in_validate <- validate_peak_chr$start[validate_peak_chr$start > min(sample_peak_chr_sort$start)]
      find_sample_start_in_validate <- min(close_start_site_in_validate)
      sample_start_in_validate_index <- which(validate_peak_chr$start == find_sample_start_in_validate)
      validate_peak_chr_in_sample <- validate_peak_chr[sample_start_in_validate_index: length(validate_peak_chr$start),]
      for (l in 1 : length(sample_peak_chr_sort$start)){
        distance_sample_to_validate <- sample_peak_chr_sort$start[l] - validate_peak_chr_in_sample$start  
        index_in_validate <- which(abs(distance_sample_to_validate) == min(abs(distance_sample_to_validate)))
        
        result <- data.frame(chr = de_chr_list[i], 
                             start_in_sample = sample_peak_chr_sort$start[l], 
                             end_in_sample = sample_peak_chr_sort$end[l], 
                             start_in_validate = validate_peak_chr_in_sample$start[index_in_validate],
                             end_in_validate = validate_peak_chr_in_sample$end[index_in_validate],
                             start_to_start_distance = distance_sample_to_validate[index_in_validate] ) %>% rbind(result,.)
      }
    }
    else if (min(sample_peak_chr_sort$start) < start_in_validate_peak) {
      count_small_start_sample <- count_small_start_sample + 1
      tail(sample_peak_chr_sort[sample_peak_chr_sort$start < start_in_validate_peak,], n=1)
      find_validate_start_in_sample <- max(sample_peak_chr_sort$start[sample_peak_chr_sort$start < start_in_validate_peak])
      validate_start_in_sample_index <- which(sample_peak_chr_sort$start == find_validate_start_in_sample)
      for (j in validate_start_in_sample_index : length(sample_peak_chr_sort$seqnames)) {
        distance_sample_to_validate <- as.numeric(sample_peak_chr_sort$start[j]) - as.numeric(validate_peak_chr$start)
        index_in_validate <- which(abs(distance_sample_to_validate) == min(abs(distance_sample_to_validate)))
        
        result <- data.frame(chr = de_chr_list[i],
                             start_in_sample = sample_peak_chr_sort$start[j],
                             end_in_sample = sample_peak_chr_sort$end[j],
                             start_in_validate = validate_peak_chr$start[index_in_validate],
                             end_in_validate = validate_peak_chr$end[index_in_validate],
                             start_to_start_distance = distance_sample_to_validate[index_in_validate]) %>% rbind(result,.)
      }
    }
  }
  message(paste("there are", count_big_start_sample, 'chromosomes that have bigger start sites in sample than in validate'))
  message(paste("there are", count_small_start_sample, 'chromosomes that have smaller start sites in sample than in validate'))
  return (result) 
  
}


# function to write top peak
find_thredshold <- function(input_peak, rep_number) {
  rep <- input_peak[[rep_number]]
  rpkm_sort_rep <- rep[order(-rep$RPKM),]
  plot(rpkm_sort_rep$RPKM)
}

filter_rep <- function(input_peak, rpkm_thredshold, rep_number) {
  rep <- input_peak[[rep_number]]
  rpkm_sort_rep <- rep[order(-rep$RPKM),]
  rpkm_rep_filter <- rpkm_sort_rep[rpkm_sort_rep$RPKM > rpkm_thredshold,]
  print (rpkm_rep_filter)
}

# write top peak from diffBind model contrast. apply for sample with 4 replicates 
write_top_peak <- function(input_peak, rpkm_thredshold, save_name) {
  rpkm_rep1_filter <- filter_rep(input_peak, rpkm_thredshold, 1)
  rpkm_rep2_filter <- filter_rep(input_peak, rpkm_thredshold, 2)
  rpkm_rep3_filter <- filter_rep(input_peak, rpkm_thredshold, 3)
  rpkm_rep4_filter <- filter_rep(input_peak, rpkm_thredshold, 4)
  combine_top <- rbind(rpkm_rep1_filter, rpkm_rep2_filter, rpkm_rep3_filter, rpkm_rep4_filter)
  combine_top_sub <- combine_top[,c(1:3,5)]
  combine_top_unique <- combine_top_sub[!duplicated(combine_top_sub),]
  print( paste("there are", length(combine_top_unique$Chr), "sites with RPKM higher than the thredshold"))
  write.table(combine_top_unique, paste0(save_name,'_top_sites_RPKM_thredshold_', rpkm_thredshold, '.bed'), sep = '\t', quote=F, row.names=F, col.names=F)
  print (combine_top_unique)
  }