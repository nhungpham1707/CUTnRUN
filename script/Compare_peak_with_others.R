# IN this script, DE site from s3norm remove low depth histone sample were
# overlap with histone, acetyl and other validation samples to check if there
# is any promoter, enhancer or they are real sfpq fusion 

# Nhung 25 05 2023 
library(dplyr)
setwd('~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/Overlap_with_other_samples/')
source('~/Dropbox/POSTDOC_MAXIMA/PROJECTS/CUT_RUN/FROM_HPC/test_diffBind_HPC/Identify_DE_peaks.R')
# all significant sites 
load('s3norm_all_sig_sites_GO_enrich_ensembl_annotation.RData')

# extract peak with fold > 2
peak_fold2 <- peak_annotation[abs(peak_annotation$Fold) > 2, ] 
head(peak_fold2)
length(peak_fold2$seqnames) # 2653

order <- peak_fold2[order(peak_fold2$Fold),]
order[,c(1:3, 12:21)]
colnames(order)
plot(order$seqnames)
### overlap with H3K27ac
# read in validation samples
fusion_ac <- read.csv("peak_files/SCC-bulkChIC-PMC-DRO-025_normalize.narrowPeak", 
                      sep = '\t', header = FALSE)
head(fusion_ac)
fusion_ac_df <- fusion_ac[,c(1:3,5)]
colnames(fusion_ac_df) <- c('chr', 'start', 'end', 'score')
head(fusion_ac_df)

luc_ac <- read.csv("peak_files/SCC-bulkChIC-PMC-DRO-022_normalize.narrowPeak", 
                   sep = '\t', header = FALSE)
head(luc_ac)
luc_ac_df <- luc_ac[,c(1:3,5)]
colnames(luc_ac_df) <- c('chr', 'start', 'end', 'score')
head(luc_ac_df)

# overlay with sfpq samples
fusion_sfpq <- read.csv("peak_files/SCC-bulkChIC-PMC-DRO-023_normalize.narrowPeak", 
                        sep = '\t', header = FALSE)
fusion_sfpq_df <- fusion_sfpq[,c(1:3,5)]
colnames(fusion_sfpq_df) <- c('chr', 'start', 'end', 'score')
head(fusion_sfpq_df)
length(fusion_sfpq_df$chr) # 2772 

luc_sfpq <- read.csv("peak_files/SCC-bulkChIC-PMC-DRO-020_normalize.narrowPeak", 
                     sep = '\t', header = FALSE)
luc_sfpq_df <- luc_sfpq[,c(1:3,4)]
colnames(luc_sfpq_df) <- c('chr', 'start', 'end', 'score')
head(luc_sfpq_df)
length(luc_sfpq_df$chr) # 5156

sample_peak <- peak_fold2
validate_peak <- fusion_sfpq_df

overlay_peak_result <- overlay_peak_start_site(sample_peak, validate_peak)
print(paste('the biggest distance between sample and validate is', max(abs(overlay_peak_result$start_to_start_distance))))
print(paste('the smallest distance between sample and validate is', min(abs(overlay_peak_result$start_to_start_distance))))
plot(overlay_peak_result$start_to_start_distance)
plot(overlay_peak_result$start_to_start_distance, ylim= c(-50000,50000))

distance_threshold <- c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
peak_number <- c()
for (i in 1:length(distance_threshold)){
overlay_peak_result_filter <- overlay_peak_result[abs(overlay_peak_result$start_to_start_distance) < distance_threshold[i],]
    
overlay_peak_result_filter[order(abs(overlay_peak_result_filter$start_to_start_distance)),]
plot(overlay_peak_result_filter$start_to_start_distance)
peak_number <- data.frame(distance_threshold = distance_threshold[i],
                          number_peak = length(overlay_peak_result_filter$chr))%>% rbind(peak_number,.)
}

plot(peak_number$distance_threshold, peak_number$number_peak)
# # chekc why did not recover chr2 30979825  on validate
# [1] "chr"                     "start_in_sample"         "end_in_sample"           "start_in_validate"       "end_in_validate"        
# [6] "start_to_start_distance"
# 411    2        30979825      30980224          30980000        30980200                    -175

sample_peak <- peak_fold2
validate_peak <- fusion_ac_df
overlay_ac <- overlay_peak_start_site(sample_peak, validate_peak )

distance_threshold <- seq(2000, 1000000, by=1000)
peak_number <- c()
for (i in 1:length(distance_threshold)){
  overlay_peak_result_filter <- overlay_ac[abs(overlay_ac$start_to_start_distance) < distance_threshold[i],]
  
 peak_number <- data.frame(distance_threshold = distance_threshold[i],
                            number_peak = length(overlay_peak_result_filter$chr))%>% rbind(peak_number,.)
}

plot(peak_number$distance_threshold, peak_number$number_peak)


