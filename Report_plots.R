# report and visualize result from alignment, duplication, and peaks
# Script is adapted from https://yezhengstat.github.io/CUTTag_tutorial/#32_Report_sequencing_mapping_summary_[required]
# Nhung 06 04 2023

library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)

################# Alignment report ###############
#################################################
# Define paths where the reports from alignment and duplication are
align_dir = "/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/alignment"
rm_dup_dir = "/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/rm_dup"
peak_dir =  "/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/peakCalling"
# Load samples and generate alignment result

sample_info = read.csv("sample_info.csv", sep = ";")
sample_info
alignResult = c()

for (i in 1:length(sample_info$SampleID)){
alignRes = read.table(file = paste(align_dir, "/", sample_info$SampleID[i], "/", sample_info$SampleID[i],".txt", sep =""), header = FALSE, fill = TRUE)

alignResult = data.frame(SampleID = sample_info$SampleID[i],
                         Antibody = sample_info$antibody[i],
                         Factor = sample_info$Factor[i],
                         Replicate = sample_info$Replicate[i],
                          SequencingDepth = alignRes$V1[1] %>% as.character %>% as.numeric,
                          MappedFragNum_hg38 = alignRes$V1[4] %>% as.character %>% as.numeric + alignRes$V1[5] %>% as.character %>% as.numeric,
                         AlignmentRate_hg38 = as.numeric(substr(alignRes$V1[nrow(alignRes)], 1, nchar(as.character(alignRes$V1[nrow(alignRes)]))-1)))  %>% rbind(alignResult, .)
}
alignResult$Factor = factor(alignResult$Factor, levels = unique(sample_info$Factor))

alignResult %>% mutate(AlignmentRate_hg38 = paste0(AlignmentRate_hg38, "%"))

# Visualize the seq depth and alignment rate
# figA = alignResult %>% ggplot(aes(x = Factor, y = SequencingDepth/1000000, fill = Factor)) +
#   geom_boxplot() +
#   geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
#   scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
#   scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
#   theme_bw(base_size = 18) +
#   ylab("Sequencing Depth per Million") +
#   xlab("") + 
#   ggtitle("A. Sequencing Depth")
# figB = alignResult %>% ggplot(aes(x = Factor, y = MappedFragNum_hg38/1000000, fill = Factor)) +
#   geom_boxplot() +
#   geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
#   scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
#   scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
#   theme_bw(base_size = 18) +
#   ylab("Mapped Fragments per Million") +
#   xlab("") +
#   ggtitle("B. Alignable Fragment (hg38)")

# figC = alignResult %>% ggplot(aes(x = Factor, y = AlignmentRate_hg38, fill = Factor)) +
#   geom_boxplot() +
#   geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
#   scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
#   scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
#   theme_bw(base_size = 18) +
#   ylab("% of Mapped Fragments") +
#   xlab("") +
#   ggtitle("C. Alignment Rate (hg38)")

# plot horizontal
figA = alignResult %>% ggplot(aes(x = SequencingDepth/1000000, y = Factor , fill = Factor)) +
  geom_boxplot() +
  geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("") +
  xlab("Sequencing Depth per Million") + 
  ggtitle("A. Sequencing Depth")
figB = alignResult %>% ggplot(aes(x = MappedFragNum_hg38/1000000, y = Factor, fill = Factor)) +
  geom_boxplot() +
  geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("") +
  xlab("Mapped Fragments per Million") +
  ggtitle("B. Alignable Fragment (hg38)")


figC = alignResult %>% ggplot(aes(x = AlignmentRate_hg38, y = Factor , fill = Factor)) +
  geom_boxplot() +
  geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("") +
  xlab("% of Mapped Fragments") +
  ggtitle("C. Alignment Rate (hg38)")

reso <- 1200

png(filename = "alignment_report.png" , width = 1000 * reso/72, height = 700 * reso/72, units ="px", res = reso)
plot = ggarrange(figA, figB, figC, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")
annotate_figure(plot, top = text_grob("Alignment report", 
                                      color = "black", face = "bold", size = 30))
dev.off()
##################Duplication rate ######################
#########################################################
dupResult = c()

for (i in 1:length(sample_info$SampleID)){
  dupRes = read.table(file = paste(rm_dup_dir, "/", sample_info$SampleID[i],"/", sample_info$SampleID[i],"_marked_dup_metrics.txt", sep =""), header = TRUE, fill = TRUE)
#dupRes = read.table("SCC-ChIC-PMC-DRO-F1_marked_dup_metrics.txt", header = TRUE, fill = TRUE)
  dupRes$READ_PAIRS_EXAMINED[1]
  dupResult = data.frame(SampleID = sample_info$SampleID[i],
                       Antibody = sample_info$antibody[i],
                       Factor = sample_info$Factor[i],
                       Replicate = sample_info$Replicate[i],
                       MappedFragNum_dup_hg38 = dupRes$READ_PAIRS_EXAMINED[1] %>% as.character %>% as.numeric, 
                       DuplicationRate = dupRes$PERCENT_DUPLICATION[1] %>% as.character %>% as.numeric * 100, 
                       EstimatedLibrarySize = dupRes$ESTIMATED_LIBRARY_SIZE[1] %>% as.character %>% as.numeric) %>% mutate(UniqueFragNum = MappedFragNum_dup_hg38 * (1-DuplicationRate/100))  %>% rbind(dupResult, .)
}

  dupResult

  dupResult$Factor = factor(dupResult$Factor, levels = unique(dupResult$Factor))
  alignDupSummary = left_join(alignResult, dupResult, by = c("SampleID", "Antibody", "Factor", "Replicate")) %>% mutate(DuplicationRate = paste0(DuplicationRate, "%"))
  alignDupSummary

  # Visualize dup rate
  # fig4A = dupResult %>% ggplot(aes(x = Factor, y = DuplicationRate, fill = Factor)) +
  #   geom_boxplot() +
  #   geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  #   scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  #   scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  #   theme_bw(base_size = 18) +
  #   ylab("Duplication Rate (*100%)") +
  #   xlab("") 
  # 
  # fig4B = dupResult %>% ggplot(aes(x = Factor, y = EstimatedLibrarySize, fill = Factor)) +
  #   geom_boxplot() +
  #   geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  #   scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  #   scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  #   theme_bw(base_size = 18) +
  #   ylab("Estimated Library Size") +
  #   xlab("") 
  # 
  # fig4C = dupResult %>% ggplot(aes(x = Factor, y = UniqueFragNum, fill = Factor)) +
  #   geom_boxplot() +
  #   geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  #   scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  #   scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  #   theme_bw(base_size = 18) +
  #   ylab("# of Unique Fragments") +
  #   xlab("")
  
  # plot horizontal 
  fig4A = dupResult %>% ggplot(aes(x = DuplicationRate, y = Factor, fill = Factor)) +
    geom_boxplot() +
    geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("") +
    xlab("Duplication Rate (*100%)") 
  
  fig4B = dupResult %>% ggplot(aes(x = EstimatedLibrarySize, y = Factor , fill = Factor)) +
    geom_boxplot() +
    geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("") +
    xlab("Estimated Library Size") 
  
  fig4C = dupResult %>% ggplot(aes(x = UniqueFragNum, y = Factor, fill = Factor)) +
    geom_boxplot() +
    geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("") +
    xlab("# of Unique Fragments")
  
  reso <- 1200
  
  png(filename = "duplication_rate_report.png" , width = 1300 * reso/72, height = 700 * reso/72, units ="px", res = reso)
  plot = ggarrange(fig4A, fig4B, fig4C, ncol = 3, common.legend = TRUE, legend="bottom")
  annotate_figure(plot, top = text_grob("Sequence duplication report", 
                                        color = "black", face = "bold", size = 30))
  dev.off()
  
  ## Access mapped fragment size distribution
  fragLen = c()
  for (i in 1:length(sample_info$SampleID)){
    fragLen = read.table(paste0(align_dir, "/", sample_info$SampleID[i], "/", sample_info$SampleID[i],"_fragmentLen.txt"), header = FALSE) %>% mutate(fragLen = V1 %>% as.numeric, fragCount = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)), 
                                                                                                                                     SampleID = sample_info$SampleID[i],                                                    
                Antibody = sample_info$antibody[i],
                Factor = sample_info$Factor[i],
                Replicate = sample_info$Replicate[i],) %>% rbind(fragLen, .) 
  }

  
fragLen$sampleID = factor(fragLen$SampleID, levels = unique(fragLen$SampleID) )
fragLen$Factor = factor(fragLen$Factor, levels = unique(fragLen$Factor))
fragLen$Replicate = factor(fragLen$Replicate, levels = unique(fragLen$Replicate))

fragLen
## Generate the fragment size density plot (violin plot)
# customize the y limit (scale_y_continuous(limits= c(0,600),breaks = seq(0, 600, 50)))
#. try plot with bigger range to see where the seq line and then define the ranges for plotting

fig5A = fragLen %>% ggplot(aes(x = Factor, y = fragLen, weight = Weight, fill = Factor, group = Factor)) +
  geom_violin(bw = 5) +
  scale_y_continuous(limits= c(0,600),breaks = seq(0, 600, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length") +
  xlab("")

reso <- 1200

png(filename = "alignment_length_report.png" , width = 1300 * reso/72, height = 700 * reso/72, units ="px", res = reso)
plot <- ggarrange(fig5A, ncol = 1)
annotate_figure(plot, top = text_grob("Alignment fragment length", 
                                      color = "black", face = "bold", size = 30))
dev.off()


# fig5A = fragLen %>% ggplot(aes(x = Factor, y = fragLen, weight = Weight, fill = Factor, group = Factor)) +
#   geom_violin(bw = 5) +
#   scale_y_continuous(limits= c(0,600),breaks = seq(0, 600, 50)) +
#   scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
#   scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
#   theme_bw(base_size = 20) +
#   ggpubr::rotate_x_text(angle = 20) +
#   ylab("Fragment Length") +
#   xlab("")
# 
# fig5B = fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = as.factor(Factor), group = Factor, linetype = Replicate)) +
#   geom_line(size=1) +
#   scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
#   theme_bw(base_size = 20) +
#   xlab("Fragment Length") +
#   ylab("Count") +
#   coord_cartesian(xlim = c(0, 1000))
# 
# fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = as.factor(Factor), group = interaction(Factor, Replicate)))+
#   geom_line(size=1) +
#   scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
#   theme_bw(base_size = 20) +
#   xlab("Fragment Length") +
#   ylab("Count") +
#   coord_cartesian(xlim = c(0, 600))
# 
# ggarrange(fig5A, fig5B, ncol = 2)

#### peaks statistic
# the summit bed file This is in BED format, which contains the peak summits locations for every peak. The 5th column in this file is -log10pvalue the same as NAME_peaks.bed https://miawang113.wordpress.com/2019/01/30/output-files-from-macs-2/
peakN = c()
peakWidth = c()
# peakInfo = read.table(paste0(peak_dir, "/", sample_info$SampleID[i], "/narrow/", sample_info$SampleID[i], "_paired_control_summits.bed"), header = FALSE, fill = TRUE) %>% mutate(width = abs(V3-V2))

for (i in 1:length(sample_info$SampleID)){
peakInfo = read.table(paste0(peak_dir, "/", sample_info$SampleID[i], "/narrow/", sample_info$SampleID[i], "_paired_control_peaks.narrowPeak"), header = FALSE, fill = TRUE) %>% mutate(width = abs(V3-V2))
peakN = data.frame(peakN = nrow(peakInfo),
                   Factor = sample_info$Factor[i],
                   Replicate = sample_info$Replicate[i],
                   SampleID = sample_info$SampleID[i]) %>% rbind(peakN,.)
peakWidth = data.frame(width = peakInfo$width, 
                       Factor = sample_info$Factor[i],
                       Replicate = sample_info$Replicate[i],
                       sampleID = sample_info$SampleID[i]) %>% rbind(peakWidth, .)
}
peakN
peakWidth

fig7A = peakN %>% ggplot(aes(x = Factor, y = peakN, fill = Factor)) +
  geom_boxplot() +
  geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("")

fig7B = peakWidth %>% ggplot(aes(x = Factor, y = width, fill = Factor)) +
  geom_violin() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_y_continuous(trans = "log", breaks = c(400, 3000, max(peakWidth$width))) +
  theme_bw(base_size = 18) +
  ylab("Width of Peaks") +
  xlab("")

reso <- 1200

png(filename = "peaks_number_report.png" , width = 1300 * reso/72, height = 700 * reso/72, units ="px", res = reso)
plot <- ggarrange(fig7A, fig7B, ncol = 2, nrow=1, common.legend = TRUE, legend="bottom")

annotate_figure(plot, top = text_grob("Peaks numbers", 
                                      color = "black", face = "bold", size = 30))
dev.off()

