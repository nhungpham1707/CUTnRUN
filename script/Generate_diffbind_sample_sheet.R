# Automatically generate sample sheet for diffBind analysis 
# Nhung 31 05 2023 

library(dplyr)
sample_sheet_dir <- Sys.getenv("SAMPLE_SHEET_DIR") 
bam_dir <- Sys.getenv("BAM_DIR")
peak_dir <- Sys.getenv("PEAK_DIR")
filename <- Sys.getenv("SAVE_NAME")

sample_IDs <- c("SCC-ChIC-PMC-DRO-T1", "SCC-ChIC-PMC-DRO-T5", "SCC-bulkChIC-PMC-DRO-008", "bulkChIC-PMC-DRO-016",
                "SCC-ChIC-PMC-DRO-F1", "SCC-ChIC-PMC-DRO-F5", "SCC-bulkChIC-PMC-DRO-002", "bulkChIC-PMC-DRO-015",
                "SCC-ChIC-PMC-DRO-L1", "SCC-ChIC-PMC-DRO-L5", "SCC-bulkChIC-PMC-DRO-005", "bulkChIC-PMC-DRO-014")
Tissue <- rep('kidney', times = length(sample_IDs))
Factor <- c("tfe3", "tfe3", "tfe3", "tfe3",
            'fusion','fusion', 'fusion', 'fusion',
            'luciferase', 'luciferase', 'luciferase', 'luciferase')
Treatmen <- c("round1", "round1", "round2", "round3", 
              "round1", "round1", "round2", "round3", 
              "round1", "round1", "round2", "round3")
Condition <- c("anti-tfe3-1000", "anti-tfe3-500",  "anti-tfe3-1000", "anti-tfe3-1000",
               "anti-tfe3-1000", "anti-tfe3-500", "anti-tfe3-1000", "anti-tfe3-1000",
               "anti-tfe3-1000", "anti-tfe3-500",  "anti-tfe3-1000", "anti-tfe3-1000")
Replicate <- c(1, 2, 3, 4,
               1, 2, 3, 4,
               1, 2, 3, 4)
PeakCaller <- rep("bed", times = length(sample_IDs))

# get bam directory 
bamReads <- c()
for (i in 1:length(sample_IDs)){
  bamReads <- data.frame(list.files(
    bam_dir,
    pattern = paste0(sample_IDs[i],"_rmdup_filt.bam$"),
    recursive = TRUE,
    full.names = TRUE
  )) %>% rbind (bamReads,.)
}

colnames(bamReads) <- "bamReads"

# get peak directory
peaks <- c()
for (i in 1:length(sample_IDs)){
  peaks <- data.frame(list.files(
    peak_dir,
    pattern = paste0(sample_IDs[i],".*.narrowPeak$"),
    recursive = TRUE,
    full.names = TRUE
  )) %>% rbind (peaks,.)
}

colnames(peaks) <- "Peaks"

# col_names <- c("SampleID", "Tissue", "Factor", "Condition", "Treatment", "Replicate",
               # "bamReads", "Peaks","PeakCaller")
# write to file
sample_sheet <- data.frame(SampleID = sample_IDs,
                           Tissue = Tissue,
                           Factor = Factor,
                           Condition = Condition,
                           Treatmen = Treatmen,
                           Replicate = Replicate,
                           bamReads = bamReads,
                           Peaks = peaks,
                           PeakCaller = PeakCaller)
write.csv(sample_sheet, paste0(sample_sheet_dir, '/', filename,'.csv'), row.names = FALSE)

