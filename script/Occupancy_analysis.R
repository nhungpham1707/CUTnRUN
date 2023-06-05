# Nhung 05 06 2023
library(DiffBind)
library(VennDiagram)
library(tidyverse)
library(dplyr)

setwd('/home/pmc_research/npham/PROJECTS/CUTnRUN_Maroussia/script')
source('/home/pmc_research/npham/PROJECTS/CUTnRUN_Maroussia/script/Identify_DE_peaks_functions.R')

peak_path <- '/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/peak_s3norm_antiTfe3_new_control/'

tfe3= c("SCC-ChIC-PMC-DRO-T1", "SCC-ChIC-PMC-DRO-T5", "bulkChIC-PMC-DRO-016", "SCC-bulkChIC-PMC-DRO-008")
luciferase= c( "SCC-ChIC-PMC-DRO-L1", "SCC-ChIC-PMC-DRO-L5", "bulkChIC-PMC-DRO-014", "SCC-bulkChIC-PMC-DRO-005")
fusion= c("SCC-ChIC-PMC-DRO-F1", "SCC-ChIC-PMC-DRO-F5", "bulkChIC-PMC-DRO-015", "SCC-bulkChIC-PMC-DRO-002")

read_reps_peak_file <- function(path, sample_list, number_of_sites){
    rep1 <- read.csv(paste0(path, '/', sample_list[1], '_normalize.narrowPeak'), sep = '\t', header = FALSE)
    rep2 <- read.csv(paste0(path, '/', sample_list[2], '_normalize.narrowPeak'), sep = '\t', header = FALSE)
    rep3 <- read.csv(paste0(path, '/', sample_list[3], '_normalize.narrowPeak'), sep = '\t', header = FALSE)
    rep4 <- read.csv(paste0(path, '/', sample_list[4], '_normalize.narrowPeak'), sep = '\t', header = FALSE)
    
    rep1 <- rep1[2:length(rep1[,1]),c(1:3,5)]
    colnames(rep1) <- c('chr', 'start', 'end', 'score')
    rep1$mergeChrStart <- paste0(rep1$chr, rep1$start)
    rep2 <- rep2[2:length(rep2[,1]),c(1:3,5)]
    colnames(rep2) <- c('chr', 'start', 'end', 'score')
    rep2$mergeChrStart <- paste0(rep2$chr, rep2$start)
    rep3 <- rep3[2:length(rep3[,1]),c(1:3,5)]
    colnames(rep3) <- c('chr', 'start', 'end', 'score')
    rep3$mergeChrStart <- paste0(rep3$chr, rep3$start)
    rep4 <- rep4[2:length(rep4[,1]),c(1:3,5)]
    colnames(rep4) <- c('chr', 'start', 'end', 'score')
    rep4$mergeChrStart <- paste0(rep4$chr, rep4$start)
    
    all_rep <- rbind(rep1, rep2, rep3, rep4)
    print(paste('Combine all replicates resulted in', length(unique(all_rep$mergeChrStart)), 'unique sites'))
    all_rep <- all_rep[order(-all_rep$score),] 
    all_rep_unique <- unique( all_rep[ , c('mergeChrStart') ] )
    print (all_rep_unique[1:10,])
    top_sites <- all_rep_unique[1:number_of_sites,]
    print(paste("there are", length(unique(top_sites$mergeChrStart)), "unique top sites"))
    result_list <- list(top_sites, rep1, rep2, rep3, rep4)
    return (result_list)
}

tfe3_peak <- read_reps_peak_file(peak_path, tfe3, 2000)
rep1 <- tfe3_peak[[1]]
rep2 <- tfe3_peak[[2]]
rep3 <- tfe3_peak[[3]]
rep4 <- tfe3_peak[[4]]

all_rep <- rbind(rep1, rep2, rep3, rep4)
    print(paste('Combine all replicates resulted in', length(unique(all_rep$mergeChrStart)), 'unique sites'))
    all_rep <- all_rep[order(-all_rep$score),1:3] 
    all_rep_unique <- all_rep[!duplicated(all_rep),]
length (all_rep_unique$chr)
    top_sites <- all_rep_unique[1:number_of_sites,]
    print(paste("there are", length(unique(top_sites$mergeChrStart)), "unique top sites"))

