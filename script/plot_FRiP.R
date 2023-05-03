# plot fraction of reads in peak (signal-to-noise ratio)
#Nhung 03 05 2023

# get paths from environment
fig_dir <- Sys.getenv("FIGURE_DIR_VARIABLE")
file_dir <- Sys.getenv("FRIP_DIR_VARIABLE")

# read data 
read_in_peak <- read.table(paste0(file_dir,'/read_in_peak.txt'))
read_in_peak

total_reads <- read.table(paste0(file_dir,'/total_reads.txt'))
total_reads

sample_IDs <- read.table(paste0(file_dir,'/sample_IDs.txt'))
sample_IDs

# calculate frip

# calculate frip

frip = read_in_peak/total_reads

frip_dataframe <- data.frame(sample_ID = t(sample_IDs), total_read = as.numeric(t(total_reads)), value = as.numeric(t(frip)))
# frip_dataframe
reso <- 600
png(filename = paste0(fig_dir, "/", date(),"FRIP.png"), width = 1200 * reso/72, height = 700 * reso/72, units ="px", res = reso)
par(mar = c(5,15,5,5))
par(mfrow=c(1,2))

figA <- barplot(frip_dataframe$total_read, main = "A. total reads",xlab = "Reads",
                ylab = "",
                xlim =c(0,max(frip_dataframe$total_read) + 0.1) ,
                names.arg = frip_dataframe$sample_ID,
                horiz = TRUE, las=1, col= 	'darkslateblue')

figB <- barplot(frip_dataframe$value, main = "B. Fraction of reads in peaks (FRiP)",
                xlab = "FRiP score",
                ylab = "",
                xlim =c(0,max(frip_dataframe$value) + 0.1) ,
                names.arg = frip_dataframe$sample_ID,
                horiz = TRUE, las=1, col= 	'cornflowerblue')
dev.off()

