# this script calculate mean from replicates in merge bedgraph files output from Merge_bedgraphs.sh

# Nhung 24 05 2023

res_dir = Sys.getenv("RESULT_DIR_VARIABLE")
save_name = Sys.getenv("SAVE_NAME_VARIABLE")
file_path = Sys.getenv("FILE_PATH_VARIABLE")

merge_file <- read.csv(file_path, header= FALSE, sep = '\t')
colnames(merge_file) <- c('chr', 'start', 'end', 'rep1', 'rep2', 'rep3', 'rep4')

message("head of merge_file is ")
head(merge_file)

message("unique merge_file chromosome")
unique(merge_file$chr) 

message("total number of sites")
length(merge_file$chr) 

merge_file$mean <- rowMeans(merge_file[,4:7])

message("head merge_file with mean column added is")
head(merge_file, n = 30)

merge_file_bed <- merge_file[,c(1:3,8)]
# keep only lines that have seq name start with chr
chr_merge_file_bed <- merge_file_bed[startsWith(merge_file_bed$chr, 'chr'),]

message("the number of seqs that have name does not start with chr is")
length(setdiff(merge_file_bed$chr, chr_merge_file_bed$chr)) 

# write to file
write.table(chr_merge_file_bed, file=paste0(res_dir,'/', Sys.Date(), save_name), sep="\t", quote=F, row.names=F, col.names=F)