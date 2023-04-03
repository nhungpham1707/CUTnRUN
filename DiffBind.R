# This script uses DiffBind to identify differentially bound sites and generate visualization 

# Nhung 29 03 2023

library(DiffBind)
library(Gviz)
# library(BiocParallel)
# multicoreParam <- MulticoreParam(workers = 3)
# register(multicoreParam)
# registered()
# importFrom(BiocParallel, bplapply)
# importFrom(BiocParallel, mclapply)

# read in sample. The sample sheet file is created manually.

samples <- read.csv('sample_sheet2.csv', sep=";")

head(samples)

# create dba object for further analysis
dbObj <- dba(sampleSheet=samples)

##bParallel has to be set to F else it will fail in cluster hPC environment 
# generate peak count - this step takes ~ 30 mins - 1hour 
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, bParallel =F) 

# not sure what they do
for(list in 1:length(dbObj[["peaks"]])){dbObj[["peaks"]][[list]]$cReads <- 1}
for(list in 1:length(dbObj[["peaks"]])){dbObj[["peaks"]][[list]]$cRPKM <- 1}
for(list in 1:length(dbObj[["peaks"]])){dbObj[["peaks"]][[list]]$RPKM <- dbObj[["peaks"]][[list]]$RPKM + 1}
for(list in 1:length(dbObj[["peaks"]])){dbObj[["peaks"]][[list]]$Reads <- dbObj[["peaks"]][[list]]$Reads + 1}

info <- dba.show(dbObj)
libsizes <- cbind(LibReads=dbObj$Reads, FRiP=dbObj$FRiP,
                  PeakReads=round(dbObj$Reads * dbObj$FRiP))
rownames(libsizes) <- info$ID

# plot PCA
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
dba.plotPCA(dbObj,  attributes=DBA_CONDITION, label=DBA_ID)

# make correlation plot
plot(dbObj)

##Normalize the data by library sizes
B.DBcount <- dba.normalize(dbObj)
norm <- dba.normalize(B.DBcount, bRetrieve=TRUE)
##Differential peaks not correcting for potential batch effect
dbOri <- dba.contrast(dbObj,design = "~Factor")
dbOri <- dba.analyze(dbOri)
DBresultsOri <- dba.report(dbOri)
plot(dbOri, contrast=1)
dba.plotPCA(dbOri, contrast=1, label = DBA_ID)
dba.show(dbOri, bContrasts=TRUE)