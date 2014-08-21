## Statistical comparisons between esets from gct files and cel files for expt. June 20th 2014. Test to see if the esets are different!

i = 2 # To pick out which experiment to run this test on
#libraries
library(affy)
library(arrayQualityMetrics)
library(limma)
library(gcrma)
#Heatmap libraries
library(pheatmap)
library(gplots)
library(pathprint)
#Correlation libraries
library(ggplot2)
#library(CCA)

#Directories & files
celdir = "/Users/admin/Projects/Pathprint/Data/celfiles/"
gctdir = "/Users/admin/Projects/Pathprint/Data/gctfiles"
setwd("/Users/admin/Dropbox/HSPH")
results.path <- "results/comparison/"
covarsfilename = "metadata/cel_zon.txt" #Zon expts from "Dropbox/HSCI_Curation/Len Zon data/Replicates_DataFiles.xls"
QCresultsDir = "results/ArrayQC/"

#Import gct data
gct<-read.delim("metadata/GCT_lenzon.txt") #can set nrows
gctlist<-gct$'Server.Path'
nidlist<-gct$'Nid'
gpllist<-gct$'GPL'
nid = nidlist[i]
gpl_platform = gpllist[i]

#Import expression data from gct file
exprs_expt<-as.matrix(read.delim(paste("/Users/admin/Projects/Pathprint/Data/gctfiles/",nid,"_",gpl_platform,".gct", sep = ''), skip = 2, row.names = 1 , header = TRUE,  as.is = TRUE))
exprsfile <- exprs_expt[,-1]

#Import CEL data and metadata; Process it
##Must have phenotype file (covars)
covarslist <- read.delim(file.path(covarsfilename), row.names = 1, sep = "\t")
covars <- covarslist[covarslist["exp_nid"] == nid,]
celFiles <- file.path(celdir, row.names(covars))
affyRaw <- ReadAffy(filenames = celFiles) # Read in cel files, create affybatch file
pData(affyRaw)<-covars #Assigning phenotype data
sampleNames(affyRaw)<-pData(affyRaw)$replicate_name #Assigning sample names
validObject(affyRaw) #check if valid object
rm(covars)
#Processing
##Background correction and normalization (using RMA/justGCRMA)
affyNorm <- rma(affyRaw, background = TRUE, normalize = TRUE)
#affyNorm <- justGCRMA(filenames = celFiles, samplenames = pData(affyRaw)$replicate_name)

eset<-exprs(affyNorm)

#affyNorm is normalized eset from cel file
#eset is the matrix derived from affyNorm
#exprsfile is gct eset read in as a matrix

# OUTPUT RESULTS TO
cel_eset <- eset
gct_eset <- exprsfile
colnames(cel_eset) <- paste("cel", colnames(cel_eset), sep= "_")
colnames(gct_eset) <- paste("gct", colnames(gct_eset), sep= "_")
merged.eset<- merge(cel_eset, gct_eset, by = "row.names", all = TRUE)
rownames(merged.eset) <- merged.eset[,1]
merged.eset <- as.matrix(merged.eset[,-1])
#run pathprints
mergedeset.fingerprint <-exprs2fingerprint(exprs = merged.eset, platform = as.character(gpl_platform), species = "Danio rerio", progressBar = TRUE)
# Write fingerprint to file
dir.create(file.path(paste(results.path, nid, "/", sep="")), showWarnings = FALSE)
write.table(mergedeset.fingerprint, file=paste(results.path, nid, "/", nid, 
                                         "_rma_jointpathprint.txt", sep=""),
            sep="\t", quote=FALSE, row.names=TRUE)
# HEATMAP
merged.fprint <- mergedeset.fingerprint

# Finding the na subset
emptydf <- subset(merged.fprint, is.na(merged.fprint[,1]))
head(emptydf)
# Excluding the na data 
merged.fprint <- na.omit(merged.fprint)
# Heatmap with sds. Derived from John (github) and Shannan's codes.
pdf(file=paste(results.path, "/","rma_heatmap_comparison.pdf", sep=""))
for (sd in seq(0.5,0.9, 0.1)){
  heatmap.data <- merged.fprint[apply(merged.fprint, 1, sd) > sd,]
  if(!is.null(heatmap.data) && !is.null(row.names(heatmap.data))){
    heatmap.data<- heatmap.data[complete.cases(heatmap.data),]
    if(nrow(heatmap.data) > 0){
      pheatmap(heatmap.data, 
               col=c("blue", "white", "red"), 
               legend_breaks=c(-1,0,1), 
               main=paste("All Samples, sd=", sd, sep=""), 
               show_colnames=T, 
               fontsize_row=5,
               fontsize_col = 5,
               cluster_cols=TRUE,
               show_rownames=T)
    }}}
write.table(heatmap.data, file=paste(results.path, nid, "/", nid, 
                                               "_jointheatmap.txt", sep=""),
            sep="\t", quote=FALSE, row.names=TRUE)
dev.off()  
