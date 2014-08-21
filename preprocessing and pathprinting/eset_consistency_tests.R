## Statistical comparisons between esets from gct files and cel files for expt. June 20th 2014. Test to see if the esets are different!

i = 2 # To pick out which experiment to run this test on
#libraries
library(affy)
library(arrayQualityMetrics)
library(limma)
library(gcrma)
#Correlation libraries
library(ggplot2)
#library(CCA)

#Directories & files
celdir = "/Users/admin/Projects/Pathprint/Data/celfiles/"
gctdir = "/Users/admin/Projects/Pathprint/Data/gctfiles"
setwd("/Users/admin/Dropbox/HSPH")
results.path <- "results/correlation/"
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
##Background correction and normalization (using RMA)
affyNorm <- rma(affyRaw, background = TRUE, normalize = TRUE)
#affyNorm <- rma(affyRaw)

eset<-exprs(affyNorm)

#affyNorm is normalized eset from cel file
#eset is the matrix derived from affyNorm
#exprsfile is gct eset read in as a matrix
# OUTPUT RESULTS TO
pdf(file=paste(results.path, "/","rma_normalized_comparison.pdf", sep=""))

#Run correlation
gcteset <- t(apply(exprsfile,1,as.numeric))
celeset <- t(apply(eset,1,as.numeric))

#plot(var(gcteset[,1], celeset[,1], na.rm = TRUE))
plot(cov(gcteset,celeset), main = "Covariance run for expression matrices Expt 7017")
plot(cor(gcteset, celeset), main = "Correlation run for expression matrices Expt 7017")
#plot(cor(gcteset[,3:4], celeset[,3:4]))
#compare.matrix(gcteset, celeset, nlevels = 5)

sample1gct <- as.vector(gcteset[,1])
sample1cel <- as.vector(celeset[,1])
sample1gct.rank <- rank(sample1gct)
sample1cel.rank <- rank(sample1cel)
plot(cor(sample1gct.rank, sample1cel.rank), main = "Correl for ranked genes vector of Repl 1 Expt 7017")
plot(cov(sample1gct.rank, sample1cel.rank), main = "Covar for ranked genes vector of Repl 1 Expt 7017")
boxplot(sample1cel.rank, sample1gct.rank, names = c("gct", "cel"), main = "Ranked genes vector of replicate 1 Expt 7017")
wilcox.test(sample1gct, sample1cel)
#correl = matcor(gcteset, celeset)
#img.matcor(correl)
#img.matcor(correl, type = 2)
dev.off()
