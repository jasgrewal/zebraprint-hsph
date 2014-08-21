# Takes in a batch of cel files (based on expt id (nid)) and affixes phenotype info, generates expression set, normalizes (previous runs by justGCRMA), runs QC,
# Then runs pathway fingerprint on expression set, generates pathway signatures txt file and a heatmap at different sd values

#libraries
library(affy)
library(arrayQualityMetrics)
library(limma)
library(gcrma)

#pathprint libraries
library(pheatmap)
library("Biobase")
library(pathprint)
library(gplots)

#Directories & files
celdir = "/Users/admin/Projects/Pathprint/Data/celfiles/"
setwd("/Users/admin/Dropbox/HSPH")
results.path <- "results/pathprints/"
covarsfilename = "metadata/cel_zon.txt" #Zon expts from "Dropbox/HSCI_Curation/Len Zon data/Replicates_DataFiles.xls"
QCresultsDir = "results/ArrayQC/"

#metadata
covarslist <- read.delim(file.path(covarsfilename), row.names = 1, sep = "\t")

# To get gpls
gct<-read.delim("metadata/GCT_lenzon.txt") #can set nrows
nidlist<-gct$'Nid'
gpllist<-gct$'GPL'

for(i in 1:length(nidlist)){ #Iterating through expt id's and fingerprinting each one

 #Import data and metadata
 ##Must have phenotype file (covars)
 nid <- nidlist[i]
 gpl_platform = (gct[gct['Nid'] == nid,])$'GPL'[1]
 dir.create(file.path(paste(results.path, nid, "/", sep="")), showWarnings = FALSE)
 dir.create(file.path(paste(results.path, nid, "/ArrayQC", sep="")), showWarnings = FALSE)
 QCresultsDir = paste(results.path, nid, "/ArrayQC", sep="")
 
 covars <- covarslist[covarslist["exp_nid"] == nid,]
 celFiles <- file.path(celdir, row.names(covars))
 affyRaw <- ReadAffy(filenames = celFiles) # Read in cel files, create affybatch file

 pData(affyRaw)<-covars #Assigning phenotype data
 sampleNames(affyRaw)<-pData(affyRaw)$replicate_name #Assigning sample names
 validObject(affyRaw) #check if valid object

 #rm(covars)
 #Processing
 ##Raw Data QC (using arrayQualityMetrics)
 normtype = "gcrma" # CHANGE HERE FOR NORMALIZATION STEP TOO (gcrma? rma?)

 ##!##arrayQualityMetrics(expressionset=(affyRaw), outdir=file.path(QCresultsDir, paste(nid, "report_raw", sep=".")), force=TRUE, do.logtransform=TRUE) # intgroup Was group in covers.desc
 ##Background correction and normalization (using RMA/gcRMA)
 #NORMALIZATION - rma or GCRMA
 # IF RMA: 
 #affyNorm <- rma(affyRaw, background = TRUE, normalize = TRUE) 
 setwd("/")
 affyNorm <- justGCRMA(filenames = celFiles, sampleNames = sampleNames(affyRaw))
 setwd("/Users/admin/Dropbox/HSPH")
 
 #Normalized data QC (using arrayQualityMetrics)
 eset<-exprs(affyNorm)
 ##!##arrayQualityMetrics(expressionset=(affyNorm), outdir=file.path(QCresultsDir, paste(nid,normtype, "report", sep=".")), force=TRUE, do.logtransform=FALSE)

 # PATHPRINT

 #read gct file for getting platform of respective expt for fingerprint
 gct<-read.delim("metadata/GCT_lenzon.txt")

 #Fingerprint new expression data
 eset <- exprs(affyNorm)
 expt.fingerprint <- exprs2fingerprint(exprs = eset, platform = as.character(gpl_platform), species = "Danio rerio", progressBar = TRUE)
 # Write eset and fingerprint to file
 write.table(as.matrix(eset), file=paste(results.path, nid, "/", nid, 
                                          "_cel_",normtype, "_eset.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)
 write.table(expt.fingerprint, file=paste(results.path, nid, "/", nid, 
        "_cel_",normtype, "_pathprint.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)

 # Finding the na subset
 emptydf <- subset(expt.fingerprint, is.na(expt.fingerprint[,1]))
# head(emptydf)
 # Excluding the na data 
 expt.fingerprintnew <- na.omit(expt.fingerprint)

 # Heatmap distribution of pathprint expression signatures
 pdf(file=paste(results.path, nid, "/", nid,"_cel_",normtype, "_heatmap_sdweight.pdf", sep=""))
 logfile = paste(results.path, "logfile_cel.txt", sep="")
 for (sd in seq(0.5,0.8, 0.1)){
  heatmap.data <- expt.fingerprintnew[apply(expt.fingerprintnew, 1, sd) > sd,]
  if(!is.null(heatmap.data) && !is.null(row.names(heatmap.data))){
    heatmap.data<- heatmap.data[complete.cases(heatmap.data),]
    if(nrow(heatmap.data) > 0){
      pheatmap(heatmap.data, 
               col=c("blue", "white", "red"), 
               legend_breaks=c(-1,0,1), 
               main=paste("All Samples, sd=", sd, sep=""), 
               show_colnames=T, 
               fontsize_row=6, 
               fontsize_col = 7,
               cluster_cols=TRUE,
               show_rownames=T)
    }
    else{
      write(paste(nid, " for sd ", sd, " has no subset", sep = ""), file = logfile,
            append = TRUE, sep = "\n")
    }
 }
  else{
    write(paste(nid, " for sd ", sd, " has unlabelled subset\n\t\t", heatmap.data, sep = ""), file = logfile,
          append = TRUE, sep = "\n")
  }
}
dev.off()
}