#Takes in a batch of gct files (based on expt id (nid)) and affixes phenotype info
#Then runs pathway fingerprint on expression set, generates pathway signatures txt file and a heatmap at different sd values


setwd("/Users/admin/Dropbox/HSPH")

results.path <- "results/pathprints/"
#gct.path <- "/n/hsphS10/hsphfs1/hsci/sccr/sccr-prod/gct/"
#gct <- list.files(path = gct.path, pattern = "gct")
gct<-read.delim("metadta/GCT_lenzon.txt") #can set nrows
gctlist<-gct$'Server.Path'
nidlist<-gct$'Nid'
gpllist<-gct$'GPL'
head(gctlist)
length(gctlist)
metadata <- read.delim("metadata/Replicate_list.txt", check.names = F)
#47 files
library(pheatmap)
library("Biobase")
library(pathprint)
library(gplots)
for(i in 1:length(gctlist)){ #Iterating through expt id's and fingerprinting each one
  
  #select only the metadata req for this expt
  nid = nidlist[i]
  gpl_platform = gpllist[i]
  #17/06/14 note: Replicate names have now changed (use old gct with mapping to new replicate names, or pathprint from cel)
  
  #Import expression data
  exprs_expt<-as.matrix(read.delim(paste("/Users/admin/Projects/Pathprint/Data/gctfiles/",nid,"_",gpl_platform,".gct", sep = ''), skip = 2, row.names = 1, header = TRUE,  as.is = TRUE))
  exprsfile <- exprs_expt[,-1] #Drop definition column
  #class(exprs_expt)
  #dim(exprs_expt)
  #colnames(exprs_expt)

  threshold <- 0.8
  
  #Fingerprint new expression data
  expt.fingerprint <- exprs2fingerprint(exprs = exprsfile, platform = as.character(gpl_platform), species = "Danio rerio", progressBar = TRUE)
  # Write fingerprint to file
  dir.create(file.path(paste(results.path, nid, "/", sep="")), showWarnings = FALSE)
  write.table(expt.fingerprint, file=paste(results.path, nid, "/scratch/", nid, 
                                      "_pathprint.txt", sep=""),
              sep="\t", quote=FALSE, row.names=TRUE)
  
  # Finding the na subset
  emptydf <- subset(expt.fingerprint, is.na(expt.fingerprint[,1]))
  head(emptydf)
  # Excluding the na data 
  expt.fingerprintnew <- na.omit(expt.fingerprint)

  # Heatmap representation of pathprint
  pdf(file=paste(results.path, nid, "/scratch/", nid,"_heatmap_sdweight.pdf", sep=""))
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
             fontsize_row=10, 
             cluster_cols=TRUE,
             show_rownames=T)
       }
       else{
          write(paste(nid, " for sd ", sd, " has no subset", sep = ""), file = paste(results.path, "logfile.txt", sep=""),
            append = TRUE, sep = "\n")
        }
    }
    else{
      write(paste(nid, " for sd ", sd, " has unlabelled subset\n\t\t", heatmap.data, sep = ""), file = paste(results.path, "logfile.txt", sep=""),
            append = TRUE, sep = "\n")
    }
  }
  dev.off()
}
