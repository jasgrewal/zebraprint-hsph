#Script used for preparing eset for combat normalization
celdir = "/Users/admin/Projects/Pathprint/Data/celfiles/"
setwd("/Users/admin/Dropbox/HSPH")
#metadata
covarslist <- read.delim(file.path(covarsfilename), row.names = 1, sep = "\t")
gct<-read.delim("metadata/GCT_lenzon.txt") #can set nrows

rm(eset)
nidlist = c(7026,7598,7895)
gpllist<-gct$'GPL'
setwd("/")
for(i in 1:length(nidlist)){ #Iterating through expt id's and fingerprinting each one
  
  #Import data and metadata
  ##Must have phenotype file (covars)
  nid <- nidlist[i]
  gpl_platform = (gct[gct['Nid'] == nid,])$'GPL'[1]
  
  covars <- covarslist[covarslist["exp_nid"] == nid,]
  celFiles <- file.path(celdir, row.names(covars))
  affyRaw <- ReadAffy(filenames = celFiles) # Read in cel files, create affybatch file
  
  pData(affyRaw)<-covars #Assigning phenotype data
  sampleNames(affyRaw)<-pData(affyRaw)$replicate_name #Assigning sample names
  validObject(affyRaw) #check if valid object
  
  #Processing
  ##Raw Data QC (using arrayQualityMetrics)
  normtype = "gcrma" # CHANGE HERE FOR NORMALIZATION STEP TOO (gcrma? rma?)
  ##Background correction and normalization (using RMA/gcRMA)
  #NORMALIZATION - rma or GCRMA
  # IF RMA: 
  #affyNorm <- rma(affyRaw, background = TRUE, normalize = TRUE) 
  affyNorm <- justGCRMA(filenames = celFiles, sampleNames = sampleNames(affyRaw))
  
  esetlast<-exprs(affyNorm)
  #append nid to col
  colnames(esetlast) <- paste(nid, colnames(esetlast), sep= "_")
  if(i==1){
    eset <- esetlast
  }
  else{
  eset<- cbind(esetlast,eset)
  }
}

eset_pheno = matrix(, nrow = ncol(eset), ncol = 3)
colnames(eset_pheno) =c("Array name", "Batch", "Covariate")
rownames(eset_pheno) = colnames(eset)
write.table(eset_pheno, file = "/Users/admin/Dropbox/HSPH/results/meta_analysis/sampleinfo_raw.txt",sep="\t", quote=FALSE, row.names=TRUE)
eset_14som <- eset
eset<-cbind(eset_5som,eset_10som,eset_14som,eset_36som)
write.table(as.matrix(eset), file = "/Users/admin/Dropbox/HSPH/results/meta_analysis/eset.txt",
            sep="\t")


# SAVING COMBAT NORMALIZED ESETS TO RESPECTIVE FOLDERS IN /RESULTS/PATHPRINTS/NID/...

library(futile.matrix)
# Read in eset file with the combat normalized experiments
setwd('/Users/admin/Dropbox/HSPH/results/meta_analysis/')
eset_5som = read.delim('eset_5som.txt',sep = "\t", dec=".", header = TRUE,  check.names=FALSE, comment.char = "")
eset_10som = read.delim('eset_10som.txt',sep = "\t", dec=".", header = TRUE,  check.names=FALSE, comment.char = "")
eset_14som = read.delim('eset_14som.txt',sep = "\t", dec=".", header = TRUE,  check.names=FALSE, comment.char = "")
eset_36som = read.delim('eset_36som.txt',sep = "\t", dec=".", header = TRUE,  check.names=FALSE, comment.char = "")
eset_ONLY<-cbind(eset_5som,eset_10som,eset_14som,eset_36som)

eset_combat_5som = read.delim('eset_combatparam_5som.txt',sep = "\t", dec=".", header = TRUE,  check.names=FALSE, comment.char = "")
#eset_combat_10som = read.delim('eset_combatparam_10som.txt',sep = "\t", dec=".", header = TRUE,  check.names=FALSE, comment.char = "")
eset_combat_14som = read.delim('eset_combatparam_14som.txt',sep = "\t", dec=".", header = TRUE,  check.names=FALSE, comment.char = "")
eset_combat_36som = read.delim('eset_combatparam_36som.txt',sep = "\t", dec=".", header = TRUE,  check.names=FALSE, comment.char = "")
eset_COMBAT<-cbind(eset_combat_5som,eset_combat_14som,eset_combat_36som)
#Set output directory to where pathprints are stored (+ justgcrma normalized esets)
setwd("/Users/admin/Dropbox/HSPH")
results.path <- "results/pathprints/"
#Select subset of normalized esets by column name (experiment)
eset<-eset_COMBAT
nidlist = c(7017,7864,7112,7224,7237,7341,7350,7886,7026,7598,7895,7158,7169,7180,7248,7618,7765,7784,7906,7917,7928,8030,8041,8052)
nidlist<-c(7864,7017,7895,7026,8052,8041,8030,7928,7917,7906,7784,7765,7618,7248,7180,7169,7158)
for(i in 1:length(nidlist)){  
  nid <- nidlist[i]
  eset.subset <- select(eset, row.pat = NULL, col.pat = paste("^",nid, sep=""))
  colnames(eset.subset) <- str_replace_all(colnames(eset.subset), paste("^",nid, sep=""),"")
  colnames(eset.subset) <- str_replace_all(colnames(eset.subset), "^_","")
  #Output combat eset to file
  write.table(eset.subset, file=paste(results.path, nid, "/", nid, 
                                           "_gcrma_combatparam", "_eset.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)
  # Process and write pathway fingerprint to file
  expt.fingerprint <- exprs2fingerprint(exprs = eset.subset, platform = as.character(gpl_platform), species = "Danio rerio", progressBar = TRUE)
  write.table(expt.fingerprint, file=paste(results.path, nid, "/", nid, 
             "_gcrma_combatparam", "_pathprint.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)
}