# Phylogenetic analysis (Both initial and for meta-analysis of sub-clusters)
# Battleplots, nj tree, max parsimony trees, hclust, heatmap distribution with sd cut-offs
# correlation matrices cast as heatmaps

#pathprint libraries
library(pheatmap)
library("Biobase")
library(pathprint)
library(gplots)
library(ape)
library(phangorn)
library(stringr)
try(library(multicore))
library(plyr)
library(ggplot2)
library(plotrix)
library(cluster)

#Directories & files
celdir = "/Users/admin/Projects/Pathprint/Data/celfiles/"
setwd("/Users/admin/Dropbox/HSPH")
dir.create(results.path <- "results/phylogeny/")
covarsfilename = "metadata/cel_zon.txt" #Zon expts from "Dropbox/HSCI_Curation/Len Zon data/Replicates_DataFiles.xls"
QCresultsDir = "results/ArrayQC/"
pathprintDir = "results/pathprints/"
# To get gpls and id list
gct<-read.delim("metadata/GCT_lenzon.txt") #can set nrows
nidlist<-gct$'Nid'
gpllist<-gct$'GPL'
normtype = "gcrma"
#gfp list = c(7007, 7143, 7211, 7224, 7237, 7341, 7350,7609,7618)
#clust1:  c(7017,7864,7112,7224,7237,7341,7350,7886,7026,7598,7895,7101)
#bls list 7101,7017,7864,7211,7112,7224,7237,7341,7350,7886,7754,7026,7598,7895,7121
#10 somites comparison for gata: 7112,7224,7237,7341,7350,7886,7017, 
#bls new list = 7101,7895,7598,7026,7886,7350,7341,7237,7224,7112,7864,7017
#bls-scl group = 7895, 7026,7864,7017
nidlist = c(7895, 7026,7864,7017) 

type = "pathprint"
rowcount = 633
if(type == "eset"){
  rowcount = 15617
}
expt.pathprints = matrix(, nrow = rowcount, ncol =0) #nrow = 15617 for esets, nrow = 633 for pathprints
setwd("/Users/admin/Dropbox/HSPH")
for(i in 1:length(nidlist)){  #for (i in 1:length(nidlist))
  nid <- nidlist[i]
  gpl_platform = (gct[gct['Nid'] == nid,])$'GPL'[1]
  #import pathprint
  if(file.exists(paste(pathprintDir,"/", nid, "/", nid, "_gcrma_combatparam_", type,".txt", sep=""))){
    pathprint <- as.matrix(read.delim(file=paste(pathprintDir,"/", nid, "/", nid, "_gcrma_combatparam_", type,".txt", sep=""),
                                      row.names = 1, header = TRUE, as.is = TRUE))  
  }else{
    pathprint <- as.matrix(read.delim(file=paste(pathprintDir,"/", nid, "/", nid, 
                        "_cel_",normtype, "_",type,".txt", sep=""), row.names = 1, header = TRUE, as.is = TRUE))
 }
    # Removing arrayQC failed replicates as they come
  if(i>0){
  if(nid == "7026"){pathprint <- pathprint[,-5]}
  if(nid == "7143"){pathprint <- pathprint[,-2]}
  if(nid == "7169"){pathprint <- pathprint[,-3]}
  if(nid == "7259"){pathprint <- pathprint[,-4]}
  if(nid == "7627"){pathprint <- pathprint[,-10]}
  if(nid == "7765"){pathprint <- pathprint[,-1]}
  if(nid == "7864"){pathprint <- pathprint[,-4]}
  if(nid == "7906"){pathprint <- pathprint[,-6]}
  if(nid == "7917"){pathprint <- pathprint[,-2]}
  if(nid == "7928"){pathprint <- pathprint[,-5]}
  if(nid == "8030"){pathprint <- pathprint[,-3]}
  if(nid == "8041"){pathprint <- pathprint[,-3]}
  if(nid == "8061"){pathprint <- pathprint[,-6]}
  if(nid == "8061"){pathprint <- pathprint[,-7]}
  }
  #append nid to col
  colnames(pathprint) <- paste(nid, colnames(pathprint), sep= "_")
  #Concatenate to existing matrix holder
  rownames(expt.pathprints) <- rownames(pathprint)
  expt.pathprints <- cbind(pathprint,expt.pathprints)
} #length(nidlist)
# Finding the na subset
emptydf <- subset(expt.pathprints, is.na(expt.pathprints[,1]))
# Excluding the na data 
expt.pathprints <- na.omit(expt.pathprints)

#Shorten replicate names for ease of visualization
colnames(expt.pathprints) <- str_replace_all(colnames(expt.pathprints), "\\.positive[s]?","+")
colnames(expt.pathprints)<-str_replace_all(colnames(expt.pathprints), "\\.negative[s]?", "-")
colnames(expt.pathprints)<-str_replace_all(colnames(expt.pathprints), "\\.morphant[s]?", "-")
colnames(expt.pathprints)<-str_replace_all(colnames(expt.pathprints), "\\.morpholino[s]?", "-")
colnames(expt.pathprints)<-str_replace_all(colnames(expt.pathprints), "cells\\.", "")
colnames(expt.pathprints)<-str_replace_all(colnames(expt.pathprints), "Rep", "")
colnames(expt.pathprints)<-str_replace_all(colnames(expt.pathprints), "_whole.embryo", "")
colnames(expt.pathprints)<-str_replace_all(colnames(expt.pathprints), "\\.{2}", "\\.")
colnames(expt.pathprints)<-str_replace_all(colnames(expt.pathprints), "\\.{2}", "\\.")
colnames(expt.pathprints)<-str_replace_all(colnames(expt.pathprints), "\\.{2}", "\\.")

#Define cost matrix
CM<-matrix(c(0,1,2,1,0,1,2,1,0), ncol = 3)
dimnames(CM) <- list(c(-1,0,1), c(-1,0,1))

#Clustering analysis + plot output to file
setwd('/Users/admin/Dropbox/HSPH/results/meta_analysis/')
results.path="newresults/"
pdf(file=paste(results.path, "/","heatmap_",type,"_phyl.pdf", sep=""))

#Battleship plot
#distance matrix
dist.mat<-dist(t(expt.pathprints))
labellength = 0.4
battleship.plot(as.matrix(dist.mat), cex.labels = labellength, maxxspan = 0.45, maxyspan = 0.45, col = "red", border = NA)

#Trees
#phydat object constructed
common.dat <- phyDat(t(expt.pathprints), type = "USER", levels = c(-1,0,1))
#Creating neighbour joining tree
drer.NJ.tree <- nj(dist(t(expt.pathprints)))
#Get max parsimony tree
drer.p.tree <- pratchet(common.dat, start = drer.NJ.tree, k = 50, method = "sankoff", cost= CM, trace = 0)
parsimony(drer.NJ.tree, common.dat, method = "sankoff")

#Print trees to nexus format
write.nexus(drer.NJ.tree, file = paste(results.path, "/", "nj_tree.nxs", sep = ""))

#Print trees
plot(drer.NJ.tree, cex = labellength, main = "NJ tree")
#plot(drer.NJ.tree, type = "f", cex = labellength, main = "NJ tree circular")
plot(root(drer.NJ.tree,1), cex = labellength, main = "Rooted NJ tree")
plot(drer.p.tree, cex = labellength, main="Maximum Parsimony tree (no bootstrapping)")

#Plotting hclust
hc <- hclust(dist(t(expt.pathprints)), method = "average")
plot(hc, main = "Default from hclust average", cex = 0.6, hang = 0.1, axes = TRUE, frame.plot = FALSE)

#kmedoids approach
for(i in seq(3,5,1)){
  k= i
  pamm = pam(dist.mat, k)
  clusplot(pamm, main = paste("K-medoids clusplot for k = ", k,sep=""), labels = 4)
  plot(pamm, col = pamm$medoids, main = paste("K-medoids Silhouette for k = ", k, sep = ""))
 # pamm$medoids
}
#heatmaps
expt.fingerprintnew <- expt.pathprints 
if(type == "pathprint"){
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
               fontsize_col = 7, cex = 0.6,
               cluster_cols=TRUE,
               show_rownames=T)
    }}}}else{
  heatmap.data <- expt.fingerprintnew
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
}}}
dev.off()

# Additional phlyogenetic clustering for meta analysis of sub clusters of interest
setwd('/Users/admin/Dropbox/HSPH/results/meta_analysis/')
results.path="newresults/"
normal <- normal_path
combat <- combat_path
pdf(file=paste(results.path, "/","test2_phyl.pdf", sep=""))
plot(hclust(dist(t(normal)),method="average"), main = "Normal pathprints hclust average", cex = 0.4, hang = 0.1, axes = TRUE, frame.plot = FALSE)
plot(hclust(dist(t(combat)),method="average"), main = "Combat pathprints hclust average", cex = 0.4, hang = 0.1, axes = TRUE, frame.plot = FALSE)
#plot(hclust(dist(t(combatonly)),method="average"), main = "Combat only pathprints hclust average", cex = 0.4, hang = 0.1, axes = TRUE, frame.plot = FALSE)

plot(hclust(dist(t(cbind(normal,combat))),method="average"), main = "Combined pathprints hclust average", cex = 0.3, hang = 0.1, axes = TRUE, frame.plot = FALSE)

dev.off()

corel.eset <- (round(cor(x=normal,y=combat,use="everything", method = "pearson")))
plot(im(corel.eset[nrow(corel.eset):1,]), main = "Correlation Matrix Map")
corRaw <- corel.eset
dissimilarity<-1-cor(Raw)
distance<-as.dist(dissimilarity)
plot(hclust(distance),  main="Dissimilarity = 1 - Correlation")
                        