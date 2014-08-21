#NodePathway for initial clusters
#Get NJ tree and fingerprint matrix
#library
library(gplots)
library(ape)
library(phangorn)
library(stringr)
library(pheatmap)
library(cluster)
library(pathprint)
#Directories & files
setwd("/Users/admin/Dropbox/HSPH")
dir.create(results.path <- "results/phylogeny/branch/")
covarsfilename = "metadata/cel_zon.txt" #Zon expts from "Dropbox/HSCI_Curation/Len Zon data/Replicates_DataFiles.xls"
QCresultsDir = "results/ArrayQC/"
pathprintDir = "results/pathprints/"

#nodePathways function
nodePathways<-function(tree, matrix, threshold){
  # Function to extract the informative pathways from a tree
  # built from pathway fingerprints
  # tree - tree built using ape
  # matrix - fingerprint matrix with colnames equal to tree tip names
  # threshold - threshold at which to determine informative pathways
  # Returns list of
  # 1) informative pathways for each node
  # 2) branch tip names for each node
  library(pathprint)
  library(ape)
  nodeTipList <- lapply(1:tree$Nnode, node2tips, tree)
  pathwayList <- vector("list", length(nodeTipList))
  # split into all possible 2 fold combinations
  for (i in 1:length(nodeTipList)){
    combinations<-combn(nodeTipList[[i]], 2)
    for (j in 1:ncol(combinations)){
      submatrix <- matrix[,unlist(combinations[,j])]
      fac <- c(rep(1, length(unlist(combinations[1,j]))),
               rep(2, length(unlist(combinations[2,j])))
      )
      pathwayList[[i]] <- unique(c(
        pathwayList[[i]],
        diffPathways(submatrix, fac, threshold)
      ))
    }
  }
  names(pathwayList)<-tree$node.label
  names(nodeTipList)<-tree$node.label
  return(list(InformPathways = pathwayList, BranchTipNames = nodeTipList))
}
node2tips<-function(nodeNumber, tree){
  # function to determine all the branch tip names from node numbers
  node.types1 = rep("tip", length(tree$tip.label))
  node.types2 = rep("internal", length(tree$node.label))
  node.types = c(node.types1, node.types2)
  internalNodes <- findall("internal", node.types)
  nodeRef = nodeNumber + length(tree$tip.label)
  searchList<-recursiveSearch(nodeRef, internalNodes, tree)
  return(lapply(searchList, function(x){tree$tip.label[unlist(x)]}))
}
recursiveSearch<-function(num, internalIndex, tr){
  # function to recursively search trees to extract tip names
  daughter_edgenums = findall(num, tr$edge[,1])
  daughter_nodenums = tr$edge[,2][daughter_edgenums]
  if (length(daughter_nodenums) > 0){
    side <- vector('list', length(daughter_nodenums)) 
    for (i in 1:length(daughter_nodenums)){
      side[[i]]<-daughter_nodenums[i]
      if (side[[i]] %in% internalIndex){
        side[[i]]<-recursiveSearch(side[[i]], internalIndex, tr)
      } }
    return(side)
  }
  else {return(num)}
}
findall <- function(what, inlist){
  # find all matches
  TFmatches = inlist == what
  indices = 1:length(inlist)
  matching_indices = indices[TFmatches]
  return(matching_indices)
}

# To get gpls and id list
gct<-read.delim("metadata/GCT_lenzon.txt") #can set nrows
nidlist<-gct$'Nid'
gpllist<-gct$'GPL'
normtype = "gcrma"
#nidlist for scl, bls, cd41 is c(7765, 7906, 7917,7627,7169,8061,7259,7143,7864,7026,7928,8030,8041)
#nidlist for gata with scl is c(8052,7341,7350, 7864,7895,7906)
#nidlist for both bls, cdx4, and scl is c(7765, 7928,7627,7143,7026,7017,7864,7101,7169,8061,7259,7917,7906)
#nidlist for cluster 2 is c(7928,7906,7917,7765,7627,7169,8030,8041)
#nidlist for 3 dpf and greater is 8061,7259,8082,7627,7037,7007,7202,7440,7589
#nidlist for 14 som and 36 hpf diff is c(7765, 7906,7917,7906,7627,8061,7259,7440,7143,7864,7026,7928,8030,8041)
#nidlist for 5-14 soms for gata expt clustering c(7017,7026,7112,7121,7211,7224,7237,7341,7350,7598,7754,7864)
nidlist = c(7017,7864) #uncomment and use for considering subsets of experiments
expt.pathprints = matrix(, nrow = 633, ncol =0)
threshold = 0.8
for(i in 1:length(nidlist)){  #for (i in 1:length(nidlist))
  nid <- nidlist[i]
  gpl_platform = (gct[gct['Nid'] == nid,])$'GPL'[1]
  #import pathprint
  pathprint <- as.matrix(read.delim(file=paste(pathprintDir,"/", nid, "/", nid, 
                                               "_cel_",normtype, "_pathprint.txt", sep=""), row.names = 1, header = TRUE, as.is = TRUE))
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
  #Make consensus for expt & rename column name
  temp.cons <- as.matrix(consensusFingerprint(pathprint, threshold = threshold))
  colnames(temp.cons) <- paste(nid, "consensus", threshold, sep = "_")
  #Concatenate to existing matrix holder
  rownames(temp.cons) <- rownames(temp.cons)
  expt.pathprints <- cbind(temp.cons,expt.pathprints)
} #length(nidlist)
# Finding the na subset
emptydf <- subset(expt.pathprints, is.na(expt.pathprints[,1]))
# Excluding the na data 
expt.pathprints <- na.omit(expt.pathprints)
#Trees
#phydat object constructed
common.dat <- phyDat(t(expt.pathprints), type = "USER", levels = c(-1,0,1))
#Define cost matrix
CM<-matrix(c(0,1,2,1,0,1,2,1,0), ncol = 3)
dimnames(CM) <- list(c(-1,0,1), c(-1,0,1))

#nodePathways analysis
drerio.fingerprint = expt.pathprints

#Creating neighbour joining tree
drerio.nj <- nj(dist(t(drerio.fingerprint)))
drerio.nj$node.label = 1:drerio.nj$Nnode #letters[1:drerio.nj$Nnode]
#use nodePathways function
drerio.node.pathways <- nodePathways(tree = drerio.nj, matrix = drerio.fingerprint[], threshold = .8)
#phydat object constructed
common.dat <- phyDat(t(drerio.fingerprint), type = "USER", levels = c(-1,0,1))
#Get max parsimony tree
drer.p.tree <- pratchet(common.dat, start = drerio.nj, k = 50, method = "sankoff", cost= CM, trace = 0)

#for node a (informative pathways)
pdf(file=paste(results.path, "/","test1.pdf", sep=""))
plot((drerio.nj), show.node.label = TRUE, cex = 0.7)

branchnode1 <- as.matrix(drerio.fingerprint[drerio.node.pathways[["InformPathways"]]$'7',unlist(drerio.node.pathways[["BranchTipNames"]]$'7')])
branchnode2 <- as.matrix(drerio.fingerprint[drerio.node.pathways[["InformPathways"]]$'8',unlist(drerio.node.pathways[["BranchTipNames"]]$'8')])
branchnode3 <- as.matrix(drerio.fingerprint[drerio.node.pathways[["InformPathways"]]$'9',unlist(drerio.node.pathways[["BranchTipNames"]]$'9')])
branchnode4 <- as.matrix(drerio.fingerprint[drerio.node.pathways[["InformPathways"]]$'10',unlist(drerio.node.pathways[["BranchTipNames"]]$'10')])

sd = 0.9
heatmap.data <- branchnode1[apply(branchnode1, 1, sd) > sd,]
pheatmap(heatmap.data, 
         col=c("blue", "white", "red"), 
         legend_breaks=c(-1,0,1), 
         main=paste("All Samples, sd=", sd, sep=""), 
         show_colnames=T, cex = .7,

         cluster_cols=TRUE,
         show_rownames=T)

pheatmap(branchnode1,
        scale = "none", col = c("blue", "white", "red"),
        legend_breaks=c(-1,0,1), 
        #main=paste("All Samples, sd=", sd, sep=""), 
        show_colnames=T, 
        cluster_cols=TRUE,
        show_rownames=T,
        cex = 0.8  )
sd = 0.5
heatmap.data <- branchnode2[apply(branchnode2, 1, sd) > sd,]
pheatmap(heatmap.data, 
         col=c("blue", "white", "red"), 
         legend_breaks=c(-1,0,1), 
         main=paste("All Samples, sd=", sd, sep=""), 
         show_colnames=T, cex = .8,
         
         cluster_cols=TRUE,
         show_rownames=T)
pheatmap(branchnode2,
         scale = "none", col = c("blue", "white", "red"),
         legend_breaks=c(-1,0,1), 
         #main=paste("All Samples, sd=", sd, sep=""), 
         show_colnames=T, 
         cluster_cols=TRUE,
         show_rownames=T,
         cex = 0.9   )

sd = 0.6
heatmap.data <- branchnode3[apply(branchnode3, 1, sd) > sd,]
pheatmap(heatmap.data, 
         col=c("blue", "white", "red"), 
         legend_breaks=c(-1,0,1), 
         main=paste("All Samples, sd=", sd, sep=""), 
         show_colnames=T, cex = .8,
         
         cluster_cols=TRUE,
         show_rownames=T)
pheatmap(branchnode3,
         scale = "none", col = c("blue", "white", "red"),
         legend_breaks=c(-1,0,1), 
         #main=paste("All Samples, sd=", sd, sep=""), 
         show_colnames=T, 
         cluster_cols=TRUE,
         show_rownames=T,
         cex = 0.9   )
sd = 0.6
heatmap.data <- branchnode4[apply(branchnode4, 1, sd) > sd,]
pheatmap(heatmap.data, 
         col=c("blue", "white", "red"), 
         legend_breaks=c(-1,0,1), 
         main=paste("All Samples, sd=", sd, sep=""), 
         show_colnames=T, cex = .8,
         
         cluster_cols=TRUE,
         show_rownames=T)
pheatmap(branchnode4,
         scale = "none", col = c("blue", "white", "red"),
         legend_breaks=c(-1,0,1), 
         #main=paste("All Samples, sd=", sd, sep=""), 
         show_colnames=T, 
         cluster_cols=TRUE,
         show_rownames=T,
         cex = 0.9   )

dev.off()
