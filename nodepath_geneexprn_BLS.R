#NodePathways for bloodless and scl mutant experiments (meta-analysis)
#Get NJ tree and fingerprint matrix
#Also contains script for extraction of significant genes using limma/samr 
#library
library(gplots)
library(ape)
library(phangorn)
library(stringr)
library(pheatmap)
library(cluster)
library(pathprint)
library(ape)
#Directories & files
setwd("/Users/admin/Dropbox/HSPH")
#dir.create(results.path <- "results/phylogeny/branch/")
dir.create(results.path <- "results/meta_analysis/newresults/bls/branch/")
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
nidlist = c(7895, 7026,7864,7017,7101) # bls and scl experiments, with one outgroup experiment (7101)
nidlist = c(7864,7017)#,7026)#,7895)
threshold = 0.8
type = "pathprint"
rowcount = 633
if(type == "eset"){
  rowcount = 15617
}
# Importing expression set/pathprint (depending on type)
expt.pathprints = matrix(, nrow = rowcount, ncol =0)
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
  # Removing arrayQC failed replicates
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
  #Concatenate to existing matrix holder
  colnames(pathprint) <- paste(nid, colnames(pathprint), sep= "_")
  rownames(expt.pathprints) <- rownames(pathprint)
  expt.pathprints <- cbind(pathprint,expt.pathprints)
} #length(nidlist)
# Finding the na subset
emptydf <- subset(expt.pathprints, is.na(expt.pathprints[,1]))
# Excluding the na data 
expt.pathprints <- na.omit(expt.pathprints)

#PART 1: EXTRACT PATHWAYS IMPORTANT FOR BRANCHES (FROM SOURCE: PATHPRINTS)
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
tree = root(drerio.nj,1)
drerio.node.pathways <- nodePathways(tree = tree, matrix = drerio.fingerprint[], threshold = .8)
#for nodes (informative pathways)
png(file=paste(results.path,"/",type, "_njlabel.png",sep=""), res=180, height=900, width=1600) 
plot(tree, show.node.label = TRUE, cex = 0.6, main = paste(type, " Rooted NJ tree", sep = ""))
dev.off()
pdf(file=paste(results.path, "/",type,"eset_allheatmaps_test1.pdf", sep=""))
#Branchnodes decided after looking at the nj tree plot
branchnode1 <- as.matrix(drerio.fingerprint[drerio.node.pathways[["InformPathways"]]$'1',unlist(drerio.node.pathways[["BranchTipNames"]]$'1')])
branchnode2 <- as.matrix(drerio.fingerprint[drerio.node.pathways[["InformPathways"]]$'7',unlist(drerio.node.pathways[["BranchTipNames"]]$'7')])
branchnode3 <- as.matrix(drerio.fingerprint[drerio.node.pathways[["InformPathways"]]$'12',unlist(drerio.node.pathways[["BranchTipNames"]]$'12')])
branchnode4 <- as.matrix(drerio.fingerprint[drerio.node.pathways[["InformPathways"]]$'13',unlist(drerio.node.pathways[["BranchTipNames"]]$'13')])
branchnode5 <- as.matrix(drerio.fingerprint[drerio.node.pathways[["InformPathways"]]$'11',unlist(drerio.node.pathways[["BranchTipNames"]]$'11')])

for(sd in seq(0,0.8,0.2)){
heatmap.data <- branchnode1[apply(branchnode1, 1, sd) > sd,]
write.table(heatmap.data, file = paste(results.path,type,"_bls5_param_",sd,".txt", sep=""),sep="\t", quote=FALSE, row.names=TRUE)
pheatmap(heatmap.data, 
         col=c("blue", "white", "red"), 
         legend_breaks=c(-1,0,1), 
         main=paste("All Samples, sd=", sd, sep=""), 
         show_colnames=T, cex = .9,

         cluster_cols=TRUE,
         show_rownames=T)
}
#for(sd in seq(0,.8,.2)){
heatmap.data <- branchnode2[apply(branchnode2, 1, sd) > sd,]
write.table(heatmap.data, file = paste(results.path, "/branchlist/",type,"_0sd_branch2.xls", sep=""),sep="\t", quote=FALSE, row.names=TRUE)
pheatmap(heatmap.data, 
         col=c("blue", "white", "red"), 
         legend_breaks=c(-1,0,1), 
         main=paste("All Samples, sd=", sd, sep=""), 
         show_colnames=T, cex = .8,
         cluster_cols=TRUE,
         show_rownames=T)
}
#for(sd in seq(0,.8,.2)){
heatmap.data <- branchnode3[apply(branchnode3, 1, sd) > sd,]
write.table(heatmap.data, file = paste(results.path, "/branchlist/",type,"_0sd_branch3.xls", sep=""),sep="\t", quote=FALSE, row.names=TRUE)
pheatmap(heatmap.data, 
         col=c("blue", "white", "red"), 
         legend_breaks=c(-1,0,1), 
         main=paste("All Samples, sd=", sd, sep=""), 
         show_colnames=T, cex = .8,
         cluster_cols=TRUE,
         show_rownames=T)
}
#for(sd in seq(0,.8,.2)){
heatmap.data <- branchnode4[apply(branchnode4, 1, sd) > sd,]
write.table(heatmap.data, file = paste(results.path, "/branchlist/",type,"_0sd_branch4.xls", sep=""),sep="\t", quote=FALSE, row.names=TRUE)
pheatmap(heatmap.data, 
         col=c("blue", "white", "red"), 
         legend_breaks=c(-1,0,1), 
         main=paste("All Samples, sd=", sd, sep=""), 
         show_colnames=T, cex = .8,
         cluster_cols=TRUE,
         show_rownames=T)
}
sd=.9
heatmap.data <- branchnode5[apply(branchnode5, 1, sd) > sd,]
write.table(heatmap.data, file = paste(results.path, "/branchlist/",type,"_0sd_branch5.xls", sep=""),sep="\t", quote=FALSE, row.names=TRUE)
pheatmap(heatmap.data, 
         col=c("blue", "white", "red"), 
         legend_breaks=c(-1,0,1), 
         main=paste("All Samples, sd=", sd, sep=""), 
         show_colnames=T, cex = .8,        
         cluster_cols=TRUE,
         show_rownames=T)
dev.off()


#PART 2: EXTRACT GENE LISTS FOR PATHWAY ENRICHMENT FROM EXPRESSION SETS
###Get differentially expressed genes from expressionset (eset) - using samr
#x=expression set with two groups - branchpoint being studied
#y = list = c(1,1,1,2,2,2,2,2) (for 1s and 2s distinguishing the two branches' replicates)
eset_orig = eset
eset = eset_scl14som # subset from original eset with branch cluster of interest
y = y_scl14som #Branch is mutants of interest (=2), versus rest (controls, = 1)
data=list(x=eset, y=y, genenames = rownames(eset))
data$logged2 = TRUE #sometimes get error logged2 is invalid arg type, so have to run this
samr.obj <- samr(data,resp.type = "Two class unpaired", nperms = 100)
delta.table <- samr.compute.delta.table(samr.obj, min.foldchange=1,nvals=200) #Minimum fold change is 1
siggenes.table <- samr.compute.siggenes.table(samr.obj, del=0.5, data, delta.table,all.genes=TRUE) #delta = 0.5
upgenes<-siggenes.table$genes.up
upgenes <- (upgenes[as.numeric(upgenes[,8])<5,]) #Subselect list of genes with q-value (%) < 5
downgenes<-siggenes.table$genes.lo
downgenes <- (downgenes[as.numeric(downgenes[,8])<5,])  #Subselect list of genes with q-value (%) < 5
samr.dir = "results/meta_analysis/newresults/bls/branch/samr/"
#write to file
write.csv(upgenes,paste(samr.dir, "/DEGS_14scl_up_samr.csv", sep=""))
write.csv(downgenes,paste(samr.dir, "/DEGS_14scl_down_samr.csv", sep=""))

###Get differentially expressed genes from expressionset (eset) - using limma
#eset_14som; eset_5som;
library(limma)
#scl wt-mut and bls wt-mut comparisons
eset<-eset_14som[,6:11] #scl eset
#design group means parametrization
design <- cbind(WT=c(1,1,1,0,0,0),MU=c(0,0,0,1,1,1))
colnames(design) <- c("wt_scl14som", "mut_scl14som")
fit<-lmFit(eset,design)
#make all pariwise comparisons between the three groups, make appropriate contrast matrix
contrast.matrix<-makeContrasts(MUvsWT=mut_scl14som-wt_scl14som,levels=design)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
#list of top genes differentially expressed in group2 versus group1
topTable(fit2, coef=1, adjust = "BH")
#assign outcome of each hypothesis test
results<-decideTests(fit2)
vennDiagram(results, cex = 0.5)

#somite stage comparison
eset_14som<-eset[,16:26]
eset_5som<-eset[,7:15]
limma.dir = "/Users/admin/Dropbox/HSPH/results/meta_analysis/newresults/bls/branch/limma/"
eset<-eset_14som
design <- model.matrix(~ 0+factor(c(2,2,2,1,1,0,0,0,3,3,3)))
colnames(design) <- c("wt_scl14som", "wt_bls14som", "blsmut14som", "sclmut14som")
#-****-
eset<-eset_5som
design <- model.matrix(~ 0+factor(c(2,2,1,1,3,3,3,0,0)))
colnames(design) <- c("wt_scl5som", "wt_bls5som", "blsmut5som", "sclmut5som")

fit<-lmFit(eset,design)
#make all pairwise comparisons between the three groups, make appropriate contrast matrix
contrast.matrix<-makeContrasts(blsmut14som-wt_bls14som, sclmut14som-wt_scl14som, levels = design)
contrast.matrix<-makeContrasts(blsmut5som-wt_bls5som, sclmut5som-wt_scl5som, levels = design)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
#list of top genes differentially expressed in group2 versus group1
#topTable(fit2, coef=1, adjust = "BH")
#assign outcome of each hypothesis test
results<-decideTests(fit2)
vennDiagram(results, cex = 0.9,include=c("up","down"))
dev.print(png, paste(limma.dir,"rplot01.png",sep=""),width=900,pointsize = 12)
vennDiagram(results, cex = 0.9)
dev.print(png, paste(limma.dir,"rplot02.png",sep=""),width=900,pointsize = 12)

dim(results[results[,1]!=0 & results[,2]!=0, ])
#print results to file
write.table(results[results[,1]!=0 | results[,2]!=0, ], paste(limma.dir,"resultlist.txt",sep=""), quote = FALSE, sep = "\t", row.names = TRUE)
#print names of intersecting genes
write.table(rownames(results[results[,1]!=0 & results[,2]!=0, ]), paste(limma.dir,"intersect_genelist.txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rownames(results[results[,1]==1 & results[,2]==1, ]), paste(limma.dir,"intersect_up.txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rownames(results[results[,1]==-1 & results[,2]==-1, ]), paste(limma.dir,"intersect_down.txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE)

#print gene names from each of the two sets
write.table(rownames(results[results[,1]!=0, ]), paste(limma.dir,"1genelist.txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rownames(results[results[,1]==1, ]), paste(limma.dir,"1upgenelist.txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rownames(results[results[,1]==-1, ]), paste(limma.dir,"1downgenelist.txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE)

write.table(rownames(results[results[,2]!=0, ]), paste(limma.dir,"2genelist.txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rownames(results[results[,2]==1, ]), paste(limma.dir,"2upgenelist.txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(rownames(results[results[,2]==-1, ]), paste(limma.dir,"2downgenelist.txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE)


#fit2$F and fit2$F.p combine the 3 pariwise comparisons into one F-test. 
#To find genes which vary b/w the 3 targets in any way, look for genes with small p values
#Find top 30 genes
topTableF(fit2, number = 30)

sel.diif<-p.adjust(fit2$F.p.value,method="fdr")<0.05
deg<-eset[,1:6][sel.diif,]
resultsLimmaCommon <- limma:::topTable(fit2, coef = 1, adjust.method = "fdr", number = nrow(eset),sort.by="P")
resultsLimmaCommon_Sig <- resultsLimmaCommon[abs(resultsLimmaCommon$logFC)> 2 & resultsLimmaCommon$adj.P.Val < 0.05, ]

summary(resultsLimmaCommon_Sig)


#Checking pathprint repository for gene matches to pathways
pathprint.Dr.gs[grep("799805", pathprint.Dr.gs[], ignore.case = TRUE)]
grep("799805", pathprint.Dr.gs[], ignore.case = TRUE)
