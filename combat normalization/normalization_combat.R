#This script reads in an eset, sample info distinguishing batches and covariates, and conducts CoMBat normalization
#Also has notes on pca analysis after combat normalization
library(sva)
library(affy)
library(affyPLM)
library(genefilter)
library(limma)
library(Biobase)
library(gplots)
library(ape)
library(phangorn)
setwd('/Users/admin/Dropbox/HSPH/results/meta_analysis/')
combat.output = "/Users/admin/Dropbox/HSPH/results/meta_analysis/eset_combatparam_36som.txt"
combatpathprints.result = "/Users/admin/Dropbox/HSPH/results/meta_analysis/ComBat normalization analysis/"
eset = read.table('eset_36som.txt',sep = "\t", dec=".", header = TRUE,  check.names=FALSE, comment.char = "")
sample_info = read.delim('sample_info_wt.txt')
edata <- eset[,colnames(eset)%in%(sample_info$Array.name)]
pheno = sample_info
edata <- ExpressionSet(as.matrix(edata))
batch = pheno$Batch
mod = model.matrix(~as.factor(Covariate), data=pheno)
mod0 = model.matrix(~1, data = pheno)

#filter genes with highest variance
#assess variance
vars<-apply(edata,1,var)
#Remove genes with lowest variance
edata.var <- varFilter(edata, var.func = IQR, var.cutoff = 0.1, filterByQuantile = TRUE)
edata.var <- as.matrix(edata.var)
edatanew<-edata.var
#Get rows dropped by varFilter
edata.droppedvar <- as.matrix(edata[setdiff(rownames(edata), rownames(edata.var)),])
#run combat
combatresult <-ComBat(dat = as.matrix(edatanew),batch = batch, mod = mod,numCovs = NULL,par.prior = TRUE)
#Append rows dropped by varFilter
combatresult.varadded <- rbind(edata.droppedvar, combatresult)
write.table(combatresult.varadded, file = combat.output,sep="\t", quote=FALSE, row.names=TRUE)

#Comparing pathprints of original eset, combat normalized, and combat normalized + var drops added
pdf(file=paste(combatpathprints.result, "36som_combat_pathprint.pdf", sep = ""))
#Original pathprint
expt.fingerprint <- exprs2fingerprint(exprs =edatanew, platform = as.character(gpl_platform), species = "Danio rerio", progressBar = TRUE)
# Finding the na subset
emptydf <- subset(expt.fingerprint, is.na(expt.fingerprint[,1]))
# Excluding the na data 
expt.fingerprintnew <- na.omit(expt.fingerprint)
for (sd in seq(0.6,0.8, 0.1)){
heatmap.data <- expt.fingerprintnew[apply(expt.fingerprintnew, 1, sd) > sd,]
pheatmap(as.matrix(heatmap.data), 
         col=c("blue", "white", "red"), 
         legend_breaks=c(-1,0,1), 
         main=paste("Non-combat normalized pathprint, sd=", sd, sep=""),
         show_colnames=T, 
         cex = 0.6,
         cluster_cols=TRUE,
         show_rownames=T)
}
#Combat fingerprint
expt.fingerprint <- exprs2fingerprint(exprs =combatresult, platform = as.character(gpl_platform), species = "Danio rerio", progressBar = TRUE)
# Finding the na subset
emptydf <- subset(expt.fingerprint, is.na(expt.fingerprint[,1]))
# Excluding the na data 
expt.fingerprintnew <- na.omit(expt.fingerprint)
for (sd in seq(0.6,0.8, 0.1)){
  heatmap.data <- expt.fingerprintnew[apply(expt.fingerprintnew, 1, sd) > sd,]
  pheatmap(as.matrix(heatmap.data), 
         col=c("blue", "white", "red"), 
         legend_breaks=c(-1,0,1), 
         main=paste("Combat pathprint, sd=", sd, sep=""), 
         show_colnames=T, 
         cex=0.6,
         cluster_cols=TRUE,
         show_rownames=T)
}
#Combat fingerprint with complete list of genes
expt.fingerprint <- exprs2fingerprint(exprs =combatresult.varadded, platform = as.character(gpl_platform), species = "Danio rerio", progressBar = TRUE)
# Finding the na subset
emptydf <- subset(expt.fingerprint, is.na(expt.fingerprint[,1]))
# Excluding the na data 
expt.fingerprintnew <- na.omit(expt.fingerprint)
for (sd in seq(0.6,0.8, 0.1)){
  heatmap.data <- expt.fingerprintnew[apply(expt.fingerprintnew, 1, sd) > sd,]
  pheatmap(as.matrix(heatmap.data), 
         col=c("blue", "white", "red"), 
         legend_breaks=c(-1,0,1), 
         main=paste("Combat + var filter genes pathprint, sd=", sd, sep=""), 
         show_colnames=T, 
         cex=0.6,
         cluster_cols=TRUE,
         show_rownames=T)
}
dev.off()

pdf(file=paste(combatpathprints.result, "36som_combat_compare.pdf", sep=""))

#Phylogenetics (just nj tree based on eigen distance)
#phydat object constructed

exptset <- as.matrix(eset_36som)
drer.NJ.tree <- nj(dist(t(exptset)))
#Get max parsimony tree
#Print trees
#plot(drer.NJ.tree, cex = labellength, main = "NJ tree non comBat normalized")
plot(root(drer.NJ.tree,1), cex = labellength, main = "Rooted NJ tree non comBat normalized")
common.dat <- phyDat(t(expt.pathprints), type = "USER", levels = c(-1,0,1))
drer.p.tree <- pratchet(common.dat, start = drer.NJ.tree, k = 50, method = "sankoff", cost= CM, trace = 0)         
plot(drer.p.tree, cex = labellength, main="Maximum Parsimony tree (no bootstrapping)")

pheatmap(as.matrix(exptset), 
         col=c("blue", "white", "red"), 
         legend_breaks=c(-1,0,1), 
         main="All Samples non Combat, sd=0.8", 
         show_colnames=T, 
         fontsize_row=6, 
         fontsize_col = 7,
         cluster_cols=TRUE,
         show_rownames=T)
#Creating neighbour joining tree
exptset <- as.matrix(eset_combat_36som)
drer.NJ.tree <- nj(dist(t(exptset)))
#Get max parsimony tree
#Print trees
#plot(drer.NJ.tree, cex = labellength, main = "NJ tree comBat normalized")
plot(root(drer.NJ.tree,1), cex = labellength, main = "Rooted NJ tree comBat normalized")
# HEATMAP with sds. Derived from John (github) and Shannan's codes.
pheatmap(as.matrix(exptset), 
               col=c("blue", "white", "red"), 
               legend_breaks=c(-1,0,1), 
               main="All Samples Combat, sd=0.8", 
               show_colnames=T, 
               fontsize_row=6, 
               fontsize_col = 7,
               cluster_cols=TRUE,
               show_rownames=T)
#Creating neighbour joining tree
exptset <- as.matrix(eset_combatonly_36som)
drer.NJ.tree <- nj(dist(t(exptset)))
#Get max parsimony tree
#Print trees
#plot(drer.NJ.tree, cex = labellength, main = "NJ tree comBat normalized + var set")
plot(root(drer.NJ.tree,1), cex = labellength, main = "Rooted NJ tree comBat normalized + var set")
# HEATMAP with sds. Derived from John (github) and Shannan's codes.
pheatmap(as.matrix(exptset), 
         col=c("blue", "white", "red"), 
         legend_breaks=c(-1,0,1), 
         main="All Samples Combat complete, sd=0.8", 
         show_colnames=T, 
         fontsize_row=6, 
         fontsize_col = 7,
         cluster_cols=TRUE,
         show_rownames=T)
dev.off()

#PCA analysis
#performs pca after scaling the data
#Use prcomp for numerical stability. Can center and scale if magnitude of the numbers in the matrix are not of a comparable size
pca <- prcomp(t(combatresult), scale=F, center = F) 
pca$loadings
summary(pca)
pcs<-data.frame(pca$x)
#add in the annotation data by merging for the 1st 2 principal components
pcs<-merge(pcs[,1:2],pheno, by = 0)
#plot first two principal components
library(lattice)
#plot by sample - selecting which pc's to keep
#pick where greatest change occurs to screeplot
screeplot(pca, type = "lines")
#xyplot(PC2~PC1,pcs,group=sample,auto.key=T)
#or use Kaiser's criterion to retain pc's with variance above 1 only (applied pca to standardized data)
(pca$sdev)^2
#or keep pc's that account for total 80% (benchmark) of variance
summary(pca)
library(car)
scatterplotMatrix(combatresult[,4:7])
#2) matrix of eigenvectors (rotation)
#3) principal component data (x)
#4) centering (center)
#5) scaling (scale)
summary(pca)$importance[,1:5]

#BIOCONDUCTOR VIGNETTE
x11(height=6, width=12, pointsize=12); par(mfrow=c(1,2)) # Set plotting parameters.
mycolors <- c("red", "green", "blue", "magenta", "black") # Define plotting colors.
plot(pca$x, pch=20, col=mycolors[sort(rep(1:5, 500))]) 
# Plots scatter plot for the first two principal components that are stored in pca$x[,1:2].
plot(pca$x, type="n"); text(pca$x, rownames(pca$x), cex=0.8, col=mycolors[sort(rep(1:5, 500))])
# Same as above, but prints labels.
library(geneplotter); smoothScatter(pca$x) # Same as above, but generates a smooth scatter plot that shows the density of the data points.
pairs(pca$x[,1:8], pch=20, col=mycolors[sort(rep(1:5, 500))]) 
# Plots scatter plots for all combinations between the first four principal components.
biplot(pca, cex = 0.5) 

