zebraprint
==========
This project aimed to analyze microarray gene expression data at a pathway level, using PathPrint (http://compbio.sph.harvard.edu/hidelab/pathprint/Pathprint.html)
Starting with raw cel files (or alternatively, gct files), justGCRMA normalization was done and the expression sets pathprinted to get the pathway signatures.
Initial clustering analysis (with battleplots, nj clustering, max-parsimony trees, and k-medoid clustering) revealed some interesting clustering. We took a look at what pathways are driving the node separation for these subclusters.
We also did some combat normalization to remove batch effects confounding wildtype signatures, and subsequent meta-analysis gave us greater resolution into changes between samples at the pathway and gene-expression levels. 
The scripts have been organized by the steps in the analysis. Please refer to the Folders Description section for a brief summary of the scripts in each folder and where the relevant results have been output.

==========
Folder Descriptions
==========
1) Preprocessing and pathprinting: Preprocessing and normalization of raw data files (cel files or gct files). Subsequently pathprinted the gene expression sets (results in results/pathprints posted on
basecamp/[hidelab]/zebraprint)
2) Initial analysis: Initial clustering analysis of experiments to detect any strong sub-clusters. Preliminary branch-based pathways extraction of some sub-clusters of interest using nodepathway (Gabriel
Altschuler’s code). (results in results/phylogeny posted on basecamp/[hidelab]/zebraprint).
3) combat normalization: There were some batch effects confounding wild-type signals for experiments conducted at the same somite stage. ComBat normalization was conducted (parametric as well as non-parametric) to get a better resolution of the mutant vs wildtype effect in the experiments of interest. (results in /results/meta_analysis/ComBat normalization analysis).
4) nodepath_geneexprn_BLS.R : meta-analysis conducted to understand what is driving the difference between samples of interest, at the pathways (analyzing the nodes), and at the gene (differential expression) level. (results in /results/meta_analysis/newresults/bls/)
