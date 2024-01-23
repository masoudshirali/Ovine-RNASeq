
# Load the libraries

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
allowWGCNAThreads()          # allow multi-threading (optional)

# 1. Read the gene counts table 
data=read.csv("5.featurecounts/Lambs.featurecounts.hisat2.Rmatrix",header=T,row.names=1,sep="\t", check.names = FALSE)
data=data[ , !names(data) %in% c("7085","7073")]
colnames(data)<-gsub("Control","",colnames(data))
colnames(data)<-gsub("Low","",colnames(data))
colnames(data)<-gsub("Medium","",colnames(data))
colnames(data)<-gsub("High","",colnames(data))

# Read the metadata
sample_metadata = read.csv(file = "metadata_with_methaneinfoadded_metadata.csv")

# 2. QC - outlier detection
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
pdf("6.deseq2/1.hclust_samples.pdf")
htree <- hclust(dist(t(data)), method = "average")
plot(htree)
dev.off()

# pca - method 2

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)
pdf("6.deseq2/2.pca.pdf")
ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
dev.off()

# exclude outlier samples
#samples.to.be.excluded <- c('GSM4615000', 'GSM4614993', 'GSM4614995')
#data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]

# 3. Normalization
# create a deseq2 dataset

# fixing rownames and colnames in data and sample_metadat
rownames(sample_metadata) <- sample_metadata$ID
sample_metadata$ID <- factor(sample_metadata$ID)
rownames(sample_metadata)<-gsub("[a-zA-Z ]", "", rownames(sample_metadata))
head(sample_metadata)

# making the rownames and column names identical
# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
data <- data[,unique(rownames(sample_metadata))]
all(colnames(data) == rownames(sample_metadata))

# create dds
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = sample_metadata,
                              design = ~ 1) # not spcifying model

## remove all genes with counts < 15 in more than 75% of samples (22*0.75=16)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 10) >= 16,]
nrow(dds75) # 16583 genes


# perform variance stabilization
dds_norm <- vst(dds75)
write.csv(assay(dds_norm),"6.deseq2/Lambs_allSamples_normalized_Counts",row.names=T)


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

# 4. Network Construction 
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                  powerVector = power,
                  networkType = "signed",
                  verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power
pdf("6.deseq2/3.power_threshold.pdf")
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()
  

grid.arrange(a1, a2, nrow = 2)
dev.off()

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

#1.
soft_power <- 14
temp_cor <- cor
cor <- WGCNA::cor

# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                 maxBlockSize = 14000,
                 TOMType = "signed",
                 power = soft_power,
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE,
                 randomSeed = 1234,
                 verbose = 3)


cor <- temp_cor

# 5. Module Eigengenes 
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# grey module = all genes that doesn't fall into other modules were assigned to the grey module

# 6A. Relate modules to traits --------------------------------------------------
# module trait associations



# create traits file - binarize categorical variables
traits <- colData %>% 
  mutate(disease_state_bin = ifelse(grepl('COVID', disease_state), 1, 0)) %>% 
  select(8)


# binarize categorical variables

colData$severity <- factor(colData$severity, levels = c("Healthy", "Convalescent", "ICU", "Moderate", "Severe"))

severity.out <- binarizeCategoricalColumns(colData$severity,
                           includePairwise = FALSE,
                           includeLevelVsAll = TRUE,
                           minCount = 1)


traits <- cbind(traits, severity.out)


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)



# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[18:22],
             y = names(heatmap.data)[1:17],
             col = c("blue1", "skyblue", "white", "pink", "red"))



module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()



# 6B. Intramodular analysis: Identifying driver genes 

# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


module.membership.measure.pvals[1:10,1:10]


# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(norm.counts, traits$data.Severe.vs.all, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)

# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.







dataExpr_deseq <- DESeqDataSetFromMatrix(countData = data0,colData = sample_metadata,design = ~ CH4production)
mcols(dataExpr_deseq)$basepairs = data0$geneLengt1
fpkm_matrix = fpm(dataExpr_deseq)
datExpr = t(log2(fpkm_matrix+1))

match(sample_metadata$ID, colnames(data0))
# Calculate sample distance and cluster the samples
sampleTree = hclust(dist(datExpr), method = "average");
# plot sample tree
pdf(file = "6.deseq2/1-n-sampleClustering.pdf", width = 40, height = 9);
par(cex = 1.3);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

#===============================================================================
#
#  Choose soft threshold parameter
#
#===============================================================================

# Choose a set of soft threshold parameters
powers = c(c(1:20), seq(from = 22, to=30, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) 
# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file = "6.deseq2/2-n-sft.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#===============================================================================
#
#  Turn data expression into topological overlap matrix
#
#===============================================================================

power=sft$powerEstimate #4

# Option 1: automatic
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = power,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3)
cor<- stats::cor

# unsigned -> nodes with positive & negative correlation are treated equally 
# signed -> nodes with negative correlation are considered *unconnected*, treated as zero

sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
pdf(file = "6.deseq2/4-module_tree_blockwise.pdf", width = 8, height = 6);
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
     
################################################################################
################################################################################
# Option 2a: step-by-step
power = power
adjacency = adjacency(datExpr, power = power)
TOM = TOMsimilarity(adjacency); # Turn adjacency into topological overlap
dissTOM = 1-TOM


# Option 2b: 
TOM = TOMsimilarityFromExpr(datExpr, power = power)
dissTOM = 1-TOM 
dim(dissTOM)


#===============================================================================
#
#  Construct modules (proceed with the genetree from option 2b)
#
#===============================================================================
# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file = "6.deseq2/3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

# Module identification using dynamic tree cut
# We like large modules, so we set the minimum module size relatively high:
# minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                            pamRespectsDendro = FALSE,minClusterSize = 30);
table(dynamicMods)
length(table(dynamicMods)) 
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file = "6.deseq2/4-module_tree.pdf", width = 8, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()

#===============================================================================
#
#  Merge modules
#
#===============================================================================
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf("6.deseq2/eigengene_cluster.pdf")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()

# Merge close modules
MEDissThres=0.40
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors = merge$colors  
mergedMEs = merge$newMEs  
# Plot merged module tree
pdf(file = "6.deseq2/5-merged_Module_Tree.pdf", width = 12, height = 9)  
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()
write.table(merge$oldMEs,file="6.deseq2/oldMEs.txt");
write.table(merge$newMEs,file="6.deseq2/newMEs.txt");

#===============================================================================
#
#  Export of networks to external software
#
#===============================================================================

# Export the gene list of old modules 
for (i in 1:length(merge$oldMEs)){
  modules = c(substring(names(merge$oldMEs)[i], 3));
  genes = colnames(datExpr)
  inModule = is.finite(match(dynamicColors,modules))
  modGenes = genes[inModule]
  modTOM=TOM[inModule,inModule]
  dimnames(modTOM)=list(modGenes,modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("output_for_cytoscape/orign_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("output_for_cytoscape/orign_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = dynamicColors[inModule]);
}
# Export the gene list of new modules 
for (i in 1:length(merge$newMEs)){
  modules = c(substring(names(merge$newMEs)[i], 3));
  genes = colnames(datExpr)
  inModule = is.finite(match(dynamicColors,modules))
  modGenes = genes[inModule]
  modTOM=TOM[inModule,inModule]
  dimnames(modTOM)=list(modGenes,modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("output_for_cytoscape/merge_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("output_for_cytoscape/merge_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = dynamicColors[inModule]);
}

#===============================================================================
#
#  PART 1: Correlate module eigen-genes and samples (or other discrete data)
#
#===============================================================================
# Heatmap of old module eigen-genes and samples
pdf(file="6.deseq2/oldMEs.pdf",heigh=80,width=20)
library("pheatmap")
rownames(merge$oldMEs)=names(data0[,-181])
pheatmap(merge$oldMEs,cluster_col=T,cluster_row=T,show_rownames=F,show_colnames=T,fontsize=6)
dev.off()


# Heatmap of new module eigen-genes and sample trait (e.g. Zone)
col_ann <- sample_metadata[,c(1,3)]
rownames(col_ann) <- col_ann[,1]
col_ann <- data.frame(col_ann)
col_ann$Zone <- as.factor(col_ann$Zone)
col_ann <- col_ann[order(col_ann$Zone),]
col_ann$sample_ID <- NULL
head(col_ann)
ann_color <- list("col_ann" = c("Z1" = "yellow",
                                "Z2" = "red",
                                "Z3" = "green"))

data <- data.frame(merge$newMEs)
data <- data[order(match(rownames(data), rownames(col_ann))),]
dim(merge$newMEs)

pdf(file="6.deseq2/newMEs.pdf",heigh=60,width=20)
rownames(merge$newMEs)=names(data0[,-181])
pheatmap(data,cluster_col=T,cluster_row=F,show_rownames=F,
         show_colnames=T,fontsize=6,
         annotation_row = col_ann, annotation_colors = ann_color)
dev.off()

#=====================================================================================
#
#  PART 2: Correlation between gene modules and methane production (continuous data)
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs0)

# Read microbial data as traits
bac_traits = read.table("b_order_234.txt", header = T, sep = "\t")
rownames(bac_traits) = bac_traits[, 1]
bac_traits = bac_traits[, -1]
# sample names should be consistent in eigen genes and traits !!!!
bac_traits = bac_traits[match(rownames(MEs), rownames(bac_traits)), ]
table(rownames(MEs) == rownames(bac_traits))

# Calculate pearson correlation coefficients between module eigen-genes and traits
moduleTraitCor = cor(MEs, bac_traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#write.table(moduleTraitCor,file="moduleTrait_correlation.txt");
#write.table(moduleTraitPvalue,file="moduleTrait_pValue.txt");


#  Plot heatmap of module-traits relationship
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("5-module-traits-bacteria-order.pdf", width = 80, height = 15)
par(mar = c(15, 12, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(bac_traits),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#  Plot heatmap of module-traits relationship
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor[1:25,1:25], 2), "\n(",
                    signif(moduleTraitPvalue[1:25,1:25], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[1:25,1:25])
pdf("5-module-traits-bacteria-order1.pdf", width = 20, height = 10)
par(mar = c(15, 12, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor[1:25,1:25],
               xLabels = colnames(bac_traits[1:25,1:25]),
               yLabels = colnames(MEs[1:25,1:25]),
               ySymbols = colnames(MEs[1:25,1:25]),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#=====================================================================================
#
#   Intramodular analysis: identifying genes with high geneModuleMembership & geneTraitSignficance
#
#=====================================================================================


# Define variable Verru containing the Verrucomicrobiales column of bac_traits
Verru = as.data.frame(bac_traits$Verrucomicrobiales);
names(Verru) = "Verrucomicrobiales"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

MET = orderMEs(cbind(MEs, Verru))

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, Verru, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Verru), sep="");
names(GSPvalue) = paste("p.GS.", names(Verru), sep="");

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(10,10,1,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


module = "lightgreen"
# Rename to moduleColors
moduleColors = mergedColors
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Verrucomicrobiales",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "blue")


## Draw bubble plot for particular module
colsum_bac_traits <- colSums(bac_traits)
colsum_bac_traits <- data.frame(colsum_bac_traits)
colsum_bac_traits$b_order <- rownames(colsum_bac_traits)
library(tidyr)
moduleTraitCor_long <- data.frame(moduleTraitCor)
moduleTraitCor_long$module <- rownames(moduleTraitCor)
moduleTraitCor_long <- moduleTraitCor_long[,c(235,1:234)]
moduleTraitCor_long <- gather(moduleTraitCor_long, b_order, PCC, Pseudomonadales:Others, factor_key = TRUE)

moduleTraitPvalue_long <- data.frame(moduleTraitPvalue)
moduleTraitPvalue_long$module <- rownames(moduleTraitPvalue)
moduleTraitPvalue_long <- moduleTraitPvalue_long[,c(235,1:234)]
moduleTraitPvalue_long <- gather(moduleTraitPvalue_long, b_order, pval, Pseudomonadales:Others, factor_key = TRUE)

moduleTrait_long <- merge(moduleTraitCor_long, moduleTraitPvalue_long, by = c("module","b_order"))

bubble_Data <- merge(moduleTrait_long, colsum_bac_traits, by = "b_order")
#just want module = "lightgreen"
bubble_Data_lightgreen <- bubble_Data[which(bubble_Data$module == "MElightgreen"),]

library(ggplot2)
ggplot(bubble_Data_lightgreen, aes(x= colsum_bac_traits, y= PCC, size = colsum_bac_traits,
                                   color = PCC, label = b_order)) +
  geom_text(hjust = 1, size=3) +
  geom_point(alpha=1) + ylab("Module-taxon correlation") + xlab("Relative abundance (sum)") +
  theme_bw()

 
############# Summary ###################################

head(datExpr)[1:5,1:5] # transcriptome data

head(sample_metadata)[1:5,] # metadata (sample info)
head(bac_traits)[1:5,1:5] # external trait






#=====================================================================================
#
#   Cytoscape
#
#=====================================================================================


#if(!"RCy3" %in% installed.packages()){
#  install.packages("BiocManager")
#  BiocManager::install("RCy3")
#}

# https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/
library(RCy3)

cytoscapePing () # make sure cytoscape is open
cytoscapeVersionInfo ()

###### for yellow module of the merged data (newMEs) #################################
edge <- read.delim("output_for_cytoscape/merge_CytoscapeInput-edges-lightgreen.txt")
colnames(edge)
colnames(edge) <- c("source", "target","weight","direction","fromAltName","toAltName")

node <- read.delim("output_for_cytoscape/merge_CytoscapeInput-nodes-lightgreen.txt")
colnames(node)  
colnames(node) <- c("id","altName","node_attributes") 

createNetworkFromDataFrames(node,edge[1:50,], title="my first network", collection="DataFrame Example")

################ customise the network visualization ##################################
# use other pre-set visual style
setVisualStyle('Marquee')

# set up my own style
style.name = "myStyle"
defaults <- list(NODE_SHAPE="diamond",
                 NODE_SIZE=30,
                 EDGE_TRANSPARENCY=120,
                 NODE_LABEL_POSITION="W,E,c,0.00,0.00")
nodeLabels <- mapVisualProperty('node label','id','p')
nodeFills <- mapVisualProperty('node fill color','node_attributes','d',c("A","B"), c("#FF9900","#66AAAA"))
arrowShapes <- mapVisualProperty('Edge Target Arrow Shape','interaction','d',c("activates","inhibits","interacts"),c("Arrow","T","None"))
edgeWidth <- mapVisualProperty('edge width','weight','p')

createVisualStyle(style.name, defaults, list(nodeLabels,nodeFills,arrowShapes,edgeWidth))
setVisualStyle(style.name)

##################################################### METHOD 2 ####################################################################################
#2.
softPower <- 14
# calling adjacency function
adjacency <- adjacency(norm.counts, power = softPower)

# TOM
TOM <- TOMsimilarity(adjacency)#This gives similarity between genes
TOM.dissimilarity <- 1-TOM # get dissimilarity matrix

# Hierarchical Clustering Analysis
#The dissimilarity/distance measures are then clustered using linkage hierarchical clustering and a dendrogram (cluster tree) of genes is constructed.
#creating the dendrogram 
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average") 
pdf("6.deseq2/5.dendrogram_gene_clustering_TOM_dissimilarity.pdf")
#plotting the dendrogram
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
labels = FALSE, hang = 0.04)
dev.off()

# identify modules
Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
table(Modules) #returns a table of the counts of factor levels in an object. In this case how many genes are assigned to each created module.

#You can now plot the module assignment under the gene dendrogram for visualization
ModuleColors <- labels2colors(Modules) #assigns each module number a color
table(ModuleColors) #returns the counts for each color (aka the number of genes within each module)

#plots the gene dendrogram with the module colors
pdf("6.deseq2/6.Gene_dendrogram_with_modulecolors.pdf")
plotDendroAndColors(geneTree, ModuleColors,"Module",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")
dev.off()

# Module eigengene identification
MElist <- moduleEigengenes(norm.counts, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)

#Module merging
#To further condense the clusters (branches) into more meaningful modules you can cluster modules based on pairwise eigengene correlations and merge the modules that have similar expression profiles.
ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity
METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
pdf("6.deseq2/7.Cluster_dendrogram.pdf")
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic
plot(METree)
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75
dev.off()

merge <- mergeCloseModules(norm.counts, ModuleColors, cutHeight = .25)
# The merged module colors, assigning one color to each module
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs

# dendrogram with original and merged modules
pdf("6.deseq2/8.original_and_merged_modules_dendrogram.pdf")
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
c("Original Module", "Merged Module"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors for original and merged modules")
dev.off()

# External trait matching
# pull out continuous traits
allTraits <- sample_metadata[,c(3:7)]
# Define numbers of genes and samples
nGenes = ncol(norm.counts)
nSamples = nrow(norm.counts)
module.trait.correlation = cor(mergedMEs, allTraits, use = "p") #p for pearson correlation coefficient 
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples) #calculate the p-value associated with the correlation

# create module-trait heatmap
# Will display correlations and their p-values
pdf("6.deseq2/9.Module-trait_relationships.pdf")
textMatrix = paste(signif(module.trait.correlation, 2), "\n(",
signif(module.trait.Pvalue, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.correlation)
par(mar = c(6, 8.5, 3, 1))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = module.trait.correlation,
xLabels = names(allTraits),
yLabels = names(mergedMEs),
ySymbols = names(mergedMEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.4,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()
#Each row corresponds to a module eigengene, and the columns correspond to a trait. 
#Each cell contains a p-value and correlation. Those with strong positive correlations are shaded a darker red while those with stronger negative correlations become more blue.

#Target gene identification

#Network Visualization of Eigengenes, to study the relationship among found modules
# Isolate desired variable
metpro = as.data.frame(sample_metadata$CH4production);
names(metpro) = "methaneproduction"
# Add thevariable to existing module eigengenes
MET = orderMEs(cbind(MEs, metpro))
# Plot the relationships among the eigengenes and the trait
pdf("6.deseq2/10.Network_eigengenes.pdf")
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(5,4,1,2), cex.lab = 0.8, xLabelsAngle
= 90)
dev.off()

# eigengene dendrogram
pdf("6.deseq2/11.eigengene_dendrogram.pdf")
# Plot the dendrogram
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
plotHeatmaps = FALSE)
dev.off()

# eigengene adjacency heatmap
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
pdf("6.deseq2/12.eigengene_adjacency_heatmap.pdf")
par(cex = 1.0, mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(5,5,2,2),
plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()



