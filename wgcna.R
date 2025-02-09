# Load the libraries

library(WGCNA)
library(flashClust)
library(curl)
library(DESeq2)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(tidyverse)
library(CorLevelPlot)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
allowWGCNAThreads()          # allow multi-threading (optional)

##################################################################################################
# Read the gene counts table and metadata
##################################################################################################

  data=read.csv("5.featurecounts/Lambs.featurecounts.hisat2.Rmatrix",header=T,row.names=1,sep="\t", check.names = FALSE)
  data=data[ , !names(data) %in% c("7085","7073")]
  colnames(data)<-gsub("Control","",colnames(data))
  colnames(data)<-gsub("Low","",colnames(data))
  colnames(data)<-gsub("Medium","",colnames(data))
  colnames(data)<-gsub("High","",colnames(data))
  
  # Read the metadata
  sample_metadata = read.csv(file = "metadata_with_methaneinfoadded_metadata.csv")
  #sample_metadata$X<-NULL
  rownames(sample_metadata) <- sample_metadata$ID
  sample_metadata$ID <- factor(sample_metadata$ID)
  rownames(sample_metadata)<-gsub("[a-zA-Z ]", "", rownames(sample_metadata))

###########################################################################################
# QC - outlier detection
###########################################################################################
  
# detect outlier genes
  gsg <- goodSamplesGenes(t(data))
  summary(gsg)
  gsg$allOK
  
  table(gsg$goodGenes) # 25030 genes passed
  table(gsg$goodSamples)
  
# if allOK returen false, remove genes that are detectd as outliers
  data <- data[gsg$goodGenes == TRUE,]
  
# detect outlier samples - hierarchical clustering - method 1
  pdf("7.wgcna/1.hclust_samples.pdf")
  sampleTree <- hclust(dist(t(data)), method = "average") #Clustering samples based on distance 
  #Setting the graphical parameters
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
# Plotting the cluster dendrogram
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
  cex.axis = 1.5, cex.main = 2)
  dev.off()
  
# detect outlier samples - pca - method 2
  pca <- prcomp(t(data))
  pca.dat <- pca$x
  
  pca.var <- pca$sdev^2
  pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
  
  pca.dat <- as.data.frame(pca.dat)
  pdf("7.wgcna/2.pca.pdf")
  ggplot(pca.dat, aes(PC1, PC2)) +
    geom_point() +
    geom_text(label = rownames(pca.dat)) +
    labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
         y = paste0('PC2: ', pca.var.percent[2], ' %'))
  dev.off()
  
# exclude outlier samples
#samples.to.be.excluded <- c('GSM4615000', 'GSM4614993', 'GSM4614995')
#data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]
  
#####################################################################################
# Normalization
###################################################################################
  
# create a deseq2 dataset
# making the rownames and column names identical
# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
  data <- data[,unique(rownames(sample_metadata))]
  all(colnames(data) == rownames(sample_metadata))

# create dds
  dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = sample_metadata,
                              design = ~ 1) # not spcifying model

## remove all genes with counts < 4 in more than 75% of samples (22*0.75=16)
## suggested by WGCNA on RNAseq FAQ

  dds75 <- dds[rowSums(counts(dds) >= 4)]
  nrow(dds75) 


# perform variance stabilization
  dds_norm <- vst(dds75)
  dds_norm <- vst(dds)
  write.csv(assay(dds_norm),"7.wgcna/Lambs_allSamples_normalized_Counts",row.names=T)

# get normalized counts
  norm.counts <- assay(dds_norm) %>% 
  t()

########################################################################################
# Network Construction 
########################################################################################

# Choose a set of soft-thresholding powers
  power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
  sft <- pickSoftThreshold(norm.counts,
                  powerVector = power,
                  networkType = "signed",
                  verbose = 5)

  power=sft$powerEstimate #14

  sft.data <- sft$fitIndices

# visualization to pick power
  pdf("7.wgcna/3.power_threshold.pdf")
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

  softPower <- 14
# calling adjacency function
  adjacency <- adjacency(norm.counts, power = softPower, type="signed")

# TOM
  TOM <- TOMsimilarity(adjacency, TOMType = "signed")#This gives similarity between genes
  TOM.dissimilarity <- 1-TOM # get dissimilarity matrix

# Hierarchical Clustering Analysis
#The dissimilarity/distance measures are then clustered using linkage hierarchical clustering and a dendrogram (cluster tree) of genes is constructed.
#creating the dendrogram 
  geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average") 

#plotting the dendrogram
  pdf("7.wgcna/4.dendrogram_gene_clustering_TOM_dissimilarity.pdf")
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
  labels = FALSE, hang = 0.04)
  dev.off()

##############################################################################
# identify modules
##############################################################################

  Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
  table(Modules) #returns a table of the counts of factor levels in an object. In this case how many genes are assigned to each created module.

#You can now plot the module assignment under the gene dendrogram for visualization
  ModuleColors <- labels2colors(Modules) #assigns each module number a color
  table(ModuleColors) #returns the counts for each color (aka the number of genes within each module) #53 modules

#plots the gene dendrogram with the module colors
  pdf("7.wgcna/5.Gene_dendrogram_with_modulecolors.pdf")
  plotDendroAndColors(geneTree, ModuleColors,"Module",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = "Gene dendrogram and module colors")
  dev.off()

# Module eigengene identification
  MElist <- moduleEigengenes(norm.counts, colors = ModuleColors) 
  MEs <- MElist$eigengenes 
  head(MEs)

##############################################################################
# Module merging
##############################################################################

#To further condense the clusters (branches) into more meaningful modules you can cluster modules based on pairwise eigengene correlations and merge the modules that have similar expression profiles.
  ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity
  METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
  pdf("7.wgcna/6.Clustering of module eigengenes.pdf")
  par(mar = c(0,4,2,0)) #seting margin sizes
  par(cex = 0.6);#scaling the graphic
  plot(METree)
  abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75
  dev.off()

  merge <- mergeCloseModules(norm.counts, ModuleColors, cutHeight = .25)

# The merged module colors, assigning one color to each module
  mergedColors = merge$colors
  table(mergedColors) #41 modules

# Eigengenes of the new merged modules
  mergedMEs = merge$newMEs

# dendrogram with original and merged modules
  pdf("7.wgcna/7.original_and_merged_modules_dendrogram.pdf")
  plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
  c("Original Module", "Merged Module"),
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = "Gene dendrogram and module colors for original and merged modules")
  dev.off()

  write.table(merge$oldMEs,file="7.wgcna/oldMEs.txt");
  write.table(merge$newMEs,file="7.wgcna/newMEs.txt");

####################################################################################
# External trait matching
####################################################################################

# pull out all continuous traits
  allTraits <- sample_metadata[,c(3:21)]

# Define numbers of genes and samples
  nGenes = ncol(norm.counts)
  nSamples = nrow(norm.counts)
  module.trait.correlation = cor(mergedMEs, allTraits, use = "p") #p for pearson correlation coefficient 
  module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples) #calculate the p-value associated with the correlation
  write.table(module.trait.correlation,file="7.wgcna/moduleTrait_correlation.txt");
  write.table(module.trait.Pvalue,file="7.wgcna/moduleTrait_pValue.txt");

# create module-trait heatmap
# Will display correlations and their p-values
  pdf("7.wgcna/9.Module-trait_relationships.pdf")
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
  cex.axis = 2.5,
  zlim = c(-1,1),
  main = paste("Module-trait relationships"))
  dev.off()

#Each row corresponds to a module eigengene, and the columns correspond to a trait. 
#Each cell contains a p-value and correlation. Those with strong positive correlations are shaded a darker red while those with stronger negative correlations become more blue.

# Which modules have biggest differences across methane production or methane yield?
# Create the design matrix from the `methane yield` variable  
  des_mat <- model.matrix(~sample_metadata$CH4yield)
# Run linear model on each module. Limma wants our tests to be per row, so we also need to transpose so the eigengenes are rows
  fit <- limma::lmFit(t(mergedMEs), design = des_mat)
# Apply empirical Bayes to smooth standard errors
  fit <- limma::eBayes(fit)
# Apply multiple testing correction and obtain stats  
stats_df <- limma::topTable(fit, number = ncol(mergedMEs)) %>%
  tibble::rownames_to_column("module")# It seems darkmagenta module has highest significance to methane yield. Lets dive into more details of this module.  
  
# heatmap with significance
  heatmap.data <- merge(mergedMEs , allTraits, by = 'row.names')
  head(heatmap.data)
  heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')
  pdf("7.wgcna/9.Module-trait_relationships_with_significance.pdf", width=14, height=10)
  CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[42:60],
             y = names(heatmap.data)[1:41],
             col = c("blue1", "skyblue", "white", "pink", "red"),
             rotLabX = 30, rotLabY = 30)
  dev.off()

#####################################################################################################
# Intramodular analysis: identifying genes with high geneModuleMembership & geneTraitSignficance
# Target gene identification
#####################################################################################################

# Define variable weight containing the weight column of datTrait
  metpro = as.data.frame(sample_metadata$CH4production);
  names(metpro) = "methane_production"
  metyield = as.data.frame(sample_metadata$CH4yield)
  names(metyield) = "methane_yield"
  modNames = substring(names(mergedMEs), 3) #extract module names

#Calculate the module membership and the associated p-values
  geneModuleMembership = as.data.frame(cor(norm.counts, mergedMEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")

#Calculate the gene significance and associated p-values for methane production
  geneTraitSignificance = as.data.frame(cor(norm.counts, metpro, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) = paste("GS.", names(metpro), sep="")
  names(GSPvalue) = paste("p.GS.", names(metpro), sep="")
  head(GSPvalue)
  GSPvalue.sig.metpro = subset(GSPvalue, p.GS.methane_production<0.05)#654 genes that have a high significance for methane production
  GSPvalue %>%
  as.data.frame() %>%
  arrange(p.GS.methane_production) %>%
  head(25) #displays top 25 significant genes associated with methane production
  write.csv(GSPvalue.sig.metpro,"7.wgcna/Genes_with_high_significance_to_methaneproduction.csv", row.names=TRUE)

#Calculate the gene significance and associated p-values for methane yield
  geneTraitSignificance_y = as.data.frame(cor(norm.counts, metyield, use = "p"))
  GSPvalue_y = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_y), nSamples))
  names(geneTraitSignificance_y) = paste("GS.", names(metyield), sep="")
  names(GSPvalue_y) = paste("p.GS.", names(metyield), sep="")
  head(GSPvalue_y)
  GSPvalue.sig.metyield = subset(GSPvalue_y, p.GS.methane_yield<0.05)#450 genes that have a high significance for methane yield
  GSPvalue_y %>%
  as.data.frame() %>%
  arrange(p.GS.methane_yield) %>%
  head(25) #displays top 25 significant genes associated with methane production
  write.csv(GSPvalue.sig.metyield,"7.wgcna/Genes_with_high_significance_to_methaneyield.csv", row.names=TRUE)

# we have highest significance for methane production in darkmagenta module
# Plot a scatter plot of gene significance vs. module membership in the darkmagenta module.
# calculate the module membership values (aka. module eigengene based connectivity kME):
  datKME = signedKME(norm.counts, mergedMEs)
  pdf("7.wgcna/10.genesignificance_vs_modulemembership_yellowgreen_modules.pdf",width=14)
  colorOfColumn = substring(names(datKME), 4)
  par(mfrow = c(2, 1))
  selectModules = c("yellowgreen")
  par(mfrow = c(1, 2))
  for (module in selectModules) {
    column = match(module, colorOfColumn)
    restModule = mergedColors == module
    verboseScatterplot(datKME[restModule, column], abs(geneTraitSignificance[restModule,1]), xlab = paste("Module Membership ", 
        module, "module"), ylab = "GS.methane_production", main = paste("kME.", module, 
        "vs. GS"), col = module)
    verboseScatterplot(datKME[restModule, column], abs(geneTraitSignificance_y[restModule,1]), xlab = paste("Module Membership ", 
        module, "module"), ylab = "GS.methane_yield", main = paste("kME.", module, 
        "vs. GS"), col = module)
  }
  dev.off()

# This indicates that the genes that are highly significantly associated with the trait (high gene significance) are also the genes that are the most connected within their module (high module membership). 
# Therefore genes in the darkmagenta module could be potential target genes when looking at methane production.

##################################################################################################
#Intramodular analysis: identifying genes with high GS and MM
##################################################################################################
#identifying genes with high GS and MM for significant module (GS > 0.2 & MM > 0.8)
# 1. methane production trait
  intra_modular_analysis=data.frame(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]))
  rownames(intra_modular_analysis) = colnames(norm.counts)[ModuleColors=="yellowgreen"] #only the yellowgreen module
  head(intra_modular_analysis)
  colnames(intra_modular_analysis)<- c("abs.geneModuleMembership.moduleGenes", "abs.geneTraitSignificance.moduleGenes")
  intra_modular_analysis.hubgene = subset(intra_modular_analysis, abs.geneModuleMembership.moduleGenes>0.8 & abs.geneTraitSignificance.moduleGenes > 0.2)
  write.csv(intra_modular_analysis.hubgene, file = "7.wgcna/hubgenes_in_yellowgreen_methaneproduction.csv")

# 1. methane yield trait
  intra_modular_analysis_y=data.frame(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance_y[moduleGenes, 1]))
  rownames(intra_modular_analysis_y) = colnames(norm.counts)[ModuleColors=="palevioletred3"] #only the palevioletred3 module
  head(intra_modular_analysis_y)
  colnames(intra_modular_analysis_y)<- c("abs.geneModuleMembership.moduleGenes", "abs.geneTraitSignificance.moduleGenes")
  intra_modular_analysis.hubgene_y = subset(intra_modular_analysis_y, abs.geneModuleMembership.moduleGenes>0.8 & abs.geneTraitSignificance.moduleGenes > 0.2)
  write.csv(intra_modular_analysis.hubgene_y, file = "7.wgcna/hubgenes_in_palevioletred3_methaneyield.csv")

#high intramodular connectivity ~ high kwithin => hub genes (kwithin: connectivity of the each driver gene in the darkmagenta module to all other genes in the darkmagenta)
  connectivity = intramodularConnectivity(adjacency, ModuleColors)
  connectivity = connectivity[colnames(norm.counts)[ModuleColors=="yellowgreen"],] #only the darkmagenta module
  order.kWithin = order(connectivity$kWithin, decreasing = TRUE)
  connectivity = connectivity[order.kWithin,] #order rows following kWithin
  #connectivity = connectivity[1:5,] #top 5 genes that have a high connectivity to other genes in the yellowgreen module
  connectivity %>% head(10)#top 10 genes that have a high connectivity to other genes in the yellowgreen module
  write.csv(connectivity,"7.wgcna/hubgenes_ordered_with_high_connectivity_yellowgreenmodule.csv")#top genes would be the one with high connectivity

#Identifying most important genes for one determined characteristic inside of the cluster
  probes = colnames(norm.counts)
  geneInfo0 = data.frame(Genes = probes,
                       moduleColor = mergedColors,
                       geneTraitSignificance,
                       GSPvalue)

  geneInfo1 = data.frame(Genes = probes,
                       moduleColor = mergedColors,
                       geneTraitSignificance_y,
                       GSPvalue_y)
#Order modules by their significance for trait
  modOrder = order(-abs(cor(mergedMEs, metpro, use = "p")))
  modOrder1 = order(-abs(cor(mergedMEs, metyield, use = "p")))
  for (mod in 1:ncol(geneModuleMembership))
  {
  oldNames = names(geneInfo0)
  oldNames1 = names(geneInfo1)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
  geneInfo1 = data.frame(geneInfo1, geneModuleMembership[, modOrder1[mod]], 
                         MMPvalue[, modOrder1[mod]]);
  names(geneInfo1) = c(oldNames1, paste("MM.", modNames[modOrder1[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder1[mod]], sep=""))
 }
  geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.methane_production))
  geneOrder1 = order(geneInfo1$moduleColor, -abs(geneInfo1$GS.methane_yield))
  geneInfo = geneInfo0[geneOrder, ]
  geneInfo_y = geneInfo1[geneOrder1, ]
  write.csv(geneInfo, file = "7.wgcna/geneInfo_methaneproduction.csv")
  write.csv(geneInfo_y, file = "7.wgcna/geneInfo_methaneyield.csv")

#Hub genes
  hub = chooseTopHubInEachModule(norm.counts, mergedColors)
  write.csv(hub, file = "7.wgcna/hub_genes_in_each_module.csv") 

#######################################################################################
# Network Visualization of Eigengenes, to study the relationship among found modules
#######################################################################################

# Isolate desired variable
  metpro = as.data.frame(sample_metadata$CH4production);
  names(metpro) = "methaneproduction"
# Add the variable to existing module eigengenes
  MET = orderMEs(cbind(MEs, metpro))
# Plot the relationships among the eigengenes and the trait
  pdf("7.wgcna/11.Network_eigengenes_methaneproduction.pdf")
  par(cex = 0.9)
  plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(5,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)
  dev.off()

# eigengene dendrogram
  pdf("7.wgcna/12.eigengene_dendrogram_methaneproduction.pdf")
# Plot the dendrogram
  par(cex = 1.0)
  plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
  plotHeatmaps = FALSE)
  dev.off()

# eigengene adjacency heatmap
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  pdf("7.wgcna/13.eigengene_adjacency_heatmap_methaneproduction.pdf", width=12, height=13)
  par(cex = 1.0, mar = c(1,1,1,1))
  plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(5,5,2,2),
  plotDendrograms = FALSE, xLabelsAngle = 90)
  dev.off()

###################### option 2 for module identification
cor <- WGCNA::cor
net = blockwiseModules(norm.counts, power = power,
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
pdf(file = "7.wgcna/option2-module_tree_blockwise.pdf", width = 8, height = 6);
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#########################################################################################################################################
# Last step is to export and save the network. Then you can import it in a software for network visualization as Cytoscape, for example.
#########################################################################################################################################

# Exporting the network to a cytoscape format
# Export the gene list of old modules 
  for (i in 1:length(merge$oldMEs)){
    modules = c(substring(names(merge$oldMEs)[i], 3));
    genes = colnames(norm.counts)
    inModule = is.finite(match(ModuleColors,modules))
    modGenes = genes[inModule]
    modTOM=TOM[inModule,inModule]
    dimnames(modTOM)=list(modGenes,modGenes)
    cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("7.wgcna/orign_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("7.wgcna/orign_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = ModuleColors[inModule]);
  }

# Export the gene list of new modules 
  for (i in 1:length(merge$newMEs)){
    modules = c(substring(names(merge$newMEs)[i], 3));
    genes = colnames(norm.counts)
    inModule = is.finite(match(ModuleColors,modules))
    modGenes = genes[inModule]
    modTOM=TOM[inModule,inModule]
    dimnames(modTOM)=list(modGenes,modGenes)
    cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("7.wgcna/merge_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("7.wgcna/merge_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = ModuleColors[inModule]);
  }

#####################################################################################################################################
#   Cytoscape
#####################################################################################################################################

#if(!"RCy3" %in% installed.packages()){
#  install.packages("BiocManager")
#  BiocManager::install("RCy3")
#}

# https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/
  library(RCy3)

  cytoscapePing () # make sure cytoscape is open
  cytoscapeVersionInfo ()

###### for yellowgreen module of the merged data (newMEs) #################################
  edge <- read.delim("7.wgcna/merge_CytoscapeInput-edges-darkmagenta.txt")
  colnames(edge)
  colnames(edge) <- c("source", "target","weight","direction","fromAltName","toAltName")

  node <- read.delim("7.wgcna/merge_CytoscapeInput-nodes-darkmagenta.txt")
  colnames(node)  
  colnames(node) <- c("id","altName","node_attributes") 

  createNetworkFromDataFrames(node,edge[1:50,], title="methane production network", collection="DataFrame Example")

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

