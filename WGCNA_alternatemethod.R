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
  
  table(gsg$goodGenes)
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

## remove all genes with counts < 10 in more than 75% of samples (22*0.75=16)
## suggested by WGCNA on RNAseq FAQ

  dds75 <- dds[rowSums(counts(dds) >= 10) >= 16,]
  nrow(dds75) # 16583 genes


# perform variance stabilization
  dds_norm <- vst(dds75)
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

#######################################################################################
# Network construction
#######################################################################################
  soft_power <- 14
  temp_cor <- cor
  cor <- WGCNA::cor


# memory estimate w.r.t blocksize
  bwnet <- blockwiseModules(norm.counts,
                 maxBlockSize = 7000,
                 TOMType = "signed",
                 power = soft_power,
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE,
                 randomSeed = 1234,
                 verbose = 3) #mergecutHeight value: #For fewer samples a larger valuse (0.25 to 0.3) may be warranted


  cor <- temp_cor

##########################################################################################
# Module eigengenes
##########################################################################################
  module_eigengenes <- bwnet$MEs

# Print out a preview
  head(module_eigengenes)

# get number of genes for each module
  table(bwnet$colors)#18 modules total

# Plot the dendrogram and the module colors before and after merging underneath
  pdf("7.wgcna/original.and.merged.dendrograms.pdf")
  plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)
  dev.off()



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
  zlim = c(-1,1),
  main = paste("Module-trait relationships"))
  dev.off()

heatmap.data <- merge(mergedMEs , allTraits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

pdf("sample.pdf")
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[17:19],
             y = names(heatmap.data)[1:16],
             col = c("blue1", "skyblue", "white", "pink", "red"))
dev.off()
