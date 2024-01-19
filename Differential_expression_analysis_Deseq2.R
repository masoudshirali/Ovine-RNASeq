# This analysis takes the counts data generated from featureCounts and the metadata with the phenotype information as the input files.
# So make sure you have both files in the correct directories. In this scenario, I have copied the metadata to the main working directory (/mnt/sda1/RNA/40-815970407/Sheep)

# Load all the necessary libraries. If not installed, install using either:
install.packages("package") or
# to install biconductor packages, 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install(c("clusterProfiler", "AnnotationHub"))

library(DESeq2)
library(dplyr)
library(tidyverse)

# set the working directory
setwd("/mnt/sda1/RNA/40-815970407/Sheep")

# create a directory to store DESeq2 outputs
system("mkdir Deseq2_results")

# Read the raw count data generated from featureCounts
countData<-read.csv("4.featurecounts/LambAllSamples.featureCounts.RmatrixNew.csv",sep=",", header=T, check.names=F)
countData<-countData[ , !names(countData) %in% c("Low7085.bam","Medium7073.bam")]# Remove 7085 and 7073 as they had poor mapping rates
colnames(countData)<-gsub(".bam","",colnames(countData))# Remove the .bam from the column names
orig_names <- names(countData) # keep a back-up copy of the original names
geneID <- countData$Geneid# Convert count data to a matrix of appropriate form that DEseq2 can read
countData <- as.matrix(countData[ , -1]) # removing first column geneID from the table
# make sure the rownames are gene ids and first column onwards should be samples. any other columns should be removed.otherwise deseq2 will throw error.
sampleIndex <- colnames(countData)
countData <- as.matrix(countData[,sampleIndex])
rownames(countData) <- geneID
head(countData)

# check library size (total mapped reads)
colSums(countData[,2:ncol(countData)])

# Read the metadata file
metaData1 <-read.csv("metadata_with_methaneinfoadded_metadata.csv",sep=",",header=T)
dim(metaData1)
head(metaData1)
rownames(metaData1) <- metaData1$ID
metaData1$ID <- factor(metaData1$ID)
head(metaData1)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
countData <- countData[,unique(rownames(metaData1))]
all(colnames(countData) == rownames(metaData1))

deseq2Data <- DESeqDataSetFromMatrix(countData=countData, colData=metaData1, design= ~CH4production)
deseq2Data <- deseq2Data[rowSums(assay(deseq2Data)) > 0, ]
dim(deseq2Data)
keep <- rowSums(counts(deseq2Data)) >= 4
deseq2Data <- deseq2Data[keep,]
dim(deseq2Data)

deseq2Data <- DESeq(deseq2Data)

#loop through results and extract significant DEGs for each model term
results = resultsNames(deseq2Data)
upresultstable = matrix(nrow = length(results), ncol = 1, dimnames = list(results,"upDEGs"))
downresultstable = matrix(nrow = length(results), ncol = 1, dimnames = list(results,"downDEGs"))

for(i in 1:length(results)){

  res = results(deseq2Data, 
                name = results[i])
  resorder <- res[order(res$padj),]
  upDEGs = (length(na.omit(which(res$padj<0.1 & res$log2FoldChange > 0))))
  downDEGs = (length(na.omit(which(res$padj<0.1 & res$log2FoldChange < 0))))
  resSig = subset(resorder, padj < 0.1 & log2FoldChange > 0 | padj < 0.1 & log2FoldChange < 0)
  write.csv(resSig , file=paste0(results[i],".0.1p.lfc0.updownDEGs.csv"), row.names = T)
  upresultstable[results[i],"upDEGs"] = upDEGs
  downresultstable[results[i],"downDEGs"] = downDEGs 

}

