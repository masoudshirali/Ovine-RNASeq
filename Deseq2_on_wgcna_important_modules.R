# extract genes from the inetersting modules
yg_pvr = geneInfo %>% filter(moduleColor %in% "yellowgreen" | moduleColor %in% "palevioletred3")
yellowgreen = geneInfo %>% filter(moduleColor %in% "yellowgreen")
palevioletred3 = geneInfo %>% filter(moduleColor %in% "palevioletred3")

# extract the raw counts for these genes in the modules
ygLfcRaw <- countData[rownames(countData) %in% yellowgreen$Genes,]
pvrLfcRaw <- countData[rownames(countData) %in% palevioletred3$Genes,]
yg_pvr_Raw = rbind(ygLfcRaw,pvrLfcRaw)

# removing "LOC" info
rownames(yg_pvr_Raw) = stringr::str_remove(rownames(yg_pvr_Raw), "LOC")

# Read the metadata file
metaData <-read.csv("metadata_with_methaneinfoadded_metadata.csv",sep=",",header=T)
#metaData1 = metaData[!grepl('Low|Medium', metaData$ID),]# To remove any rows with Low or medium in the ID column. (Only for test case)
#metaData = metaData1
dim(metaData)
head(metaData)
rownames(metaData) <- metaData$ID
metaData$ID <- factor(metaData$ID)
rownames(metaData)<-gsub("[a-zA-Z ]", "", rownames(metaData))
head(metaData)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
yg_pvr_Raw <- yg_pvr_Raw[,unique(rownames(metaData))]
all(colnames(yg_pvr_Raw) == rownames(metaData))

deseq2Data <- DESeqDataSetFromMatrix(countData=yg_pvr_Raw, colData=metaData, design= ~CH4production)

# Stringent approach where we keep only rows that have at least 10 reads total
keep <- rowSums(counts(deseq2Data)) >= 10
deseq2Data <- deseq2Data[keep,] #240 remaining

deseq2Data <- DESeq(deseq2Data)


#loop through results and extract significant DEGs for each model term
# speify the cut-offs for pval and lfc in the below variables.
# make sure to change the filenames with the cutoff values before saving the deg file (Line 60)

# STRINGENT - with reads removed in steps 32 to 34
# RELAXED - with reads not removed

pval = 0.1
lfc = 0.584
results = resultsNames(deseq2Data)
upresultstable = matrix(nrow = length(results), ncol = 1, dimnames = list(results,"upDEGs"))
downresultstable = matrix(nrow = length(results), ncol = 1, dimnames = list(results,"downDEGs"))

for(i in 1:length(results)){

  res = results(deseq2Data, 
                name = results[i])
  resorder <- res[order(res$padj),]
  upDEGs = (length(na.omit(which(res$padj<pval & res$log2FoldChange > lfc))))
  downDEGs = (length(na.omit(which(res$padj<pval & res$log2FoldChange < -lfc))))
  resSig = subset(resorder, padj < pval & log2FoldChange > lfc | padj < pval & log2FoldChange < -lfc)
  write.csv(resSig , file=paste0("7.wgcna/",results[i],".0.1P.0.584LFC.updownDEGs_STRINGENT.csv"), row.names = T)
  upresultstable[results[i],"upDEGs"] = upDEGs
  downresultstable[results[i],"downDEGs"] = downDEGs 
}

pval = 0.1
lfc = 0.584
results = resultsNames(deseq2Data)
upresultstable = matrix(nrow = length(results), ncol = 1, dimnames = list(results,"upDEGs"))
downresultstable = matrix(nrow = length(results), ncol = 1, dimnames = list(results,"downDEGs"))

for(i in 1:length(results)){

  res = results(deseq2Data, 
                name = results[i])
  resorder <- res[order(res$padj),]
  upDEGs = (length(na.omit(which(res$padj<pval & res$log2FoldChange > lfc))))
  downDEGs = (length(na.omit(which(res$padj<pval & res$log2FoldChange < -lfc))))
  resSig = subset(resorder, padj < pval & log2FoldChange > lfc | padj < pval & log2FoldChange < -lfc)
  write.csv(resSig , file=paste0("7.wgcna/",results[i],".0.1P.0.584LFC.updownDEGs_RELAXED.csv"), row.names = T)
  upresultstable[results[i],"upDEGs"] = upDEGs
  downresultstable[results[i],"downDEGs"] = downDEGs 
}

