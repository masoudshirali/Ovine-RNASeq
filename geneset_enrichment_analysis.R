# Import the annotation file for the ovis aries genome

gffpath="/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0"
gfffile="/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/genomic.gff"

# Load the necessary libraries

library(stringr)
library(clusterProfiler)
library(pathview)
library(stringr)
library(AnnotationHub)
library(ggridges)
library(enrichplot)

# extract annotation DB of sheep

ah <- AnnotationHub()
AnnotationHub::query(ah, c("Ovis", "aries"))
Oaries <- ah[["AH111978"]]
columns(Oaries) # see available identifiers that can be used in this package

# create geneList, which is a list of all genes from deseq2 with p<0.1
# The variable res is the result from deseq2

entrez<-read.csv("/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/EntrezIDsfromGtf",sep="\t", header=TRUE,row.names=NULL)#35105
#modify column names
colnames(entrez) <- colnames(entrez)[2:ncol(entrez)]
#drop last column
entrez <- entrez[1:(ncol(entrez)-1)]
# # Some of the GeneIDs and symbols contain "LOC" info, we need to remove it to merge it with the res and resSig files from deseq2.
entrez$GeneID <- stringr::str_remove(entrez$GeneID, "LOC")
entrez$GeneSymbol <- stringr::str_remove(entrez$GeneSymbol, "LOC")

res_allgenes<-res
res_allgenes$GeneID<-rownames(res)
res_degs<-resSig
res_degs$GeneID<-rownames(resSig)
res_degs$GeneSymbol<-rownames(resSig)

# removing "LOC" info
res_allgenes$GeneID <- stringr::str_remove(res_allgenes$GeneID, "LOC")
rownames(res_allgenes) = stringr::str_remove(rownames(res_allgenes), "LOC")
res_degs$GeneID <- stringr::str_remove(res_degs$GeneID, "LOC")
rownames(res_degs) = stringr::str_remove(rownames(res_degs), "LOC")

# res_allgenes_entrez<-merge(as.data.frame(res_allgenes), entrez, by="GeneID") # This significantly reduced the number of mapping gene IDS from gtf file, 35,117 genes were reduced to 16,120
# res_degs_entrez<-merge(as.data.frame(res_degs), entrez, by="GeneSymbol") # same here
# Therefore trying another approach to retrieve the ENTREZ IDs

# get entrezIds for the genes with symbols
res_allgenesEntrez <- select(Oaries, keys =  res_allgenes$GeneID,
  columns = c('ENTREZID'), keytype = 'SYMBOL')

res_sigGenesEntrez <- select(Oaries, keys = res_degs$GeneID,
  columns = c('ENTREZID'), keytype = 'SYMBOL')

# Now replace the NA values in entrezid column with the values from first column. This is done because many values were not converted to ENTREZ IDs but already has it from GTF file. We will retain those.
res_allgenesEntrez$ENTREZID <- ifelse(is.na(res_allgenesEntrez$ENTREZID), res_allgenesEntrez$SYMBOL, res_allgenesEntrez$ENTREZID)
colnames(res_allgenesEntrez) = c("GeneID", "ENTREZID")
res_allgenes_with_entrez = merge(as.data.frame(res_allgenes), res_allgenesEntrez, by = "GeneID") #All entrezIDs have been retrieved for the DEGs (except for 2 genes)

res_sigGenesEntrez$ENTREZID <- ifelse(is.na(res_sigGenesEntrez$ENTREZID), res_sigGenesEntrez$SYMBOL, res_sigGenesEntrez$ENTREZID)
colnames(res_sigGenesEntrez) = c("GeneID", "ENTREZID")
res_degs_with_entrez = merge(as.data.frame(res_degs), res_sigGenesEntrez, by = "GeneID") #All entrezIDs have been retrieved for the DEGs (except for 2 genes)

# In order to asses functional enrichment, both DE gene list and gene universe must be annotated in Entrez IDs:

res_universe<-res_allgenes_with_entrez$ENTREZID
# omit any NA values 
res_universe<-na.omit(res_universe)
# sort the list in decreasing order (required for clusterProfiler)
res_universe = sort(res_universe, decreasing = TRUE)
res_sigGenes<-res_degs_with_entrez$ENTREZID
res_sigGenes = sort(res_sigGenes, decreasing = TRUE)

ans.go <- enrichGO(gene = res_sigGenes, ont = "ALL",
                   OrgDb =Oaries,
                   universe = res_universe,
                   readable=TRUE,
                   pvalueCutoff = 0.1)
tab.go <- as.data.frame(ans.go)
write.csv(tab.go,"GO_enrichments.csv")

ans.kegg <- enrichKEGG(gene = res_sigGenes,
                       organism = 'oas',
                       universe = res_universe,
                       pvalueCutoff = 0.1)
tab.kegg <- as.data.frame(ans.kegg)
write.csv(tab.kegg,"KEGG_enrichments.csv")


# Visualizations of the GO analysis

prefix1 = "DGE_GO"

pdf("GO-Barplot.pdf")
barplot(ans.go, showCategory=10)
dev.off()

pdf("upsetplot_kegg.pdf")
upsetplot(ans.kegg)
dev.off()

pdf("emapplot_kegg.pdf")
emapplot(pairwise_termsim(ans.kegg))
dev.off()






























# GO classification
# the groupGO() function is designed for gene classification based on GO distribution at a specific level.
ggo <- groupGO(gene     = res_sigGenes,
               OrgDb    = Oaries,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)

# GO over-representation analysis
ego <- enrichGO(gene          = res_sigGenes,
                universe      = names(res_universe),
                OrgDb         = Oaries,
                ont           = "ALL",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.09,
                qvalueCutoff  = 0.09,
        readable      = TRUE)
head(ego)
