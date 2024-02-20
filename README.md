# OarRNASeq
RNA-Seq analysis pipeline for Ovis aries samples.

# Pre-requisites
Make sure you have installed the necessary softwares required to run the pipeline. These are mentioned in detail below.

1. FastQC
2. Trimmomatic
3. Hisat2
4. featureCounts

After you have installed the necessary software, please proceed to running the actual pipeline. We would run them in the following order:
1. bash RNASeqPipeline.sh or /.RNASeqPipeline.sh
2. Rscript Deseq2.R
3. Rscript geneset_enrichment_analysis.R
4. Rscript wgcna.R
5. Rscript Stringdb.R
