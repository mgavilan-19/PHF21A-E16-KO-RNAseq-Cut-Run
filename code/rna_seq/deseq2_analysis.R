# Read a text file into a data frame --- deseq coutns
countData <- read.table("ctx_e16_count.txt", header = TRUE, row.names = 1, sep = "\t")
#### Creating colData
colData <- data.frame(
  Samples = c("CWT_NN_1.bam", "CWT_NN_2.bam", "CWT_NN_3.bam", "CKO_NN_1.bam", "CKO_NN_2.bam", "CKO_NN_3.bam",
              "NWT_NN_1.bam", "NWT_NN_2.bam", "NWT_NN_3.bam", "NKO_NN_1.bam", "NKO_NN_2.bam", "NKO_NN_3.bam"),
  Conditions = c("canonical", "canonical", "canonical", "canonical", "canonical", "canonical",
                "neuronal","neuronal","neuronal","neuronal","neuronal","neuronal"),
  Treatment = c("WT", "WT", "WT", "KO", "KO", "KO"),
  stringsAsFactors = FALSE  # Avoid converting strings to factors
)
# Set the Samples column as row names
rownames(colData) <- colData$Samples
# Remove the Samples column from the data frame
colData$Samples <- NULL
# Print the updated data frame
print(colData)

## Make sure that rownames and colnames match
all(colnames(countData) %in% rownames(colData))
all(colnames(countData) == rownames(colData))

# Only analyzing canonical data
canonical_subset <- colData$Condition == "canonical"
# Subset the counts and metadata for only the canonical condition
counts_canonical <- countData[, canonical_subset]
colData_canonical <- colData[canonical_subset, ]

## To set deseq analysis KO vs WT
colData_canonical$Treatment <- factor(colData_canonical$Treatment, levels = c("WT", "KO"))

library("DESeq2")
dds_canonical <- DESeqDataSetFromMatrix(
  countData = counts_canonical,
  colData = colData_canonical,
  design = ~ Treatment
)
dds_canonical

#Run DESeq
dds_canonical <- DESeq(dds_canonical)
res_canonical <- results(dds_canonical)
summary(res_canonical)

# To normalize it using lfcShrink
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm", force = TRUE)

# Shrink log2 fold changes
resLFC <- lfcShrink(dds_canonical, coef="Treatment_KO_vs_WT", type="apeglm")
summary(resLFC)

# Convert this to a DataFrame 
resLFC <- as.data.frame(resLFC)

# Added column to classify up or down regulated genes
resLFC$diffexpressed <- "NO"
# If log2FoldChange > 0 and pvalue < 0.05, set as "UP"
resLFC$diffexpressed[resLFC$log2FoldChange > 0 & resLFC$padj < 0.05] <- "UP"
# If log2FoldChange < - 0 and pvalue < 0.05, set as "DOWN"
resLFC$diffexpressed[resLFC$log2FoldChange < -0 & resLFC$padj < 0.05] <- "DOWN"
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
resLFC$gene_id <- NA
resLFC$X <- rownames(resLFC)
resLFC$gene_id[resLFC$diffexpressed != "NO"] <- resLFC$X[resLFC$diffexpressed != "NO"]

##### To add the gene ID = symbol
BiocManager::install("clusterProfiler", version = "3.18", force = TRUE)
BiocManager::install("GOSemSim", force = TRUE) #its a clusterprofiler dependency
BiocManager::install("pathview", force = TRUE)
BiocManager::install("enrichplot", force = TRUE)

library(tidygraph)
library(clusterProfiler)
library(enrichplot)

# SET THE DESIRED ORGANISM HERE
organism = "org.Mm.eg.db" # This is for mouse. Here we can find other organisms: https://bioconductor.org/packages/release/BiocViews.html#___OrgDb
BiocManager::install(organism, character.only = TRUE)

library(organism, character.only = TRUE)

if (!("msigdbr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("msigdbr", update = FALSE)
}
library(msigdbr)

convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select( 
    db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey) ) )
  if( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]   
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

resLFC$hgnc_symbol <- convertIDs( row.names(resLFC), "ENSEMBL", "SYMBOL", org.Mm.eg.db )
resLFC$entrezid <- convertIDs( row.names(resLFC), "ENSEMBL", "ENTREZID", org.Mm.eg.db )

genes_per_DEG <- table(resLFC$diffexpressed)
print(genes_per_DEG)
#DOWN    NO    UP 
  247 53247   146 

 
### Save the document
# Export file 
write.csv(as.data.frame(resLFC), file = "2024-03-13-deseq_CE16_CKO_WT.csv")




