# Read a text file into a data frame --- deseq coutns
countData <- read.table("ctx_e16_count.txt", header = TRUE, row.names = 1, sep = "\t")
# header = TRUE indicates that the first line of the TXT file contains the column names.
# row.names = 1 tells R to use the first column of the data as row names

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
## TRUE
# Find which column names in countData are not in the row names of colData
########################  missing_names <- setdiff(colnames(countData), rownames(colData))
# Print out the missing names
# print(missing_names)
all(colnames(countData) == rownames(colData))
## TRUE

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

#ids = list of IDS
#fromKey = key type; toKey = key type we want to convert to
#db = the AnnotationDb object to use.
#ifMultiple = the argument specifies what to do if one source ID maps to several target IDs:
#should the function return an NA or simply the first of the multiple IDs?
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
head(resLFC, 4)
# Save to CSV without full path (if working directory is set)
write.csv(resLFC, file = "2024-06-13_deseq2.csv", row.names = FALSE)

# #### To make the MA plot
library("ggplot2")
remove.packages("ggrepel")
install.packages("ggrepel")
library(ggrepel)
library("dplyr")
library(ggpubr)

# We need to reorder the factor so that "No DEG" is plotted first (at the back).
resLFC$diffexpressed <- factor(resLFC$diffexpressed, levels = c("NO", "UP", "DOWN"))

# Define colors for different groups of genes
cols <- c("NO" = "light grey", "UP" = "orange", "DOWN" = "blue")

# Calculate log2 of baseMean if not already done
resLFC <- resLFC %>%
  mutate(log2_baseMean = log2(baseMean))

# Sort data by `diffexpressed` status and log2FoldChange
resLFC <- resLFC %>%
  arrange(diffexpressed, -log2FoldChange)

# Select the top 10 up-regulated and down-regulated genes
top_up_genes <- resLFC %>%
  filter(diffexpressed == "UP") %>%
  slice_max(order_by = log2FoldChange, n = 10) %>%
  pull(hgnc_symbol)
top_up_genes <- c("Hbb-y","Tead3","Myh14", "Sfta3-ps", "Lrig3", "Cast", "Phldb2", "Lhx8")

top_down_genes <- resLFC %>%
  filter(diffexpressed == "DOWN") %>%
  slice_min(order_by = log2FoldChange, n = 10) %>%
  pull(hgnc_symbol)
top_down_genes <- c("Xlr3b","Cpne7", "Atp5me", "Nptxr", "Cdk5r2", "Adap1", "L3mbtl1")
# Combine genes to label
all_genes_to_label <- unique(c(top_up_genes, top_down_genes, "Tead2","Cdk4"))

# Subset data frame to just the genes you want to label
label_data <- resLFC %>%
  filter(hgnc_symbol %in% all_genes_to_label)

# Use the modified data to plot the MA plot
resLFC$diffexpressed <- factor(resLFC$diffexpressed, levels = c("NO", "DOWN", "UP"))

# Assuming `cols` contains the color definitions for each `diffexpressed` status
cols <- c("NO" = "lightgrey", "DOWN" = "blue", "UP" = "dark orange")

# Use the modified data to plot the MA plot
# Sort the data by `diffexpressed` column
resLFC <- resLFC[order(resLFC$diffexpressed), ]

# Create the MA plot using the diffexpressed column for colors
ma_plot <- ggplot() +
  # Plot non-differentially expressed genes first
  geom_point(data = subset(resLFC, diffexpressed == "NO"),
             aes(x = log2(baseMean), y = log2FoldChange, color = diffexpressed), alpha = 1) +
  # Plot down-regulated genes second
  geom_point(data = subset(resLFC, diffexpressed == "DOWN"),
             aes(x = log2(baseMean), y = log2FoldChange, color = diffexpressed), alpha = 3) +
  # Plot up-regulated genes last
  geom_point(data = subset(resLFC, diffexpressed == "UP"),
             aes(x = log2(baseMean), y = log2FoldChange, color = diffexpressed), alpha = 3) +
  # Add gene labels with ggrepel
  geom_text_repel(data = label_data, aes(x = log2_baseMean, y = log2FoldChange, label = hgnc_symbol),
                  size = 4, color = "black", 
                  arrow = arrow(length = unit(0.01, "npc")),
                  point.padding = 0.5, segment.color = "black",
                  force = 2.0, max.overlaps = Inf) + 
  xlim(0, max(log2(resLFC$baseMean)) + 1) +
  ylim(-2, 2) +
  
  # Add scale_color_manual to handle colors for the legend
  scale_color_manual(values = c("UP" = "dark orange", "DOWN" = "blue", "NO" = "light grey"),
                     labels = c("UP" = "Up-Regulated", "DOWN" = "Down-Regulated", "NO" = "Not Differentially Expressed"),
                     name = "Expression Status") +
  
  theme_classic() +
  labs(title = "MA Plot", x = "Log2 Mean Expression", y = "Log2 Fold Change") +
  theme(legend.position = "right", 
        legend.title = element_blank())


# Plot the MA plot
print(ma_plot)

# Save the plot if you wish
 dev.copy(pdf, file="MAplot.pdf")
 dev.off()

 ggsave("MAplot.pdf", ma_plot, width = 17.73, height = 17.98, units = "cm")


genes_per_DEG <- table(resLFC$diffexpressed)
print(genes_per_DEG)
#DOWN    NO    UP 
  247 53247   146 

 
### Save the document
# Export file 
write.csv(as.data.frame(resLFC), file = "2024-03-13-deseq_CE16_CKO_WT.csv")




