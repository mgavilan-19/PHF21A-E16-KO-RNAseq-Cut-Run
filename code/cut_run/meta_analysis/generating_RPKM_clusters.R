# Read a text file into a data frame --- deseq coutns
data <- read.table("ctx_e16_count.txt", header = TRUE, sep = "\t")
# Change column name to merge with the other file later on
colnames(data)[colnames(data) == "Geneid"] <- "gene_id"

### Get gene length from gtf file 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rtracklayer", force = TRUE)

gtf <- rtracklayer::import('gencode.vM21.annotation.gtf')
gtf_df=as.data.frame(gtf)
head(gtf_df)

#Change width = gene length to Kbp and add to the data frame
gtf_df$width_kb <- gtf_df$width / 1000

# Before merging gotta remove duplicates and keep only the first occurrence
gtf_df <- gtf_df[!duplicated(gtf_df$gene_id), ]

# Merge gene length information with the main data based on GeneID
# To do so is necessary to clarify only to merge on the ENSMUG common numbers=
# Extract common part of gene_id before the last dot
data$common_gene_id <- sub("\\.[0-9]+$", "", data$gene_id)
gtf_df$common_gene_id <- sub("\\.[0-9]+$", "", gtf_df$gene_id)

# Merge data frames based on the common_gene_id
merged_data <- merge(data, gtf_df, by = "common_gene_id", all.x = TRUE)

# Remove the additional column used for merging
merged_data$common_gene_id <- NULL

###################### Calculate RPKM  ############################################

# Calculate the total reads for each condition
merged_data$total_reads_CWT_NN <- rowSums(merged_data[, c("CWT_NN_1.bam", "CWT_NN_2.bam", "CWT_NN_3.bam")])
merged_data$total_reads_CKO_NN <- rowSums(merged_data[, c("CKO_NN_1.bam", "CKO_NN_2.bam", "CKO_NN_3.bam")])

################### Calculate total reads for each condition / 1,000,000 = RPM
merged_data$RPKM_CWT_NN <- merged_data$total_reads_CWT_NN / 1000000
merged_data$RPKM_CKO_NN <- merged_data$total_reads_CKO_NN / 1000000

############### RPKM = RPM/ gene lenght in kilobase pair
merged_data$RPKM_CWT_NN <- (merged_data[, "total_reads_CWT_NN"] / merged_data[, "width_kb"])
merged_data$RPKM_CKO_NN <- (merged_data[, "total_reads_CKO_NN"] / merged_data[, "width_kb"])

### Create a new data frame for simplycity
# Extracting columns 
new_df <- merged_data[,c('gene_id.x','start', 'end', 'gene_type','gene_name', 'width_kb', 'RPKM_CWT_NN', 'RPKM_CKO_NN')]

## Calculate log2 of RPKM values
new_df$log2_RPKM_CWT_NN <- log2(new_df$RPKM_CWT_NN)
new_df$log2_RPKM_CKO_NN <- log2(new_df$RPKM_CKO_NN)

## Adding a pseudo-count of 1 to avoid taking a log of zero
new_df$log2_RPKM_CWT_NN_1 <- new_df$log2_RPKM_CWT_NN + 1
# Check for NA, NaN, or Inf values and remove rows that contain them
clean_df <- new_df[!is.infinite(new_df$log2_RPKM_CWT_NN_1) & !is.nan(new_df$log2_RPKM_CWT_NN_1) & !is.na(new_df$log2_RPKM_CWT_NN_1), ]

# Carry out k-means clustering
set.seed(123) # Setting seed for reproducibility
# Now perform k-means clustering on the cleaned data
kmeans_result <- kmeans(clean_df$log2_RPKM_CWT_NN_1, centers = 4, nstart = 25)

# Add the cluster assignment to your original data frame
clean_df$cluster <- kmeans_result$cluster

# Merge data frames based on the common_gene_id
clean_merged_data <- merge(clean_df, gtf_df, by = "gene_name", all.x = TRUE)

cluster_assignments <- kmeans_result$cluster
gene_cluster_df <- data.frame(
  Gene = rownames(clean_df),  # Assuming gene names are row names
  Cluster = cluster_assignments
)

genes_per_cluster <- table(gene_cluster_df$Cluster)
print(genes_per_cluster)
#1     2     3     4 
#9257 11457  5072  9226  


for (cluster in names(average_per_cluster)) {
  cat("Cluster", cluster, ": ", average_per_cluster[[cluster]], " ± ", sd_per_cluster[[cluster]], "\n")
}
#Cluster 1 :  9.001674  ±  1.29724 
#Cluster 2 :  6.041491  ±  0.9158448 
#Cluster 3 :  -1.44136  ±  1.546656 
#Cluster 4 :  2.616268  ±  1.06956 

###### make a data frame for each cluster
df_1 <- subset(clean_merged_data , cluster == '1')
df_2 <- subset(clean_merged_data , cluster == '2')
#switch so it goes in order
df_4 <- subset(clean_merged_data , cluster == '3')
df_3 <- subset(clean_merged_data , cluster == '4')


##### Make gtf file and then bed files for each data_frame. Here only showing for cluster 4
# Remove unwanted columns
df_4 <- df_4[ ,-c(1,5:12,14:16,18:48) ]
#Save gtf file
write.table(df_4 , "cluster_4_genes.gtf", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

library(data.table)
#Read gtf file
gtf_cluster_4 <- fread('cluster_4_genes.gtf', header = FALSE)
head(gtf_cluster_4)
new_names <- c("gene_id", "start", "end", "seqnames", "strand")
colnames(gtf_cluster_4) <- new_names
#gtf_cluster_3 <- gtf_cluster_3[,-c(8,11:26)]
gtf_cluster_4$score <- '0' 

# To save it as a bed file
#Install and load necessary packages if not already installed
install.packages("BiocManager")
BiocManager::install("rtracklayer", version = "3.18", force = TRUE)
BiocManager::install("GenomicRanges", version = "3.18", force = TRUE)
library(rtracklayer)
library(GenomicRanges)

gr_cluster_4 <- makeGRangesFromDataFrame(gtf_cluster_4,
                                         keep.extra.columns=TRUE,
                                         ignore.strand=FALSE,
                                         seqinfo=NULL,
                                         seqnames.field=c("seqnames", "seqname",
                                                          "chromosome", "chrom",
                                                          "chr", "chromosome_name",
                                                          "seqid"),
                                         start.field="start",
                                         end.field=c("end", "stop"),
                                         strand.field="strand",
                                         starts.in.df.are.0based=FALSE,
                                         na.rm=FALSE)

# Extract relevant information and create a BED file
bed_df_cluster_4 <- data.frame(
  chr = seqnames(gr_cluster_4),
  start = start(gr_cluster_4),
  end = end(gr_cluster_4),
  name = mcols(gr_cluster_4)$gene_id,
  score = score(gr_cluster_4),
  strand = strand(gr_cluster_4)
)

# Save the BED file
write.table(bed_df_cluster_4, "RPKM_gene_cluster_4_WT.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

