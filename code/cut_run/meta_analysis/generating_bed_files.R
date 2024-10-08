###### Starting from the files generating after Deseq2 analysis

##### Generating BED files
### Get gene length from gtf file 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rtracklayer", force = TRUE)

gtf <- rtracklayer::import('gencode.vM21.annotation.gtf')
gtf_df=as.data.frame(gtf)
head(gtf_df)

# Before merging gotta remove duplicates and keep only the first occurrence
gtf_df <- gtf_df[!duplicated(gtf_df$gene_id), ]

### Make a new data frame for both up and down DEG
up_resLFC<- resLFC[resLFC$diffexpressed == "UP",]
down_resLFC <- resLFC[resLFC$diffexpressed == "DOWN",]

# Merge gene length information with the main data based on GeneID
# To do so is necessary to clarify only to merge on the ENSMUG common numbers=
# Extract common part of gene_id before the last dot
up_resLFC$common_gene_id <- sub("\\.[0-9]+$", "", up_resLFC$gene_id)
gtf_df$common_gene_id <- sub("\\.[0-9]+$", "", gtf_df$gene_id)
# Merge data frames based on the common_gene_id
up_merged_data <- merge(up_resLFC, gtf_df, by = "common_gene_id", all.x = TRUE)
# Remove the additional column used for merging
up_merged_data$common_gene_id <- NULL
# DOWN DEG
down_resLFC$common_gene_id <- sub("\\.[0-9]+$", "", down_resLFC$gene_id)
# Merge data frames based on the common_gene_id
down_merged_data <- merge(down_resLFC, gtf_df, by = "common_gene_id", all.x = TRUE)
# Remove the additional column used for merging
down_merged_data$common_gene_id <- NULL


#### Make gtf file and then bed files ---> Example for the Up DEG.
# Remove unwanted columns
up_merged_data <- up_merged_data[ ,c(7,9,11,12,13,15) ]
#Save gtf file
write.table(up_merged_data , "2024-06-19_up_DEG.gtf", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

library(data.table)
#Read gtf file
up_DEG_gtf <- fread('2024-06-19_up_DEG.gtf', header = FALSE)
head(up_DEG_gtf)
new_names <- c("gene_id","seqnames","start","end","strand","gene_symbol")
colnames(up_DEG_gtf) <- new_names
up_DEG_gtf$score <- '0' 

##### To save it as a bed file
#Install and load necessary packages if not already installed
install.packages(c("rtracklayer", "GenomicRanges"))
library(rtracklayer)
library(GenomicRanges)
install.packages("BiocManager")
BiocManager::install("rtracklayer", version = "3.18", force = TRUE)
BiocManager::install("GenomicRanges", version = "3.18", force = TRUE)

gr_up_DEG <- makeGRangesFromDataFrame(up_DEG_gtf,
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
bed_df_up_DEG <- data.frame(
  chr = seqnames(gr_up_DEG),
  start = start(gr_up_DEG),
  end = end(gr_up_DEG),
  name = mcols(gr_up_DEG)$gene_id,
  score = score(gr_up_DEG),
  strand = strand(gr_up_DEG)
)

# Save the BED file
write.table(bed_df_up_DEG , "2024_06_19_UP_DEG.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


##### To make bed_file including the entire genome
gtf_df$score <- '0' # Adding the score 0 will allow generating the bed files 
gr_all <- makeGRangesFromDataFrame(gtf_df,
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
bed_df_all <- data.frame(
  chr = seqnames(gr_all),
  start = start(gr_all),
  end = end(gr_all),
  name = mcols(gr_all)$gene_id,
  score = score(gr_all),
  strand = strand(gr_all)
)

# Save the BED file
write.table(bed_df_all , "2024_09_25_all_genes.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



