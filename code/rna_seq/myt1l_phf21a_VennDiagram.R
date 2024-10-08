# Read the Excel sheet into a data frame
weigel_2023 <- read.table(file = "Myt1l_mouse_bulk_deseq2.txt",sep="\t", header=TRUE)
colnames(weigel_2023)

## Generate files
# E18.5 Myt1l -/- data
hom_E18 <- weigel_2023[,c(1,2,7:10)]
# P0 Myt1l -/- data
hom_P0 <- weigel_2023[,c(1,2,15:18)]
# AdMyt1l +/- data
ad_het <- weigel_2023[,c(1,2,23:26)]

####### Continuing analysis from here.....
new_names <- c("gene_name","ENSEMBL", "baseMean","log2fc_shrunk","pvalue","p_adj")
colnames(hom_E18) <- new_names

Phf21a_E16 <- read.csv(file = "2024-03-13-deseq_CE16_CKO_WT.csv")
new_names <- c("X.1","baseMean","log2fc_shrunk","lfcSE","pvalue","p_adj","diffexpressed", "delabel", "X", 
            "gene_name", "entrezid")
colnames(Phf21a_E16) <- new_names

### Select Myt1l DEG based on our analysis
# Added column to classify up or down regulated genes
hom_E18$diffexpressed <- "NO"
# If log2FoldChange > 0 and pvalue < 0.05, set as "UP"
hom_E18$diffexpressed[hom_E18$log2fc_shrunk > 0 & hom_E18$p_adj < 0.05] <- "UP"
# If log2FoldChange < - 0 and pvalue < 0.05, set as "DOWN"
hom_E18$diffexpressed[hom_E18$log2fc_shrunk < -0 & hom_E18$p_adj < 0.05] <- "DOWN"
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
hom_E18$gene_id <- NA
hom_E18$ENSEMBL <- rownames(hom_E18)
hom_E18$gene_id[hom_E18$diffexpressed != "NO"] <- hom_E18$ENSEMBL[hom_E18$diffexpressed != "NO"]

# Verify that there are no duplicates
print(any(duplicated(hom_E18$gene_name))) # Should return FALSE
#TRUE
print(any(duplicated(Phf21a_E16$gene_name))) # Should return FALSE 
#TRUE

# Remove duplicate rows by gene_name
hom_E18_unique <- hom_E18 %>% distinct(gene_name, .keep_all = TRUE)
print(any(duplicated(hom_E18_unique$gene_name)))
#FALSE
Phf21a_E16_unique <- Phf21a_E16 %>% distinct(gene_name, .keep_all = TRUE)
print(any(duplicated(Phf21a_E16_unique$gene_name)))
#FALSE

###### Generate vendiagram
#Find common and unique DEGs in weigel 2023:
all_Up_weigel <- hom_E18_unique[hom_E18_unique$diffexpressed == "UP",]
all_Down_weigel <- hom_E18_unique[hom_E18_unique$diffexpressed == "DOWN",]
# Merging in one data frame 
all_DEGs_weigel <- rbind(all_Up_weigel, all_Down_weigel)

#Find common and unique DEGs in our data:
all_Up <- Phf21a_E16_unique[Phf21a_E16_unique$diffexpressed == "UP",]
all_Down <- Phf21a_E16_unique[Phf21a_E16_unique$diffexpressed == "DOWN",]
# Merging
all_DEGs <- rbind(all_Up, all_Down)
View(all_DEGs)

# Extract gene names
genes_all_DEGs_weigel <- all_DEGs_weigel$gene_name
genes_all_DEGs <- all_DEGs$gene_name

# Find common and unique DEGs
common_genes <- intersect(genes_all_DEGs_weigel, genes_all_DEGs)
unique_genes_all_DEGs_weigel<- setdiff(genes_all_DEGs_weigel, genes_all_DEGs)
unique_genes_all_DEGs <- setdiff(genes_all_DEGs,genes_all_DEGs_weigel)

# Total number of unique DEGs
total_DEGs <- length(unique(c(genes_all_DEGs_weigel, genes_all_DEGs)))
total_DEGs #1779
# Contingency table calculation
myt1l_count <- length(genes_all_DEGs_weigel)
myt1l_count #1544
phf21a_count <- length(genes_all_DEGs)
phf21a_count #378
overlap_count <- length(common_genes)
overlap_count #143

only_myt1l_count <- myt1l_count - overlap_count
only_phf21a_count <- phf21a_count - overlap_count
neither_count <- total_DEGs - (myt1l_count + phf21a_count - overlap_count)

# Create the contingency table
contingency_table <- matrix(
  c(overlap_count, overlap_count, only_phf21a_count, only_myt1l_count),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(Mutation = c("Overlap", "No Overlap"), Genes = c("Phf21a","Myt1l"))
)
print(contingency_table)

# Perform Fisher's exact test
fisher_result <- fisher.test(contingency_table)
print(fisher_result)
#p-value < 2.2e-16 

# Install and load the VennDiagram package if not already done
library(ggplot2)
install.packages("ggVennDiagram")
library(ggVennDiagram)

# Create Venn diagram
venn.plot <- ggVennDiagram(
  x=list(Myt1l = genes_all_DEGs_weigel, Phf21a = genes_all_DEGs),
  category.names = c("Myt1l", "Phf21a"),
  label_alpha = 0,
  label_size = 10,
  set_color = c("dark blue","dark green")
  ) 
library(scales)

a <- venn.plot + 
  scale_fill_gradientn(colors = c(
    alpha("orange", 0.5),  # 60% transparent grey
    alpha("dark green", 0.6),  # 60% transparent dark green
    alpha("dark blue", 0.6)  # 60% transparent dark blue
  ), values = c(0, 0.5, 0.9, 1))  # Adjust these values to control gradient

a
# Save the plot as a PDF
 dev.copy(pdf, file="E18_Myt1L_null_PHF21A_mut.pdf")
 dev.off()

 ggsave("E18_Myt1L_null_PHF21A_mut.pdf", a , width = 17.73, height = 17.98, units = "cm")
