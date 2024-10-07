### With files generated after deseq2_analysis

### Make a new data frame for both up and down DEG
up_resLFC<- resLFC[resLFC$diffexpressed == "UP",]
down_resLFC <- resLFC[resLFC$diffexpressed == "DOWN",]

########## Gene ontology 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler", force = TRUE)
library(clusterProfiler)
library(enrichplot)
if (!("msigdbr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("msigdbr", update = FALSE)
}
library(msigdbr)

### To add the pathways and hallmarks of associated genes. Will be nedeed to run GSEA function
mm_hallmark_sets <- msigdbr(
  species = "Mus musculus", # Replace with species name relevant to your data
  category = "H"
)
# We will need this so we can use the pipe: %>%
library(magrittr)

gene_list <- as.character(up_resLFC$entrezid)

# To get the universe of all genes, you can use the org.Mm.eg.db package
all_genes <- keys(org.Mm.eg.db, keytype = "ENTREZID")

go_results <- enrichGO(
  gene         = gene_list,
  universe     = all_genes,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "ALL", # "BP", "MF", "CC", or "ALL"
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

#Summarize and view the results:
summary(go_results)

#Visualize the results:
#Dotplot
dotplot_up <- dotplot(go_results, showCategory=5)
dotplot_up <- dotplot_up + theme(axis.title.y = element_text(size = 6))
print(dotplot_up)
ggsave("GO_UP_dotplot.pdf", dotplot_up, width = 15.98, height = 20.73, units = "cm")

#### DOWN

down_gene_list <- as.character(down_resLFC$entrezid)

go_results_down <- enrichGO(
  gene         = down_gene_list,
  universe     = all_genes,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "ALL", # "BP", "MF", "CC", or "ALL"
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

#Summarize and view the results:
summary(go_results_down)

#Visualize the results:
#Dotplot
dotplot_down <- dotplot(go_results_down, showCategory=5)
dotplot_down <- dotplot_down + theme(axis.title.y = element_text(size = 6))
print(dotplot_down)
ggsave("GO_DOWN_dotplot.pdf", dotplot_down, width = 15.98, height = 20.73, units = "cm")

