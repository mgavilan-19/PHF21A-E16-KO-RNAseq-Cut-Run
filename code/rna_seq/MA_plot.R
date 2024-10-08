#### Using the files generated in Deseq2

# #### To make the MA plot
library("ggplot2")
install.packages("ggrepel")
library(ggrepel)
library("dplyr")
library(ggpubr)

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
