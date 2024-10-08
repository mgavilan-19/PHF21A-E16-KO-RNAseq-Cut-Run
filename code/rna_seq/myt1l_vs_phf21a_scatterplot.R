### Using the files generated in myt1l_vs_phf21a_VennDiagram.R
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

######### Scatterplot - transcriptomic correlation

# Merge the data frames by gene_name
combined_df <- inner_join(hom_E18_unique, Phf21a_E16_unique, by = "gene_name")

# Identify the unique classifications
combined_df <- combined_df %>%
  mutate(consistency = case_when(
    diffexpressed.x != "NO" & diffexpressed.y == "NO" ~ "Only_Myt1l",
    diffexpressed.x == "NO" & diffexpressed.y != "NO" ~ "Only_Phf21a",
     diffexpressed.x != "NO" & diffexpressed.y != "NO" ~ "Both", 
  ))

# Remove the other genes
combined_df_only <- combined_df[!is.na(combined_df$consistency),]

# Remove Rows with NA Values for correlation:
combined_df_clean <- combined_df %>% 
  filter(!is.na(log2fc_shrunk.y) & !is.na(log2fc_shrunk.x))

table(combined_df_clean$consistency)
# Both Only_Myt1l  Only_Phf21a 
# 143      1315         229 

correlation <- cor(combined_df_clean$log2fc_shrunk.y, combined_df_clean$log2fc_shrunk.x, method = "pearson")
correlation #0.33
R_squared <- correlation ^ 2
R_squared # 0.1134708

library(ggplot2)
# Plot the data
E16_vs_E18_null <- ggplot(combined_df_only, aes(x = log2fc_shrunk.y, y = log2fc_shrunk.x, color = consistency)) +
  geom_point(alpha = 5, size = 3) +
  labs(
    title = "Scatter plot of log2 fold changes",
    x = "Phf21a mutant Log2FoldChange ",
    y = "Myt1l +/- Log2FoldChange "
  ) +
  scale_color_manual(
    values = c(
      "Both" = "darkorange",
      "Only_Myt1l" = "#5252dc",
      "Only_Phf21a" = "darkgreen"
    )
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +  # Optional: Add regression line
  theme_minimal() +
  scale_x_continuous(limits = c(-1,1), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(limits = c(-1,1), expand = expansion(mult = c(0, 0.05))) +
  annotate(
    geom = "text",
    x = Inf,
    y = Inf,
    label = paste("R² = ", round(R_squared, 2), "\n", "Corr = ", round(correlation, 2)),
    hjust = 1.1, vjust = 1.1,
    size = 5,
    color = "black"
  )

print(E16_vs_E18_null)


# Save the plot if you wish
 dev.copy(pdf, file="E16phf_vs_E18myt1l_null.pdf")
 dev.off()

 ggsave("E16phf_vs_E18myt1l_null.pdf", E16_vs_E18_null , width = 17.73, height = 17.98, units = "cm")

#### Determine only phf21a DEG correlation with myt1l 
only_phf21a <- combined_df_only[c(combined_df_only$consistency == "Only_Phf21a"),]
# Plot the data
E16_vs_E18_null <- ggplot(only_phf21a, aes(x = log2fc_shrunk.y, y = log2fc_shrunk.x, color = consistency)) +
  geom_point(alpha = 5) +
  labs(
    title = "Scatter plot of log2 fold changes",
    x = "Phf21a mutant Log2FoldChange ",
    y = "Myt1l +/- Log2FoldChange "
  ) +
  scale_color_manual(
    values = c(
      "Only_Phf21a" = "darkgreen"
    )
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +  # Optional: Add regression line
  theme_minimal() +
  scale_x_continuous(limits = c(-1,1), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(limits = c(-1,1), expand = expansion(mult = c(0, 0.05))) +
  annotate(
    geom = "text",
    x = Inf,
    y = Inf,
    label = paste("R² = ", round(R_squared, 2), "\n", "Corr = ", round(correlation, 2)),
    hjust = 1.1, vjust = 1.1,
    size = 5,
    color = "black"
  )

print(E16_vs_E18_null)

# Save the plot if you wish
 dev.copy(pdf, file="E16phf_vs_E18myt1l_null_only_PHF21A.pdf")
 dev.off()

 ggsave("E16phf_vs_E18myt1l_null_only_PHF21A.pdf", E16_vs_E18_null , width = 17.73, height = 17.98, units = "cm")

#### Determine only myt1l DEG correlation with phf21a
only_myt1l <- combined_df_only[c(combined_df_only$consistency == "Only_Myt1l"),]
# Plot the data
E16_vs_E18_null <- ggplot(only_myt1l, aes(x = log2fc_shrunk.y, y = log2fc_shrunk.x, color = consistency)) +
  geom_point(alpha = 5) +
  labs(
    title = "Scatter plot of log2 fold changes",
    x = "Phf21a mutant Log2FoldChange ",
    y = "Myt1l +/- Log2FoldChange "
  ) +
  scale_color_manual(
    values = c(
      "Only_Myt1l" = "blue"
    )
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +  # Optional: Add regression line
  theme_minimal() +
  scale_x_continuous(limits = c(-1,1), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(limits = c(-1,1), expand = expansion(mult = c(0, 0.05))) +
  annotate(
    geom = "text",
    x = Inf,
    y = Inf,
    label = paste("R² = ", round(R_squared, 2), "\n", "Corr = ", round(correlation, 2)),
    hjust = 1.1, vjust = 1.1,
    size = 5,
    color = "black"
  )

print(E16_vs_E18_null)

# Save the plot if you wish
 dev.copy(pdf, file="E16phf_vs_E18myt1l_null_only_Myt1l.pdf")
 dev.off()

 ggsave("E16phf_vs_E18myt1l_null_only_Myt1l.pdf", E16_vs_E18_null , width = 17.73, height = 17.98, units = "cm")


