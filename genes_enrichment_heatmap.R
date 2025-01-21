library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(reshape2)

#Read filtered gene lists
txi_3p_res = read.csv("/Users/rusiq/Desktop/IGR/txi_3p_res_output_filtered.csv")
txi_9p_res = read.csv("/Users/rusiq/Desktop/IGR/txi_9p_res_output_filtered.csv")
txi_17p_res = read.csv("/Users/rusiq/Desktop/IGR/txi_17p_res_output_filtered.csv")

#Read enrichment results
vec_3p_unique = c("positive regulation of lymphocyte activation", "T cell mediated immunity")
vec_9p_unique = c("glutamatergic synapse", "postsynaptic specialization", "thyroid hormone metabolic process", "celular modified amino acid metabolic process")
vec_17p_unique = c("dendritic tree", "postsynaptic specialization")
pathways_vec = c(vec_3p_unique, vec_9p_unique, vec_17p_unique)

gse_selected = read.csv("compare_cluster_result_total.csv")
gse_selected = gse_selected[gse_selected$Description %in% pathways_vec,]
View(gse_selected)

split_genes <- function(df) {
  df %>%
    mutate(geneID = strsplit(geneID, "/")) %>%
    unnest(geneID)
}

# Split genes in gse_selected
gse_long <- split_genes(gse_selected)

# Function to merge log2FoldChange data
merge_log2fc <- function(genes, res_df, suffix) {
  genes %>%
    left_join(res_df, by = c("geneID" = "gene_name")) %>%
    rename_with(~paste0("log2FoldChange_", suffix), .cols = log2FoldChange)
}

# Merge log2FoldChange data for each cluster and create final result
result <- gse_long %>%
  merge_log2fc(txi_3p_res, "3p") %>%
  merge_log2fc(txi_9p_res, "9p") %>%
  merge_log2fc(txi_17p_res, "17p") %>%
  group_by(geneID) %>%
  summarise(
    log2FoldChange_3p = first(log2FoldChange_3p),
    log2FoldChange_9p = first(log2FoldChange_9p),
    log2FoldChange_17p = first(log2FoldChange_17p),
    Cluster = first(Cluster),
    Description = paste(unique(Description), collapse = ", ")
  ) %>%
  ungroup() %>%
  arrange(Cluster, Description)

# Create a sorting key for clusters
View(result)




result$log2FoldChange_3p <- as.numeric(result$log2FoldChange_3p)
result$log2FoldChange_9p <- as.numeric(result$log2FoldChange_9p)
result$log2FoldChange_17p <- as.numeric(result$log2FoldChange_17p)

# Convert dataframe to long format
df_long <- melt(result[, c("geneID", "log2FoldChange_3p", "log2FoldChange_9p", "log2FoldChange_17p", "Description")],
                id.vars = c("geneID", "Description"),
                variable.name = "Condition",
                value.name = "log2FoldChange")

# Plot heatmap
ggplot(df_long, aes(x = Condition, y = geneID, fill = log2FoldChange)) +
  geom_tile(color = "black", size = 0.3) +  # Tiny black borders
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", 
                       midpoint = 0, name = "Log2 Fold Change") +
  labs(x = "Condition", y = "Gene ID") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0)) +  # Rotate x-axis labels for readability
  facet_grid(rows = vars(Description), scales = "free_y", space = "free_y") +
  ggtitle("Unique genes pathway enrichment heatmap")




