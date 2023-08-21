#Fusion analysis
library(dplyr)
setwd("/Users/rusiq/Desktop/IGR/")

path_to_output = "/Users/rusiq/Desktop/IGR/star_fusion_results/"
fusion_list = list()
names_vec = c()
for (file in list.files(path_to_output)){
  fusions = read.csv(paste0(path_to_output, file), sep = "\t")
  fusion_list = append(fusion_list, list(fusions$X.FusionName))
  names_vec = c(names_vec, strsplit(strsplit(file, "\\.")[[1]][1], "_")[[1]][1])
}

names(fusion_list) = names_vec

#FUSION STATISTICS INTEGRATION
deletions_dataframe = read.csv("deletions_combined_dataframe.csv")
deletions_dataframe[deletions_dataframe$control == "1","class"] <- "control"
deletions_dataframe[deletions_dataframe$binary_3p == -1 & deletions_dataframe$binary_9p == -1  & deletions_dataframe$binary_17p == -1, "class"] <- "del3p_del9p_del17p"
View(deletions_dataframe[deletions_dataframe$binary_3p == -1 & deletions_dataframe$binary_9p == 0  & deletions_dataframe$binary_17p == 0, ])


merged_fusions_df <- data.frame(fusion = character(0), name = character(0), stringsAsFactors = FALSE)

for (name in names(fusion_list)) {
  sublist <- fusion_list[[name]]
  if (length(sublist) != 0) {
    merged_fusions_df <- rbind(merged_fusions_df, data.frame(fusion = unlist(sublist), name = name))
  }
}

merged_fusions_df = merged_fusions_df[merged_fusions_df$name %in% deletions_dataframe[deletions_dataframe$control == 0, "sample"],]
sample_class_df = deletions_dataframe[,c("sample", "class")]
merged_fusions_df = merge(merged_fusions_df, sample_class_df, by.x = "name", by.y = "sample", all.x = TRUE)

fusion_info <- merged_fusions_df %>%
  group_by(fusion) %>%
  summarize(total_fusions = n(), classes = list(class)) %>%
  
  ungroup()

class_counts <- merged_fusions_df %>%
  select(fusion, class) %>%
  unnest(class) %>%
  group_by(fusion, class) %>%
  summarize(total_fusions = n(), classes = list(class)) %>%
  group_by(fusion) %>%
  mutate(nSum = sum(total_fusions)) %>%
  filter(sum(total_fusions) > 20) %>%
  arrange(-nSum) %>%
  ungroup() 

order_counts = unique(class_counts$fusion)

ggplot(class_counts, aes(x = factor(fusion, level = order_counts), y = total_fusions, fill = class)) +
  geom_bar(stat = "identity", colour="black", size = 0.4) +
  labs(x = "Fusion type", y = "Total number of fusions") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), ) +
  guides(fill = guide_legend(title = "Class")) + scale_fill_npg() +
  ggtitle("Fusion types in cancer with different deletions")

fusion_info

test

