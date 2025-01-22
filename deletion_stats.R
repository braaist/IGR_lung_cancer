library(dplyr)
library(tidyr)
library(ggsci)
library(gplots)
library(ggplot2)
library(tidyverse)

#PLOT DELETION PROPORTIONS IN DIFFERENT GROUPS
deletions_dataframe = read.csv("deletions_combined_dataframe.csv")

deletions_dataframe$cancer_type = NA
deletions_dataframe[deletions_dataframe$name == "SRP045225","cancer_type"] <- "SCLC"
deletions_dataframe[deletions_dataframe$name != "SRP045225","cancer_type"] <- "NSCLC"

# Function to rename combinations
rename_combination <- function(b3p, b9p, b17p) {
  del <- c()
  if (b3p == -1) del <- c(del, "3p")
  if (b9p == -1) del <- c(del, "9p")
  if (b17p == -1) del <- c(del, "17p")
  if (length(del) == 0) return("No deletion")
  paste(del, collapse = ", ")
}

proportions_df <- deletions_dataframe %>%
  filter(control == 0) %>%
  group_by(binary_3p, binary_9p, binary_17p, name) %>%
  summarize(count = n(), .groups = "drop") %>%
  filter(!(binary_3p == 0 & binary_9p == 0 & binary_17p == 0)) %>%
  mutate(count_percent = count / sum(count) * 100) %>%
  rowwise() %>%
  mutate(new_combination = rename_combination(binary_3p, binary_9p, binary_17p)) %>%
  ungroup() %>%
  group_by(new_combination) %>%
  mutate(total_percent = sum(count_percent)) %>%
  ungroup() %>%
  mutate(new_combination = factor(new_combination, levels = unique(new_combination[order(total_percent, decreasing = TRUE)])))

ggplot(proportions_df, aes(x = new_combination, y = count_percent, fill = name)) +
  geom_bar(stat = "identity", colour = "black") +
  labs(x = "Deletions", y = "Percentage") +
  scale_fill_discrete(name = "Dataset") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  scale_fill_npg() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

proportions_df <- deletions_dataframe %>%
  mutate(cancer_type = if_else(name == "SRP045225", "SCLC", "NSCLC")) %>%
  filter(control == 0) %>%
  group_by(cancer_type, binary_3p, binary_9p, binary_17p) %>%
  summarize(count = n(), .groups = "drop") %>%
  mutate(count_percent = count / sum(count) * 100) %>%
  rowwise() %>%
  mutate(new_combination = rename_combination(binary_3p, binary_9p, binary_17p)) %>%
  ungroup() 

View(deletions_dataframe[deletions_dataframe$cancer_type == "SCLC",])
View(proportions_df)

ggplot(proportions_df, aes(x = cancer_type, y = count_percent, fill = new_combination)) +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  labs(x = "Cancer Type", y = "Proportion") +
  scale_fill_discrete(name = "Deletion Combination") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  scale_fill_aaas() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(proportions_df, aes(x = new_combination, y = count_percent, fill = cancer_type)) +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  labs(x = "Cancer Type", y = "Proportion") +
  scale_fill_discrete(name = "Deletion Combination") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  scale_fill_aaas() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Overall deletion stats
# Read and combine all CSV files
file_names <- c("SRP074349", "SRP162843", "SRP313282", "ERP001058", 
                "SRP090460", "EMTAB6957", "SRP045225")

large_scale <- map_dfr(file_names, ~read.csv(paste0(., "_large_scale_df.csv"), row.names = 1, stringsAsFactors = FALSE), .id = "source") %>%
  mutate(source = file_names[as.integer(source)])

# Process deletions
selected_deletions <- large_scale %>%
  filter(row.names(.) %in% combined_deletions[combined_deletions$control == 0, "sample"])

# Count arm deletions
arm_counts_df <- selected_deletions %>%
  pull(LargeScaleDel) %>%
  strsplit(" ") %>%
  unlist() %>%
  table() %>%
  as.data.frame() %>%
  setNames(c("ChromosomalArm", "Count")) %>%
  mutate(Count_percents = Count / nrow(selected_deletions) * 100) %>%
  arrange(desc(Count_percents)) %>%
  mutate(color = ifelse(ChromosomalArm %in% c("3p", "9p", "17p"), "#0fb500", "white"))

# Create plot
ggplot(arm_counts_df, aes(x = reorder(ChromosomalArm, -Count_percents), y = Count_percents, fill = color)) +
  geom_bar(stat = "identity", color = "black", size = 0.2) +
  labs(x = "Chromosomal Arm", y = "Percentage of deletions") +
  scale_fill_identity() +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

