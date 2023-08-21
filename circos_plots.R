#CIRCOS PLOTS
install.packages("circlize")
library(circlize)
library(dplyr)
library(ggsci)
library(RColorBrewer)
setwd("/Users/rusiq/Desktop/IGR/")

SRP074349 = read.csv("SRP074349_large_scale_df.csv")
SRP162843 = read.csv("SRP162843_large_scale_df.csv")
SRP313282 = read.csv("SRP313282_large_scale_df.csv")
ERP001058 = read.csv("ERP001058_large_scale_df.csv")
SRP090460 = read.csv("SRP090460_large_scale_df.csv")
EMTAB6957 = read.csv("EMTAB6957_large_scale_df.csv")
SRP045225 = read.csv("SRP045225_large_scale_df.csv")

large_scale = rbind(SRP074349, SRP162843, SRP313282, SRP045225,
                    ERP001058, SRP090460, EMTAB6957)

rownames(large_scale) = large_scale$X
combined_deletions = read.csv("deletions_combined_dataframe.csv")
selected_deletions = large_scale[rownames(large_scale) %in% combined_deletions[combined_deletions$control == 0, "sample",],]


arms_list <- strsplit(selected_deletions$LargeScaleDel, " ")

cooc_matrix <- matrix(0, nrow = length(layout$name), ncol = length(layout$name), 
                      dimnames = list(layout$name, layout$name))

for (arms in arms_list) {
  for (i in seq_along(arms)) {
    for (j in seq_along(arms)) {
      cooc_matrix[arms[i], arms[j]] <- cooc_matrix[arms[i], arms[j]] + 1
    }
  }
}

mat = cooc_matrix[,c("3p", "9p", "17p")]

all_arms = c(paste0(1:23, "p"), paste0(1:23, "q"))
grid.col = rep("grey", length(all_arms))
names(grid.col) = all_arms
grid.col["3p"] <- "red"
grid.col["9p"] <- "green"
grid.col["17p"] <- "blue"

chordDiagram(t(mat), grid.col = grid.col)

 
#GENERAL DELETIONS STATISTICS
arms_list <- strsplit(selected_deletions$LargeScaleDel, " ")
arm_counts <- table(unlist(arms_list))
arm_counts <- arm_counts[order(arm_counts, decreasing = TRUE)]

arm_counts_df <- data.frame(ChromosomalArm = names(arm_counts),
                            Count = as.numeric(arm_counts),
                            Count_percents = (as.numeric(arm_counts)/dim(selected_deletions)[1])*100)

custom_colors <- sample(colorRampPalette(brewer.pal(8, "Set2"))(40))

ggplot(arm_counts_df, aes(x = reorder(ChromosomalArm, -Count_percents), y = Count_percents, fill = ChromosomalArm)) +
  geom_bar(stat = "identity") +
  labs(x = "Chromosomal Arm", y = "Number of deletions") +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_bw()


