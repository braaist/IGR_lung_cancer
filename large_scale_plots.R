install.packages("circlize")
library(circlize)

install.packages("svglite")
library(svglite)
source("/Users/rusiq/Desktop/IGR/plotLargeScaleEvent_custom.R")

#SEPARATED

test = readRDS("SRP074349_finals.rds")
obj <- test[[9]]
plot <- plotLargeScaleEvent_custom(object=obj)
plot
write.csv(data.frame(obj@large.scale.cnv.events), "SRP074349_large_scale_df.csv")

test = readRDS("SRP162843_finals.rds")
obj <- test[[9]]
plot <- plotLargeScaleEvent_custom(object=obj)
plot
write.csv(data.frame(obj@large.scale.cnv.events), "SRP162843_large_scale_df.csv")

test = readRDS("SRP313282_finals.rds")
obj <- test[[9]]
plot <- plotLargeScaleEvent_custom(object=obj)
plot
write.csv(data.frame(obj@large.scale.cnv.events), "SRP313282_large_scale_df.csv")

test = readRDS("SRP090460_finals.rds")
obj <- test[[9]]
plot <- plotLargeScaleEvent_custom(object=obj)
plot
write.csv(data.frame(obj@large.scale.cnv.events), "SRP090460_large_scale_df.csv")

test = readRDS("ERP001058_finals.rds")
obj <- test[[9]]
plot <- plotLargeScaleEvent_custom(object=obj)
plot
write.csv(data.frame(obj@large.scale.cnv.events), "ERP001058_large_scale_df.csv")

test = readRDS("EMTAB6957_finals.rds")
obj <- test[[9]]
plot <- plotLargeScaleEvent_custom(object=obj)
plot
write.csv(data.frame(obj@large.scale.cnv.events), "EMTAB6957_large_scale_df.csv")

test = readRDS("SRP045225_finals.rds")
obj <- test[[9]]
plot <- plotLargeScaleEvent_custom(object=obj)
plot
write.csv(data.frame(obj@large.scale.cnv.events), "SRP045225_large_scale_df.csv")


#MERGE DELETIONS WITH SELECTED SAMPLES ONLY
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

p <- plotLargeScaleEvent_samps(selected_deletions)
p + ggtitle("Selected deletions samples")


selected_deletions$LargeScaleDel 
