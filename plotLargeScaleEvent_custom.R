#THE ORIGINAL VARIANT OF CASPER FUNCTION
#BUT RETURN PLOT, NOT SAVE IT
library(ggrepel)

plotLargeScaleEvent_custom <- function(object) {
  
  samps <- object@large.scale.cnv.events
  chrs <- as.vector(sapply(1:22, function(x) c(paste(x, "p", sep = ""), paste(x, "q", sep = ""))))
  chrMat <- matrix(0, ncol = 44, nrow = length(rownames(samps)))
  colnames(chrMat) <- chrs
  rownames(chrMat) <- rownames(samps)
  
  
  for (x in 1:dim(samps)[1]) {
    
    chrWithEvents <- as.vector(unlist(strsplit(as.vector(paste(samps$LargeScaleAmp[x], collapse = " ")), split = " ")))
    chrMat[x, match(intersect(chrWithEvents, chrs), colnames(chrMat))] <- 1
    
    chrWithEvents <- as.vector(unlist(strsplit(as.vector(paste(samps$LargeScaleDel[x], collapse = " ")), split = " ")))
    chrMat[x, match(intersect(chrWithEvents, chrs), colnames(chrMat))] <- (-1)
    
  }
  
  plot.data <- melt(chrMat)
  
  plot.data$value2 <- "neutral"
  plot.data$value2[plot.data$value > 0] <- "amplification"
  plot.data$value2[plot.data$value < 0] <- "deletion"
  plot.data$value2 <- factor(plot.data$value2, levels = c("amplification", "deletion", "neutral"))
  plot.data$X2 <- factor(plot.data$X2, levels = colnames(chrMat))
  
  p <- ggplot(aes(x = X2, y = X1, fill = value2), data = plot.data) + geom_tile(colour = "gray50", size = 0.01) + labs(x = "", 
                                                                                                                       y = "") + scale_fill_manual(values = c(amplification = muted("red"), deletion = muted("blue"), neutral = "white")) + 
    theme_grey(base_size = 6) + theme(legend.position = "right", legend.direction = "vertical", legend.title = element_blank(), 
                                      strip.text.x = element_blank(), legend.text = element_text(colour = "black", size = 7, face = "bold"), legend.key.height = grid::unit(0.8, 
                                                                                                                                                                            "cm"), legend.key.width = grid::unit(0.5, "cm"), axis.text.x = element_text(size = 5, colour = "black", angle = -45, 
                                                                                                                                                                                                                                                        hjust = 0), axis.text.y = element_text(size = 6, vjust = 0.2, colour = "black"), axis.ticks = element_line(size = 0.4), 
                                      plot.title = element_text(colour = "black", hjust = 0, size = 6, face = "bold"))
  return(p)
  
}


#THE ORIGINAL VARIANT OF CASPER FUNCTION
#BUT ALLOW CUSTOM SAMPS

plotLargeScaleEvent_samps <- function(samps) {
  
  chrs <- as.vector(sapply(1:22, function(x) c(paste(x, "p", sep = ""), paste(x, "q", sep = ""))))
  chrMat <- matrix(0, ncol = 44, nrow = length(rownames(samps)))
  colnames(chrMat) <- chrs
  rownames(chrMat) <- rownames(samps)
  
  
  for (x in 1:dim(samps)[1]) {
    
    chrWithEvents <- as.vector(unlist(strsplit(as.vector(paste(samps$LargeScaleAmp[x], collapse = " ")), split = " ")))
    chrMat[x, match(intersect(chrWithEvents, chrs), colnames(chrMat))] <- 1
    
    chrWithEvents <- as.vector(unlist(strsplit(as.vector(paste(samps$LargeScaleDel[x], collapse = " ")), split = " ")))
    chrMat[x, match(intersect(chrWithEvents, chrs), colnames(chrMat))] <- (-1)
    
  }
  
  plot.data <- melt(chrMat)
  
  plot.data$value2 <- "neutral"
  plot.data$value2[plot.data$value > 0] <- "amplification"
  plot.data$value2[plot.data$value < 0] <- "deletion"
  plot.data$value2 <- factor(plot.data$value2, levels = c("amplification", "deletion", "neutral"))
  plot.data$X2 <- factor(plot.data$X2, levels = colnames(chrMat))
  
  p <- ggplot(aes(x = X2, y = X1, fill = value2), data = plot.data) + geom_tile(colour = "gray50", size = 0.01) + labs(x = "", 
                                                                                                                       y = "") + scale_fill_manual(values = c(amplification = muted("red"), deletion = muted("blue"), neutral = "white")) + 
    theme_grey(base_size = 6) + theme(legend.position = "right", legend.direction = "vertical", legend.title = element_blank(), 
                                      strip.text.x = element_blank(), legend.text = element_text(colour = "black", size = 7, face = "bold"), legend.key.height = grid::unit(0.8, 
                                                                                                                                                                            "cm"), legend.key.width = grid::unit(0.5, "cm"), axis.text.x = element_text(size = 5, colour = "black", angle = -45, 
                                                                                                                                                                                                                                                        hjust = 0), axis.text.y = element_text(size = 6, vjust = 0.2, colour = "black"), axis.ticks = element_line(size = 0.4), 
                                      plot.title = element_text(colour = "black", hjust = 0, size = 6, face = "bold"))
  return(p)
  
}


plotPCA.custom <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  
  ## Select the PCAs and percentVar that you like instead of 1 and 2
  d <- data.frame(PC1 = pca$x[, 3], PC2 = pca$x[, 4], group = group, 
                  intgroup.df, name = colData(dds_norm)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
    geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + coord_fixed()
  
}


plotPCA.custom(dds_norm, intgroup = "class")
