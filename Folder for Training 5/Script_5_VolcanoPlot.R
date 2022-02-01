# Install libraries if necessary
library(ggplot2)
library(ggrepel)
library(dplyr)

####################################### Data loading and pre-processing #######################################
# Select the DEGs result files for comparison
File <- file.choose() #Human or mouse file
DEGs <- read.csv(File, header = TRUE)

res <- data.frame(DEGs)
res$Symbol <- res$X

##################################### Volcano Plot #######################################################
## Two ways to label desired genes (think about which one you prefer)

# 1. Set cutoff for log2FoldChange for label
FC_neg_cutoff = -1 # Genes with log2FoldChange less than this value will be labeled
FC_pos_cutoff = 1 # Genes with log2FoldChange more than this value will be labeled
title = "Training 5_Volcano Plot"

# 2. Label a certain sets of genes, change gene symbols based on need
geneList <-list("REG1B","MAL","PRAC1","FOXF2","NDN", "CSAG1","RTL8C","PNMA8A","RNU5B-1","MYH11",
                "IVL","LINC02582","CGB5","LNCAROD","KLK5","CALB1","KLK7","KLK6","CT83","ONECUT3")

# Make volcano plot
makeVolcanoPlot <- function(res,FC_neg_cutoff, FC_pos_cutoff, title, geneList){
  # Set threshold to have different color for different sets of genes
  res <- res %>% mutate(threshold = ifelse(padj<=0.050 & abs(log2FoldChange)>=1.5,"A", ifelse(padj>0.05 & abs(log2FoldChange)<1.5, "B", "C")))
  
  ## Two ways to label desired genes (run only one of them)
  # 1. Only label genes based on the cutoffs
  res$Label <- ifelse(res$log2FoldChange < FC_neg_cutoff & res$padj < 0.05 | res$log2FoldChange > FC_pos_cutoff & res$padj < 0.05, res$Symbol, "")
  
  # 2. Only label genes based on pre-defined gene list (Run either cutoff or pre-defined gene list)
  #res$Label <- ifelse(res$Symbol %in% geneList, res$Symbol, "")
  
  # Define plotting function, adjust parameters based on need
  ggplot(res, aes(log2FoldChange, -log10(padj), label = Label))+
    # Set color for different data points
    geom_point(aes(colour = threshold), size =2, alpha = 0.4) +
    scale_colour_manual(values = c("A"="red", "B"="grey", "C"="darkgreen"))+
    # Avoid text overlap
    geom_text_repel(box.padding = 0.5, max.overlaps = Inf)+
    # Set theme of the plot
    theme_bw()+
    xlim(c(-5,5))+ # x-axis range, adjust based on need
    ylim(c(0,5))+ # y-axis range,adjust based on need
    geom_vline(xintercept = c(-1.5, 1.5), lty=2,col="black",lwd=0.6)+ # Draw horizontal dash line
    geom_hline(yintercept = 1.301, lty=2,col="black",lwd=0.6)+ # Draw vertical dash line
    xlab(bquote(~Log[2]~ 'fold change'))+ # x-axis label
    ylab(bquote(~-Log[10]~ 'padj')) + # y-axis label
    theme(axis.text = element_text(color = "black", size = 10))+ # Font type and size
    ggtitle(title) # Add title
}

makeVolcanoPlot(res,FC_neg_cutoff, FC_pos_cutoff, title, geneList)

#Clean up variables and release RAM before loading next sets of data
rm(list=ls())
gc()


####################################### Quick Preview ##########################################
# Preview of the volcano plot, helps define some parameters quickly
# Make volcano plot (https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html)
#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
EnhancedVolcano(res, lab = res$Symbol,
                x = 'log2FoldChange', y = 'padj',
                pCutoff = 0.05, FCcutoff = 1.5, 
                pointSize =3, labSize = 4)

################################## Ensembl Gene IDs to Gene Symbols #######################################
### Sometimes you may have Ensembl Gene IDs instead of Gene Symbols
### Run the following cov_ID function to annotate the Ensembl Gene IDs to Gene Symbols
library(org.Hs.eg.db) # For Human sample, change based on need

cov_ID <- function(DEGs){
  df <- data.frame(DEGs)
  ens <- df$X #X is the column name for Ensembl Gene IDs
  # Covert Gene Nomenclature type
  symbols <- mapIds(org.Hs.eg.db, keys = ens, column = c('SYMBOL'), keytype = 'ENSEMBL')
  # Drop NA values and map the Gene Symbols to corresponding genes
  symbols <- symbols[!is.na(symbols)]
  symbols <- symbols[match(df$X, names(symbols))]
  df$Symbol <- symbols
  
  # Drop NA Gene Symbols
  df <- df[!is.na(df$Symbol),]
  # Drop duplicate Gene Symbols
  df <- df[!duplicated(df$Symbol),]
  
  return(df)
}

res <- cov_ID(DEGs)

# Save the result if you will
#write.csv(res,basename(File))

