########################### GO Over-representation Analysis ########################################
library(clusterProfiler)
library(org.Hs.eg.db)

File <- file.choose() # Choose DEGs result file
res <-read.csv(File, header=TRUE)

# Define a function to do the analysis repeatedly
GO_ORA <- function(res) {
  #### Prepare Input for analysis
  # Reading in input from deseq2
  df <- data.frame(res)
  df$ensembl_gene_id <- rownames(res)
  
  # get log2 fold change data
  original_gene_list <- df$log2FoldChange
  
  # name the vector
  names(original_gene_list) <- df$ensembl_gene_id
  
  # omit any NA values
  gene_list <- na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  # Extract significant results (padj < 0.05)
  sig_genes_df = subset(df, padj < 0.05)
  
  # Extract results with baseMean > 10
  sig_genes_df = subset(sig_genes_df, baseMean > 10)
  
  # From significant results, we want to filter on log2fold change
  genes <- sig_genes_df$log2FoldChange
  
  # Name the vector
  names(genes) <- sig_genes_df$ensembl_gene_id
  
  # omit NA values
  genes <- na.omit(genes)
  
  # filter on minimum log2fold change (log2FoldChange > 1.5)
  genes <- names(genes)[abs(genes) > 1.5]
  
  #### GO Over-representation Analysis (adjust parameter based on need)
  go_enrich_res <- enrichGO(gene = genes,
                        universe = names(gene_list),
                        OrgDb = org.Hs.eg.db, # For human analysis
                        keyType = 'ENTREZID', #can be ENSEMBL, SYMBOL, or ENTREZID
                        readable = TRUE,
                        ont = "BP", #MF, BP, or CC
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
  
  return(go_enrich_res)
}

# Perform GO ORA based on need
go_enrich_res <- GO_ORA(res)


# Save the results to a csv file
save_res <- data.frame(go_enrich_res)
write.csv(save_res, "Training6_GO_Result.csv")
