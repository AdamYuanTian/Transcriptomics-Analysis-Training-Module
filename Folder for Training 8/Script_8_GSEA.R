############################################ GSEA ##################################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(tidyr)
library(biomaRt)

File <- file.choose() # Choose DEGs result file
res <-read.csv(File, header=TRUE)

# GSEA analysis
# category defines the gene set (category: "H" for Hallmark, "C1" ... "C7"; subcategory: check website)
msigdb_GSEA <- function(res, category, subcategory) {
  # Add EntrezGene IDs based on Ensembl IDS (Needed for clusterProfiler GSEA)
  print("Converting Gene Symbol to EntrezGene IDs...")
  
  # Convert gene symbol to entrezgene IDs
  res_entrezID<-select(org.Hs.eg.db,
                  keys = res$X,
                  columns = c("ENTREZID", "SYMBOL"),
                  keytype = "SYMBOL")
  
  # Merge the EntrezGene IDs with the DESeq results
  res$SYMBOL <- res$X
  geneRes <- merge(res, res_entrezID)
  
  # Drop any NA values
  geneRes <- geneRes %>% drop_na(ENTREZID, log2FoldChange)
  
  # Remove dulicate Entrezgene IDs
  geneRes <- geneRes[!duplicated(geneRes[c("ENTREZID")]),]
  
  # Create a vector of the gene universe
  gsea_gene_list <- geneRes$log2FoldChange
  
  # Name vector with ENTREZ ids
  names(gsea_gene_list) <- geneRes$ENTREZID
  
  # sort the list in decreasing order (requiredProfiler)
  gsea_gene_list = sort(gsea_gene_list, decreasing = TRUE)
  
  # Retrieve gene sets from msigdb (human "Homo sapiens" or mouse "Mus musculus")
  msigdb_df <- msigdbr(species = "Homo sapiens", category = category, subcategory = subcategory) %>%
    dplyr::select(gs_name, entrez_gene)
  
  # GSEA analysis (https://rdrr.io/bioc/clusterProfiler/man/GSEA.html)
  gsea_res <- GSEA(gsea_gene_list, 
                   minGSSize = 15,
                   maxGSSize = 500,
                   TERM2GENE = msigdb_df,
                   eps = 1e-10,
                   pvalueCutoff = 0.05)
  
  # Transform the data to table format
  #gsea_res <- data.frame(gsea_res)
  
  return(gsea_res)
}

# Perform GSEA analysis based on need (define category and subcategory)
gsea_res_H <- msigdb_GSEA(res, "H", NULL)

# Readable core enrichment https://guangchuangyu.github.io/2016/07/leading-edge-analysis/
gsea_res_H <- setReadable(gsea_res_H, OrgDb = 'org.Hs.eg.db', keyType = "ENTREZID")
  
# Save the results
write.csv(gsea_res_H, "Training8_GSEA_Result.csv")
#################################### Bubble plot for GSEA ###########################################
library(ggplot2)

# Select the GSEA result files
res1 <- read.csv(file.choose(), header = TRUE)

# Hallmark gene category (https://www.cell.com/fulltext/S2405-4712(15)00218-5)
# Make bubble plots
ggplot(res1, aes(x = Condition, y = Description, size = NES, color = FDR)) +
  # Multi-level labels (https://dmitrijskass.netlify.app/2019/06/30/multi-level-labels-with-ggplot2/)
  facet_grid(cols = vars(Type), rows = vars(Category), 
             scales = "free", switch = "both")+
  geom_point(alpha = 0.4) +
  geom_text(aes(label = Positive_value), size = 8)+
  scale_size(range = c(0.1,10)) + #1, 20 for Hallmark similarity
  scale_color_gradient(low = "red", high = "blue")+
  theme(axis.line = element_line(),
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = "gray90"),
        # Change minor label size
        axis.text.x = element_text(size = 15),
        # Change major label size
        strip.text = element_text(size = 15),
        axis.title = element_blank())+
  ggtitle("Human Hallmark Gene Set")
