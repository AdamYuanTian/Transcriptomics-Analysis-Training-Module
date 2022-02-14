# Make Euler Diagram to show common differentially expressed genes (DEGs) between conditions
#install.packages("eulerr")
library(eulerr)

# Select the DEGs result files for comparison
FileA <- file.choose() #Select the 1st DEGs result, Exp_A vs Ctrl_A
FileB <- file.choose() #Select the 1st DEGs result, Exp_B vs Ctrl_B

DEGs_A <- read.csv(FileA, header = TRUE)
DEGs_B <- read.csv(FileB, header = TRUE)

# Define a function to get significant DEGs
sig_DEGs <- function(DEGs){
  df <- data.frame(DEGs)
  sig_gene_df = subset(df, padj < 0.05)
  sig_gene_df = subset(sig_gene_df, baseMean > 10)
  sig_gene_df = subset(sig_gene_df, abs(log2FoldChange) > 1.5)
  return(sig_gene_df$X)
}

# Get significant DEGs
A <- sig_DEGs(DEGs_A)
B <- sig_DEGs(DEGs_B)

#Get the intersection and difference of A and B
AnotB <-setdiff(A, B)
BnotA <-setdiff(B,A)
AandB <-intersect(A, B)

eulerplot <-euler(c("Exp_Ctrl_A"=length(AnotB), 
                    "Exp_Ctrl_B"=length(BnotA), 
                    "Exp_Ctrl_A&Exp_Ctrl_B"=length(AandB))
                  )

# https://rdrr.io/cran/eulerr/man/plot.euler.html
# https://cran.r-project.org/web/packages/eulerr/vignettes/gallery.html
plot(eulerplot, 
     main = "Training5_EulerDiagram",
     quantities = list(type = "counts"),
     edges = c("blue", "red"),
     fills = c("white", "white"),
     lwd = 4,
     legend = list(side = "right"))

     