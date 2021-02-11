# Make sure you have installed "BiocManager". If not run the following line for installation
install.packages("BiocManager")

# Install packages needed for this R script. 
# The same package installation is needed only once for the same computer.
BiocManager::install("DESeq2") # Note why this package installation is different from Line 2.

# Load the libraries we have installed
library(DESeq2) # An important library for differential gene analysis, but we will only need two simple functions from it this time


# Choose the files that contain the RNAseq counts and associated type (info table)
sample_counts <- file.choose()
sample_info <- file.choose()

# Read in the RNAseq counts and infotatble
sample_counts <- read.csv(sample_counts)
sample_info <- read.csv(sample_info)

# Don't worry if you don't understand these two lines, we can figure it out in the future
# This is to make sure that the names in the info table match up with the column names of the counts table
sample_counts <-sample_counts[, as.character(sample_info$Description)]
all(colnames(sample_counts) == sample_info$Description)

# Set up the dds object
sample_dds <- DESeqDataSetFromMatrix(countData = sample_counts,
                              colData = sample_info,
                              design = ~Type)


# Variance Stabilizing Transformation, for plotting purpose
sample_vst <- vst(sample_dds)

# PCA analysis and plot result
plotPCA(sample_vst, intgroup = "Type")
