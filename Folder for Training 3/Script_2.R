# Make sure you have installed "BiocManager". If not run the following line for installation
install.packages("BiocManager")

# Install packages needed for this R script. 
# The same package installation is needed only once for the same computer.
BiocManager::install("DESeq2") # Note why this package installation is different from Line 2.

# Load the libraries we have installed
library(DESeq2) # For differential gene expression analysis


# Choose the files that contain the RNAseq counts and associated type (info table)
sample_counts <- file.choose()
sample_info <- file.choose()

# Read in the RNAseq counts and infotatble
sample_counts <- read.table(sample_counts, sep = ",", header = TRUE, row.names = 1)
sample_info <- read.csv(sample_info)

# This is to make sure that the names in the info table match up with the column names of the counts table
sample_counts <-sample_counts[, as.character(sample_info$Description)]
all(colnames(sample_counts) == sample_info$Description)

# Convert the column from characters to factors (required by DESeq)
sample_info$Type <- factor(sample_info$Type)

# Set up the dds object
sample_dds <- DESeqDataSetFromMatrix(countData = sample_counts,
                                     colData = sample_info,
                                     design = ~Type)


# Run DESeq (this line of code is short, but it automatically performs all the differential expression analysis)
sample_dds <- DESeq(sample_dds)

# Check result. Differential expression is comparing difference between two conditions.
# Comparing difference between Type A and Type B (as defined in InfoTable)
sample_res <- results(sample_dds, contrast = c("Type", "A", "B")) 

# Save the results to a table
res <- data.frame(sample_res)
write.csv(res, "Training3_DESeq2_Result.csv")
