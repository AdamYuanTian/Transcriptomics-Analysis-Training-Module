###################################### Install packages ############################################
### Run only for the first time
install.packages("BiocManager")

BiocManager::install("DESeq2") # Note why this package installation is different from Line 2.

############################### Prepare data for your sample ####################################################
# This pipe line can be used for analyzing cancer samples
library(DESeq2) # For differential gene expression analysis


# Choose the files that contain the RNAseq counts and associated type (info table)
sample_counts <- file.choose()
sample_info <- file.choose()

# Read in the RNAseq counts and infotatble
sample_counts <- read.table(sample_counts, sep = ",", header = TRUE, row.names = 1)
sample_info <- read.csv(sample_info, header = TRUE)

# This is to make sure that the names in the info table match up with the column names of the counts table
sample_counts <-sample_counts[, as.character(sample_info$Description)]
all(colnames(sample_counts) == sample_info$Description)

# Convert the column from characters to factors (required by DESeq)
sample_info$Type <- factor(sample_info$Type)

############################## Getting TCGA data (normal and tumor) to compare with your sample ############
BiocManager::install("TCGAbiolinks") # For fetching data from TCGA database

library(TCGAbiolinks)

projectID <- "TCGA-COAD" # Change project ID based on cancer type, here I use colorectal cancer as an example
# Check out "https://portal.gdc.cancer.gov/" for more information about TCGA database

# Download RNA-seq data from TCGA database (run only once for the same data)
# Define a function for doing this
getRNAseqData <- function(projectID) {
  # The below function queries HT-Seq counts data from selected project 
  #in the Genomic Data Commons (GDC)
  query.RNASeq <- GDCquery(project = projectID,
                           data.category = "Transcriptome Profiling",
                           data.type = "Gene Expression Quantification",
                           workflow.type = "HTSeq - Counts")
  
  # This command downloads above query through GDC
  # Should only be done once
  GDCdownload(query.RNASeq)
  # After downloading the data to your local drive, you can create the RNAseqData file again by running the getRNAseqData function with the above line of code disabled.
  
  # Formats data together into a single object
  RNAseqData <- GDCprepare(query = query.RNASeq,
                           summarizedExperiment = TRUE)
  
  return(RNAseqData)
}
# Run the function
RNAseqData <- getRNAseqData(projectID)

# Organize RNA-seq data into a dataframe that later can be merged with PDX data
#Makes a data frame from Summarized Experiment (changes the type)
TCGACOUNTS <- data.frame(assay(RNAseqData))
colnames(TCGACOUNTS) <- colnames(assay(RNAseqData))
















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
