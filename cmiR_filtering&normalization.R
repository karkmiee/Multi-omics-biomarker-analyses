#Title: "microRNA data normalization"



#Data: Next generation sequencing serum microRNA content.
# This script is used to make a filtered and normalized c-miR count file for downstream analyses.

#Abbreviations used in script:
#microRNA = miR

#------------------------------------------------------------------------------------------

#Step 1. Install packages and import raw cmiR data and phenodata

# Load packages
library(dplyr)
library(DESeq2) # Tool for normalization
library(sva)  #tool to take batch effect into account
library(edgeR) # Tool for DE analysis

# Import raw counts file
counts <- read.csv("cmiR1_filtered.txt",header=TRUE,sep="\t")

# Import study design file
targets<- read.table("phenoData_filtered.txt", header = TRUE, sep = "\t")


##---------------------------------------------------------------------------

#Step 2. Use DEseq2 to normalize miR read counts and SVA to account the batch effect. 


#Filter out miRs where counts is less than 1 in 50% of the samples
keep <- rowSums(counts < 1) >= 57
Counts <- counts[!keep, ]  #317 miRs left from total n. 1473


Counts_matrix <- as.matrix(Counts) # Convert dataframe 'Counts' to a matrix for ComBat_seq function

#Use ComBat_seq to correct batch effects
batch <- as.character(targets$Batch) # batch variable
adjusted_counts <- ComBat_seq(Counts_matrix, batch, group=NULL)

# Use DESeq2 to create normalized read counts
# Setup design matrix for DE-analysis. Note: DE analysis is not necessary for the normalization.
condition <- as.character(targets$Sex) # condition of interest (Sex)
design <- data.frame(condition = as.factor(condition)) # DE-design for the analysis
rownames(design) <- colnames(adjusted_counts)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = adjusted_counts, colData = design, design = ~ condition)

# DESeq2 DE-analysis
dds <- DESeq(dds)
dds

# Variance stabilizing transformation (VST)
vst_data <- varianceStabilizingTransformation(dds)

# Save normalized counts
normalized_counts <- assay(vst_data)

#Save normalized miR counts in a txt-file
write.table(normalized_counts, file="normalized_miR_counts.txt", sep="\t")

