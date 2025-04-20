# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install and load required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Using version 3.16 of BiocManager which is compatible with R 4.2
BiocManager::install(version = "3.16")

if (!require("TCGAbiolinks", quietly = TRUE)) {
  BiocManager::install("TCGAbiolinks")
}

if (!require("maftools", quietly = TRUE)) {
  BiocManager::install("maftools")
}

# Load required libraries
library(TCGAbiolinks)
library(maftools)
library(dplyr)

# Query SNV mutation data for TCGA-BRCA
query_snv <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

print(query_snv)

# Download data
GDCdownload(
  query = query_snv, 
  method = "api", 
  files.per.chunk = 5,
  directory = "GDCdata"
)

# Prepare data
brca_maf <- GDCprepare(query_snv)
print("MAF data preparation complete")

# Read as MAF object for analysis with maftools
print(class(brca_maf))
print(head(colnames(brca_maf)))

if(!"Hugo_Symbol" %in% colnames(brca_maf)) {
  print("Converting to standard MAF format...")
}

# Read the MAF data 
maf_object <- read.maf(maf = brca_maf)
print("MAF object created successfully")

# Analyze co-occurrence and mutual exclusivity (performs pair-wise Fisher's Exact test) - top 25 genes
print("Analyzing co-occurrence and mutual exclusivity...")
cooc_me <- somaticInteractions(maf = maf_object, top = 25, pvalue = 0.05)

print(cooc_me)

# Export results to a CSV file
results_df <- as.data.frame(cooc_me)
write.csv(results_df, "BRCA_gene_interactions.csv", row.names = FALSE)
print("Results saved to BRCA_gene_interactions.csv")

# Visualize the results
pdf("BRCA_gene_interactions.pdf", width = 10, height = 10)
somaticInteractions(maf = maf_object, top = 25, pvalue = 0.05)
dev.off()
print("Results visualization saved to BRCA_gene_interactions.pdf")