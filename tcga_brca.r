# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install and load required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Using version 3.16 of BiocManager, compatible with R 4.2
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
