# Load  libraries
library(TCGAbiolinks)
library(maftools)
library(dplyr)

# Path to previously downloaded data
gdc_data_dir <- "GDCdata"

# Read driver genes from IntOGen file
driver_genes <- read.delim("IntOGen-DriverGenes_BRCA.tsv", 
                          stringsAsFactors = FALSE,
                          quote = "\"")

# Check structure of the file
print("Column names in the file:")
print(colnames(driver_genes))

# Extract gene symbols, handling the quotes properly
if ("Symbol" %in% colnames(driver_genes)) {
  driver_gene_list <- driver_genes$Symbol
} else if ("\"Symbol\"" %in% colnames(driver_genes)) {
  # If column names still have quotes
  driver_gene_list <- driver_genes$`"Symbol"`
} else {
  # Try to find the column with gene symbols
  print("Available columns:")
  print(colnames(driver_genes))
  stop("Could not find Symbol column in driver genes file")
}

# Remove any quotes from gene symbols if present
driver_gene_list <- gsub("\"", "", driver_gene_list)

print(paste("Number of driver genes loaded:", length(driver_gene_list)))
print("First 10 driver genes:")
print(head(driver_gene_list, 10))

# Load MAF data
query_snv <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

maf_files <- GDCprepare(query_snv)
maf_object <- read.maf(maf = maf_files)
print("MAF object loaded successfully")

# Basic summary of the data
print(maf_object)

# Create subset of the MAF object with only driver genes
driver_maf <- subsetMaf(maf = maf_object, genes = driver_gene_list)
print("Driver genes MAF subset created")
print(paste("Number of samples in driver genes subset:", length(unique(driver_maf@clinical.data$Tumor_Sample_Barcode))))
print(paste("Number of genes in driver genes subset:", length(unique(driver_maf@data$Hugo_Symbol))))

# Analyze co-occurrence and mutual exclusivity among driver genes
driver_cooc_me <- somaticInteractions(maf = driver_maf, pvalue = 0.05)

# Export results to a CSV file
driver_results_df <- as.data.frame(driver_cooc_me)
write.csv(driver_results_df, "BRCA_driver_gene_interactions.csv", row.names = FALSE)
print("Driver gene interaction results saved to BRCA_driver_gene_interactions.csv")

# Create improved driver gene interaction plot as PNG
png("BRCA_driver_interaction.png", width = 1600, height = 1400, res = 150)

# Increase margins substantially to prevent label cutting
# Format: c(bottom, left, top, right)
par(mar=c(10, 12, 4, 2) + 0.1)  

# Call somaticInteractions with only the supported parameters
somaticInteractions(maf = driver_maf, pvalue = 0.05)

dev.off()
print("Driver interaction plot saved to BRCA_driver_interaction.png")

# Compare oncoplot for all genes vs driver genes
# Create the oncoplot and save as PNG
png("BRCA_oncoplot.png", width = 2400, height = 1600, res = 150)

# Set margins to ensure everything fits properly
par(mar = c(2, 6, 2, 2))

# Generate the oncoplot for driver genes
oncoplot(
  maf = driver_maf,
  top = 25,  # Show top 25 genes
  titleText = "Oncoplot of Driver Genes in Breast Cancer",
  drawRowBar = TRUE,
  drawColBar = TRUE
)

# Close the device
dev.off()
print("Oncoplot saved to BRCA_oncoplot.png")

# Additional analysis: Lollipop plots for key driver genes
# Select top 5 most frequently mutated driver genes
gene_counts <- table(driver_maf@data$Hugo_Symbol)
top_drivers <- names(sort(gene_counts, decreasing = TRUE)[1:5])
print("Creating individual lollipop plots for top 5 driver genes:")
print(top_drivers)

# Create lollipop plots with improved label handling
# Alternative approach without unsupported parameters
# Create lollipop plots with alternative approach for domain visualization
for(gene in top_drivers) {
  print(paste("Creating domain-focused lollipop plot for", gene))
  
  png_file <- paste0("BRCA_lollipop_domains_", gene, ".png")
  
  tryCatch({
    # Extra wide format
    png(png_file, width = 3200, height = 1000, res = 150)
    
    # Basic plot
    lollipopPlot(
      maf = maf_object, 
      gene = gene, 
      showDomainLabel = FALSE
    )
    
    dev.off()
    print(paste("Successfully created domain-focused plot for", gene))
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    print(paste("Error creating plot for", gene, ":", e$message))
  })
}