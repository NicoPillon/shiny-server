#==================================================================================================
# Preprocessing Script for Shiny App
#==================================================================================================

# Set working directory to the script's location (RStudio only)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

#==================================================================================================
# Load Required Libraries
#==================================================================================================

library(readxl)     # For reading Excel files (not used here, but possibly needed upstream)
library(arrow)      # For writing Parquet files (used for efficient data access in the app)
library(tidyverse)  # For data manipulation (dplyr, purrr, etc.)


#==================================================================================================
# SECTION 1: PROCESS AND EXPORT EXPRESSION MATRIX
#==================================================================================================

# Load full expression matrix (genes x samples)
datamatrix <- readRDS("../../rawdata/fiber_types/data_out/transcriptome_matrix.Rds") %>%
  data.frame()

# Determine how many genes to include per chunk to divide the data in ~3 parts
n_genes <- nrow(datamatrix)
chunk_size <- ceiling(n_genes / 8)  # Number of genes per file
splits <- split(1:n_genes, ceiling(seq_along(1:n_genes) / chunk_size))  # Split row indices
file_names <- paste0("datamatrix_transcriptome_", seq_along(splits), ".parquet")      # Output file names

# Process each chunk: attach gene names, reorder columns, save as Parquet
gene_list <- do.call(rbind, lapply(seq_along(splits), function(i) {
  rows <- splits[[i]]
  chunk <- datamatrix[rows, ]
  
  # Add gene symbol as a dedicated column to enable filtering in `arrow::open_dataset()`
  chunk$TARGET <- rownames(chunk)
  chunk <- chunk[, c("TARGET", setdiff(colnames(chunk), "TARGET"))]  # Ensure TARGET is the first column
  
  # Write chunk as a Parquet file (efficient, columnar format)
  write_parquet(chunk, file.path("data", file_names[i]))
  
  # Return mapping: which gene is in which file
  data.frame(
    TARGET = rownames(chunk),
    file = file_names[i],
    stringsAsFactors = FALSE
  )
}))

# Save lookup table to match genes to their file location
saveRDS(gene_list, "data/list_gene_transcriptome.Rds")

#==================================================================================================
# SECTION 2: EXPORT METADATA
#==================================================================================================

# Load sample-level metadata (sample ID, treatment, cell type, etc.)
metadata <- readRDS("../../rawdata/fiber_types/data_out/transcriptome_metadata.Rds")

# Save to app directory
saveRDS(metadata, file = "data/metadata_transcriptome.Rds")

#==================================================================================================
# SECTION 3: EXPORT REFERENCE TABLE
#==================================================================================================
# Load reference table containing study-level information (e.g. publication info)
references <- read_xlsx("../../rawdata/fiber_types/references.xlsx")
saveRDS(references, "data/references.Rds")