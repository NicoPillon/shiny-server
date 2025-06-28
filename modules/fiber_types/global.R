#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle cells - Global
#
#--------------------------------------------------------------------------------------------------------

#======================
# Load Required Packages
#======================
# Efficient disk-based dataset access (used to load parquet files)
library(arrow)

# For advanced shiny options
library(shinycssloaders)
library(shinyWidgets)

# Data wrangling and transformation (includes dplyr, tidyr, readr, etc.)
library(tidyverse)

# For advanced ggplot geoms like geom_sina (used for scatter overlay)
library(ggforce)

# Interactive tables
library(DT)

#======================
# Load Static Data Files
#======================

# metadata
metadata_proteome <- readRDS("data/metadata_proteome.Rds")
metadata_transcriptome <- readRDS("data/metadata_transcriptome.Rds")

metadata_proteome <- data.frame(
  study = metadata_proteome$study,
  sample_id = metadata_proteome$sample_id,
  fiber_type = metadata_proteome$FiberType,
  OMICS = "Proteome"
)

metadata_transcriptome <- data.frame(
  study = metadata_transcriptome$GEO,
  sample_id = metadata_transcriptome$geo_accession,
  fiber_type = metadata_transcriptome$FiberType,
  OMICS = "Transcriptome"
)

metadata <- rbind(metadata_proteome, metadata_transcriptome)


# List of genes
gene_to_file_proteome <- readRDS("data/list_gene_proteome.Rds")
gene_to_file_proteome$OMICS <- "Proteome"

gene_to_file_transcriptome <- readRDS("data/list_gene_transcriptome.Rds")
gene_to_file_transcriptome$OMICS <- "Transcriptome"

gene_to_file <- rbind(gene_to_file_proteome, gene_to_file_transcriptome)


# references
references <- readRDS("data/references.Rds")
