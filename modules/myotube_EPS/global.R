#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic myotube response to electrical pulse stimulation
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

# Sample-level metadata (used for filtering, plotting, and statistics)
metadata <- readRDS("data/metadata.Rds")

# References table with dataset descriptions
references <- readRDS("data/references.Rds")

# Lookup table that maps gene symbols to files and rows
gene_list <- readRDS("data/list_gene.Rds")
