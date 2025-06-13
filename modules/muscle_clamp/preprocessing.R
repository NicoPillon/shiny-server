library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)
library(feather)
library(tidyverse)

#--------------------------------------------------------------------------------------------------
# MATRIX
#--------------------------------------------------------------------------------------------------
# load data, split and save as feather files
datamatrix <- readRDS("../../rawdata/muscle_clamp/data_out/datamatrix.Rds") %>%
  data.frame()
n_genes <- nrow(datamatrix)
chunk_size <- ceiling(n_genes / 3)
splits <- split(1:n_genes, ceiling(seq_along(1:n_genes) / chunk_size))
file_names <- paste0("datamatrix_", seq_along(splits), ".feather")

# Save each chunk and create the mapping
gene_list <- do.call(rbind, lapply(seq_along(splits), function(i) {
  rows <- splits[[i]]
  chunk <- datamatrix[rows, ]
  write_feather(chunk, file.path("data", file_names[i]))
  
  data.frame(
    row = seq_len(nrow(chunk)),  # row index within this chunk
    TARGET = rownames(chunk),
    file = file_names[i],
    stringsAsFactors = FALSE
  )
}))
saveRDS(gene_list, "data/list_gene.Rds")

#--------------------------------------------------------------------------------------------------
# METADATA
#--------------------------------------------------------------------------------------------------
metadata <- readRDS("../../rawdata/muscle_clamp/data_out/metadata.Rds")
saveRDS(metadata, file="data/metadata.Rds")


#--------------------------------------------------------------------------------------------------
# REFERENCES
#--------------------------------------------------------------------------------------------------
references <- readRDS("../../rawdata/muscle_clamp/data_out/references.Rds")
saveRDS(references, "data/references.Rds")