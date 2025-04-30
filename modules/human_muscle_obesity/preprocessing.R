library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------------------------
library(readxl)
library(feather)
library(tidyverse)

# load data, split and save as feather files
datamatrix <- readRDS("../../../R_databases/obesity_t2d/human_muscle_transcriptomics/data_out/datamatrix_normalized_muscle.Rds") %>%
  data.frame()
n_genes <- nrow(datamatrix)
chunk_size <- ceiling(n_genes / 4)
splits <- split(1:n_genes, ceiling(seq_along(1:n_genes) / chunk_size))
file_names <- paste0("datamatrix_", seq_along(splits), ".feather")

# Save each chunk and create the mapping
gene_list <- do.call(rbind, lapply(seq_along(splits), function(i) {
  rows <- splits[[i]]
  chunk <- datamatrix[rows, ]
  write_feather(chunk, file.path("data", file_names[i]))
  
  data.frame(
    row = seq_len(nrow(chunk)),  # row index within this chunk
    SYMBOL = rownames(chunk),
    file = file_names[i],
    stringsAsFactors = FALSE
  )
}))
saveRDS(gene_list, "data/gene_list.Rds")

# Reference datasets
references <- readRDS("../../../R_databases/obesity_t2d/human_muscle_transcriptomics/data_out/references.Rds")
saveRDS(references, "data/references.Rds")

# metadata
metadata <- readRDS("../../../R_databases/obesity_t2d/human_muscle_transcriptomics/data_out/metadata_muscle.Rds")
saveRDS(metadata, "data/metadata.Rds")

# statistics - transcriptomics
stats.transcriptomics <- read.csv("../../rawdata/human_muscle_obesity/transcriptomics/data_out/stats_obesity_DEG.csv")
saveRDS(stats.transcriptomics, "data/stats.transcriptomics.Rds")
