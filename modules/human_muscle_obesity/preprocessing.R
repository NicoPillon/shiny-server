library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------------------------
library(readxl)
library(feather)
library(tidyverse)


#--------------------------------------------------------------------------------------------------
# TRANSCRIPTOMICS
#--------------------------------------------------------------------------------------------------
# load data, split and save as feather files
datamatrix <- readRDS("../transcriptomics/data_out/datamatrix_normalized_muscle.Rds") %>%
  data.frame()
n_genes <- nrow(datamatrix)
chunk_size <- ceiling(n_genes / 4)
splits <- split(1:n_genes, ceiling(seq_along(1:n_genes) / chunk_size))
file_names <- paste0("transcriptomics_datamatrix_", seq_along(splits), ".feather")

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

# statistics - transcriptomics
stats.transcriptomics <- read.csv("../transcriptomics/data_out/stats_obesity_DEG.csv")
saveRDS(stats.transcriptomics, "data/transcriptomics_stats.Rds")



#--------------------------------------------------------------------------------------------------
# METABOLOMICS
#--------------------------------------------------------------------------------------------------
# load data, split and save as feather files
datamatrix <- readRDS("../metabolomics/data_out/datamatrix_normalized_muscle.Rds") %>%
  data.frame()
n_genes <- nrow(datamatrix)
chunk_size <- ceiling(n_genes / 1)
splits <- split(1:n_genes, ceiling(seq_along(1:n_genes) / chunk_size))
file_names <- paste0("metabolomics_datamatrix_", seq_along(splits), ".feather")

# Save each chunk and create the mapping
metabolite_list <- do.call(rbind, lapply(seq_along(splits), function(i) {
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
saveRDS(metabolite_list, "data/list_metabolite.Rds")

# statistics - metabolomics
stats.metabolomics <- read.csv("../metabolomics/data_out/stats_obesity_DEM.csv")
saveRDS(stats.metabolomics, "data/metabolomics_stats.Rds")



#--------------------------------------------------------------------------------------------------
# Metadata
#--------------------------------------------------------------------------------------------------
metadata_metabolomics <- readRDS("../metabolomics/data_out/metadata_muscle.Rds")
metadata_transcriptomics <- readRDS("../transcriptomics/data_out/metadata_muscle.Rds")
metadata <- full_join(metadata_transcriptomics,
                      metadata_metabolomics)
saveRDS(metadata, "data/metadata.Rds")


#--------------------------------------------------------------------------------------------------
# REFERENCES
#--------------------------------------------------------------------------------------------------
references_transcriptomics <- readRDS("../transcriptomics/data_out/references.Rds")
references_metabolomics <- readRDS("../metabolomics/data_out/references.Rds")
references <- full_join(references_transcriptomics,
                        references_metabolomics)

saveRDS(references, "data/references.Rds")

