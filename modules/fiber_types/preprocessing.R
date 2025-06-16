#-----------------------------------------------------------------------
#
# Profile of skeletal muscle fibers
#
#----------------------------------------------------------------------

# set working directory to source file location
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

# load libraries
library(readxl)
library(feather)
library(tidyverse)

#---------------------------------------------------------------------------
# PROTEOME
#---------------------------------------------------------------------------
# metadata
metadata_proteome <- readRDS("../../rawdata/fiber_types/data_out/proteome_metadata.Rds")
saveRDS(metadata_proteome[,c(1,2,13)], "data/metadata_proteome.Rds")

# load data, split and save as feather files
data_proteome <- readRDS("../../rawdata/fiber_types/data_out/proteome_matrix.Rds") %>%
  data.frame()
n_genes <- nrow(data_proteome)
chunk_size <- ceiling(n_genes / 1)
splits <- split(1:n_genes, ceiling(seq_along(1:n_genes) / chunk_size))
file_names <- paste0("data_proteome_", seq_along(splits), ".feather")

# Save each chunk and create the mapping
gene_list <- do.call(rbind, lapply(seq_along(splits), function(i) {
  rows <- splits[[i]]
  chunk <- data_proteome[rows, ]
  write_feather(chunk, file.path("data", file_names[i]))
  
  data.frame(
    row = seq_len(nrow(chunk)),  # row index within this chunk
    SYMBOL = rownames(chunk),
    file = file_names[i],
    stringsAsFactors = FALSE
  )
}))
saveRDS(gene_list, "data/gene_list_proteome.Rds")


#---------------------------------------------------------------------------
# TRANSCRIPTOME
#---------------------------------------------------------------------------
# metadata
metadata_transcriptome <- readRDS("../../rawdata/fiber_types/data_out/transcriptome_metadata.Rds")
saveRDS(metadata_transcriptome[,c(1:8,19)], "data/metadata_transcriptome.Rds")

# load data, split and save as feather files
data_transcriptome <- readRDS("../../rawdata/fiber_types/data_out/transcriptome_matrix.Rds") %>%
  data.frame()
n_genes <- nrow(data_transcriptome)
chunk_size <- ceiling(n_genes / 5)
splits <- split(1:n_genes, ceiling(seq_along(1:n_genes) / chunk_size))
file_names <- paste0("data_transcriptome_", seq_along(splits), ".feather")

# Save each chunk and create the mapping
gene_list <- do.call(rbind, lapply(seq_along(splits), function(i) {
  rows <- splits[[i]]
  chunk <- data_transcriptome[rows, ]
  write_feather(chunk, file.path("data", file_names[i]))
  
  data.frame(
    row = seq_len(nrow(chunk)),  # row index within this chunk
    SYMBOL = rownames(chunk),
    file = file_names[i],
    stringsAsFactors = FALSE
  )
}))
saveRDS(gene_list, "data/gene_list_transcriptome.Rds")



#---------------------------------------------------------------------------
# References
#---------------------------------------------------------------------------
references <- read_xlsx("../../rawdata/fiber_types/references.xlsx")
saveRDS(references, "data/references.Rds")
