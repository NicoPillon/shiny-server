library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------------------------
library(readxl)
library(feather)
library(tidyverse)

# matrix
datamatrix <- readRDS("../../../R_databases/exercise/mouse_exercise_tissues_transcriptomics/data_out/mouse_exercise_tissues_matrix.Rds")
write_feather(datamatrix[1:7000,], "data/datamatrix_1.feather")
write_feather(datamatrix[7001:14000,], "data/datamatrix_2.feather")
write_feather(datamatrix[14001:nrow(datamatrix),], "data/datamatrix_3.feather")

# metadata
metadata <- readRDS("../../../R_databases/exercise/mouse_exercise_tissues_transcriptomics/data_out/mouse_exercise_tissues_metadata.Rds")
saveRDS(metadata, file="data/metadata.Rds")

# list of datasets
datasets <- readRDS("../../../R_databases/exercise/mouse_exercise_tissues_transcriptomics/data_out/mouse_exercise_tissues_references.Rds")
saveRDS(datasets, file="data/references.Rds")

# list of genes
saveRDS(rownames(datamatrix), file="data/genelist.Rds")
