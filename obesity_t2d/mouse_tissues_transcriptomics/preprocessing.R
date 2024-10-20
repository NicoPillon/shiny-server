library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------------------------
library(readxl)

# matrix
datamatrix <- readRDS("../../../R_databases/obesity_t2d/mouse_tissues_transcriptomics/data_out/datamatrix.Rds")
saveRDS(datamatrix, file="data/datamatrix.Rds")

# metadata
metadata <- readRDS("../../../R_databases/obesity_t2d/mouse_tissues_transcriptomics/data_out/metadata.Rds")
saveRDS(metadata, file="data/metadata.Rds")

# list of datasets
datasets <- readRDS("../../../R_databases/obesity_t2d/mouse_tissues_transcriptomics/data_out/references.Rds")
saveRDS(datasets, file="data/references.Rds")
