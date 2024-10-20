library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------------------------
library(readxl)
library(feather)
library(tidyverse)

datamatrix <- readRDS("../../../R_databases/obesity_t2d/human_transcriptomics/muscle/data_out/datamatrix_normalized_muscle.Rds") %>%
  data.frame()

metadata <- readRDS("../../../R_databases/obesity_t2d/human_transcriptomics/muscle/data_out/metadata_muscle.Rds")

references <- readRDS("../../../R_databases/obesity_t2d/human_transcriptomics/muscle/data_out/references.Rds")

write_feather(datamatrix[1:7000,], "data/human_datamatrix_1.feather")
write_feather(datamatrix[7001:nrow(datamatrix),], "data/human_datamatrix_2.feather")

saveRDS(rownames(datamatrix), "data/human_genelist.Rds")
saveRDS(metadata, "data/human_metadata.Rds")
saveRDS(references, "data/human_references.Rds")