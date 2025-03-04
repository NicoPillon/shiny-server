library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------------------------
library(readxl)
library(feather)
library(tidyverse)

datamatrix <- readRDS("../../../R_databases/muscle_composition_remodelling/human_muscle_aging/data_out/AgingDataset_datamatrix_normalized_muscle.Rds") %>%
  data.frame()

metadata <- readRDS("../../../R_databases/muscle_composition_remodelling/human_muscle_aging/data_out/AgingDataset_metadata_muscle.Rds")

references <- readRDS("../../../R_databases/muscle_composition_remodelling/human_muscle_aging/data_out/references.Rds")

write_feather(datamatrix[1:5000,], "data/human_datamatrix_1.feather")
write_feather(datamatrix[5001:10000,], "data/human_datamatrix_2.feather")
write_feather(datamatrix[10001:15000,], "data/human_datamatrix_3.feather")
write_feather(datamatrix[15001:nrow(datamatrix),], "data/human_datamatrix_4.feather")

saveRDS(rownames(datamatrix), "data/human_genelist.Rds")
saveRDS(metadata, "data/human_metadata.Rds")
saveRDS(references, "data/human_references.Rds")

