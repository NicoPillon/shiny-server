library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------------------------
library(readxl)

saveRDS(readRDS("../../../R_databases/inflammation/human_adipose_tissue_immune_populations/data_out/matrix.Rds"),
        file="data/matrix.Rds")

saveRDS(readRDS("../../../R_databases/inflammation/human_adipose_tissue_immune_populations/data_out/metadata.Rds"),
        file="data/metadata.Rds")

