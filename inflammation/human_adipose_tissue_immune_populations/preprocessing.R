library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------------------------
library(readxl)

saveRDS(readRDS("../../R_databases/Adipose_tissue_FACS/FACS_WAT.Rds"),
        file="data/data_raw.Rds")

saveRDS(read_xlsx("../../R_databases/Adipose_tissue_FACS/FACS_WAT.xlsx"),
        file="data/samples.Rds")

