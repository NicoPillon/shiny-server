library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------------------------
library(readxl)

saveRDS(readRDS("../../R_databases/Mouse_obesity/RNAseq/RNAseq_SYMBOL.Rds"),
        file="data/data_raw.Rds")



datasets <- read_xlsx("../../R_databases/Mouse_obesity/Mouse obesity.xlsx")
datasets <- datasets[1:9,]            
saveRDS(datasets, file="data/datasets.Rds")
