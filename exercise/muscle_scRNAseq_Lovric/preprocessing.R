library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)

# data
saveRDS(readRDS("../../R_databases/SingleCellRNAseq/Lovric_exercise/Lovric_celltypes.Rds"),
        file="data/data.Rds")

