library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)

# data
saveRDS(readRDS("../../../R_databases/other/Muscle_AMPK/data_out/muscleAMPK_datamatrix.Rds"),
        file="data/AMPK_datamatrix.Rds")

# metadata
saveRDS(readRDS("../../../R_databases/other/Muscle_AMPK/data_out/muscleAMPK_metadata.Rds"),
        file="data/AMPK_metadata.Rds")