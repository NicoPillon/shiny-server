library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)

# data
saveRDS(readRDS("../../R_databases/Myotube_Palmitate/C2C12_fattyAcids_norm.Rds"),
        file="data/data.Rds")


#stats
saveRDS(readRDS("../../R_databases/Myotube_Palmitate/C2C12_palmitate_stats.Rds"),
        file="data/stats.Rds")
