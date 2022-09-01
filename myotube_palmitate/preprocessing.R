library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)

# data
saveRDS(readRDS("../../R_databases/Myotube_Palmitate/Myotube_Palmitate_norm.Rds"),
        file="data/data.Rds")


#stats
saveRDS(readRDS("../../R_databases/Myotube_Palmitate/Myotube_Palmitate_stats.Rds"),
        file="data/stats.Rds")
