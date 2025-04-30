library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)

# data
saveRDS(readRDS("../../R_databases/Myotube_EPS/Myotube_EPS_norm.Rds"),
        file="data/data.Rds")

