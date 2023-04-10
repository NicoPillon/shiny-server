library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)

# data
saveRDS(readRDS("../../R_databases/Exercise_circulating_EV_miRNA/data/Exercise_circulating_EV_miRNA_stats.Rds"),
        file="data/data.Rds")


#stats
saveRDS(read.csv("../../R_databases/Exercise_circulating_EV_miRNA/data/metaMiR.csv"),
        file="data/stats.Rds")
