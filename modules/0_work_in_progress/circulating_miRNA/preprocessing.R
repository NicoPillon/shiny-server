library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)

#datasets
saveRDS(readRDS("../../../R_databases/exercise/human_exercise_circulating_miRNA/data/datasets.Rds"),
        file="data/datasets.Rds")

# data
saveRDS(readRDS("../../../R_databases/exercise/human_exercise_circulating_miRNA/data/Exercise_circulating_EV_miRNA_stats.Rds"),
        file="data/data.Rds")

#stats
saveRDS(read.csv("../../../R_databases/exercise/human_exercise_circulating_miRNA/data/metaMiR.csv"),
        file="data/stats.Rds")
