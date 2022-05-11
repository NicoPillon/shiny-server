library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)

# data
saveRDS(readRDS("../../R_databases/Muscle_AMPK/Muscle_AMPK.Rds"),
        file="data/AMPK_data.Rds")


#stats
saveRDS(readRDS("../../R_databases/Muscle_AMPK/GSE107212_stats_A1.Rds"),
        file="data/AMPK_stats_a1.Rds")
saveRDS(readRDS("../../R_databases/Muscle_AMPK/GSE107212_stats_A2.Rds"),
        file="data/AMPK_stats_a2.Rds")
saveRDS(readRDS("../../R_databases/Muscle_AMPK/GSE61904_stats.Rds"),
        file="data/AMPK_stats_a1a2.Rds")
saveRDS(readRDS("../../R_databases/Muscle_AMPK/GSE4067_stats_KO.Rds"),
        file="data/AMPK_stats_g3KO.Rds")
saveRDS(readRDS("../../R_databases/Muscle_AMPK/GSE4067_stats_TG.Rds"),
        file="data/AMPK_stats_g3TG.Rds")
saveRDS(readRDS("../../R_databases/Muscle_AMPK/AICAR_stats.Rds"),
        file="data/AMPK_stats_AICAR.Rds")