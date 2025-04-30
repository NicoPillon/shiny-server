library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################

#----------------------------------------------------------------------------------------------
# Mouse data
#----------------------------------------------------------------------------------------------
saveRDS(readRDS("../../R_databases/Macrophage_Polarization/mouse/data_batch_corrected.Rds"),
        file="data/mouse_data.Rds")

# Mouse stats
saveRDS(readRDS("../../R_databases/Macrophage_Polarization/mouse/Stats_M1vsM0.Rds"),
        file="data/mouse_M1vsM0.Rds")
saveRDS(readRDS("../../R_databases/Macrophage_Polarization/mouse/Stats_M2vsM0.Rds"),
        file="data/mouse_M2vsM0.Rds")
saveRDS(readRDS("../../R_databases/Macrophage_Polarization/mouse/Stats_M1vsM2.Rds"),
        file="data/mouse_M1vsM2.Rds")


#----------------------------------------------------------------------------------------------
# Human data
#----------------------------------------------------------------------------------------------
saveRDS(readRDS("../../R_databases/Macrophage_Polarization/human/data_batch_corrected.Rds"),
        file="data/human_data.Rds")

# Human stats
saveRDS(readRDS("../../R_databases/Macrophage_Polarization/human/Stats_M1vsM0.Rds"),
        file="data/human_M1vsM0.Rds")
saveRDS(readRDS("../../R_databases/Macrophage_Polarization/human/Stats_M2vsM0.Rds"),
        file="data/human_M2vsM0.Rds")
saveRDS(readRDS("../../R_databases/Macrophage_Polarization/human/Stats_M1vsM2.Rds"),
        file="data/human_M1vsM2.Rds")


#----------------------------------------------------------------------------------------------
# datasets
#----------------------------------------------------------------------------------------------
saveRDS(read.csv("../../R_databases/Macrophage_Polarization/Macrophage_polarization_datasets_2021-08-31.csv"),
        file="data/datasets.Rds")

