library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)
# Mouse data
saveRDS(readRDS("../../R_databases/Muscle_Injury/data_processed/MuscleInjury_data.Rds"),
        file="data/MuscleInjury_data.Rds")
saveRDS(readRDS("../../R_databases/Muscle_Injury/data_processed/MuscleInjury_stats.Rds"),
        file="data/MuscleInjury_stats.Rds")
saveRDS(readRDS("../../R_databases/Muscle_Injury/data_processed/MuscleInjury_samples.Rds"),
        file="data/MuscleInjury_samples.Rds")


#---------------------------------------------------------------------
#dataset
#---------------------------------------------------------------------
datasets <- read_xlsx("../../R_databases/Muscle_Injury/Datasets.xlsx")

#obese mice are excluded
datasets <- datasets[!datasets$Weight %in% "OBE",]

#remove weight column
datasets$Weight <- NULL

datasets$Treatment
datasets$Treatment <- c("Eccentric contraction", "Freeze injury", "Blunt injury", "Treadmill exercise", "Incision injury")
saveRDS(datasets, file="data/MuscleInjury_datasets.Rds")
