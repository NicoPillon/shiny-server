library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------------------------
library(readxl)

# human data
saveRDS(readRDS("../../R_databases/Muscle_Obesity_Transcriptomics/human/data_out/datamatrix_normalized_muscle.Rds"),
        file="data/human_datamatrix.Rds")
saveRDS(readRDS("../../R_databases/Muscle_Obesity_Transcriptomics/human/data_out/metadata_muscle.Rds"),
        file="data/human_metadata.Rds")

# mouse data