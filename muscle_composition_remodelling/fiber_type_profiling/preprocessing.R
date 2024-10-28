library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)

# data
saveRDS(readRDS("../../../R_databases/muscle_composition_remodelling/fiber_type_profile/proteome_matrix.Rds"),
        file="data/data_proteome.Rds")
saveRDS(readRDS("../../../R_databases/muscle_composition_remodelling/fiber_type_profile/transcriptome_matrix.Rds"),
        file="data/data_transcriptome.Rds")


# metadata
saveRDS(readRDS("../../../R_databases/muscle_composition_remodelling/fiber_type_profile/proteome_metadata.Rds"),
        file="data/metadata_proteome.Rds")
saveRDS(readRDS("../../../R_databases/muscle_composition_remodelling/fiber_type_profile/transcriptome_metadata.Rds"),
        file="data/metadata_transcriptome")
