library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)

# C2C12
saveRDS(readRDS("../../R_databases/Hypoxia/Mouse_C2C12_hypoxia/C2C12_SYMBOL_norm.Rds"),
        file="data/C2C12_data.Rds")
saveRDS(readRDS("../../R_databases/Hypoxia/Mouse_C2C12_hypoxia/C2C12_stats_hypoxia.Rds"),
        file="data/C2C12_stats.Rds")

# Human endothelial cells
  saveRDS(readRDS("../../R_databases/Hypoxia/Human_endothelial_cells/EC_hypoxia_SYMBOL_norm.Rds"),
        file="data/EC_data.Rds")
saveRDS(readRDS("../../R_databases/Hypoxia/Human_endothelial_cells/EC_hypoxia_stats.Rds"),
        file="data/EC_stats.Rds")

# Human monocye-derived macrophages
saveRDS(readRDS("../../R_databases/Hypoxia/Human_monocyte-derived_macrophages/HMDM_hypoxia_SYMBOL_norm.Rds"),
        file="data/HMDM_data.Rds")
saveRDS(readRDS("../../R_databases/Hypoxia/Human_monocyte-derived_macrophages/HMDM_hypoxia_stats.Rds"),
        file="data/HMDM_stats.Rds")

# Mouse_acute_hypoxia/GSE81286
saveRDS(readRDS("../../R_databases/Hypoxia/Mouse_acute_hypoxia/GSE81286_data.Rds"),
        file="data/GSE81286_data.Rds")
saveRDS(readRDS("../../R_databases/Hypoxia/Mouse_acute_hypoxia/GSE81286_stats2h.Rds"),
        file="data/GSE81286_stats2h.Rds")
saveRDS(readRDS("../../R_databases/Hypoxia/Mouse_acute_hypoxia/GSE81286_stats6h.Rds"),
        file="data/GSE81286_stats6h.Rds")




