library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------------------------

#data
saveRDS(readRDS("../../R_databases/HSMC_Differentiation_T2D/data/BlastDiffT2D_data_norm.Rds"),
        file="data/data_raw.Rds")

#stats, overall effect
saveRDS(readRDS("../../R_databases/HSMC_Differentiation_T2D/stats/BlastDiffT2D_stats_differentiation.Rds"),
        file="data/data_stats_differentiation.Rds")

saveRDS(readRDS("../../R_databases/HSMC_Differentiation_T2D/stats/BlastDiffT2D_stats_interaction.Rds"),
        file="data/data_stats_interaction.Rds")

saveRDS(readRDS("../../R_databases/HSMC_Differentiation_T2D/stats/BlastDiffT2D_stats_T2D.Rds"),
        file="data/data_stats_T2D.Rds")

#stats, pair-wise comparison
saveRDS(readRDS("../../R_databases/HSMC_Differentiation_T2D/stats/BlastDiffT2D_stats_diffNGT.Rds"),
        file="data/data_stats_differentiation_NGT.Rds")

saveRDS(readRDS("../../R_databases/HSMC_Differentiation_T2D/stats/BlastDiffT2D_stats_diffT2D.Rds"),
        file="data/data_stats_differentiation_T2D.Rds")

saveRDS(readRDS("../../R_databases/HSMC_Differentiation_T2D/stats/BlastDiffT2D_stats_T2Dblast.Rds"),
        file="data/data_stats_T2D_blast.Rds")

saveRDS(readRDS("../../R_databases/HSMC_Differentiation_T2D/stats/BlastDiffT2D_stats_T2DTubes.Rds"),
        file="data/data_stats_T2D_tubes.Rds")

