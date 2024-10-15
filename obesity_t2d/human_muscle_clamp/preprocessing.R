library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------------------------

#meta-analysis
saveRDS(readRDS("../../R_databases/Muscle_Clamp/MuscleClamp_stats_MetaAnalysis.Rds"),
        file="data/stats_MetaAnalysis.Rds")

#insulin resistance
saveRDS(readRDS("../../R_databases/Muscle_Clamp/GSE22309/GSE22309_SYMBOL.Rds"),
        file="data/GSE22309_data.Rds")

saveRDS(readRDS("../../R_databases/Muscle_Clamp/GSE22309/GSE22309_stats_int.DIA.Rds"),
        file="data/GSE22309_stats_int.DIA.Rds")

saveRDS(readRDS("../../R_databases/Muscle_Clamp/GSE22309/GSE22309_stats_int.IR.Rds"),
        file="data/GSE22309_stats_int.IR.Rds")

saveRDS(readRDS("../../R_databases/Muscle_Clamp/GSE22309/GSE22309_stats_IR.Rds"),
        file="data/GSE22309_stats_IR.Rds")

saveRDS(readRDS("../../R_databases/Muscle_Clamp/GSE22309/GSE22309_stats_IS.Rds"),
        file="data/GSE22309_stats_IS.Rds")

saveRDS(readRDS("../../R_databases/Muscle_Clamp/GSE22309/GSE22309_stats_DIA.Rds"),
        file="data/GSE22309_stats_DIA.Rds")

#time course
saveRDS(readRDS("../../R_databases/Muscle_Clamp/MuscleClamp_stats_TimeCourse.Rds"),
        file="data/stats_TimeCourse.Rds")