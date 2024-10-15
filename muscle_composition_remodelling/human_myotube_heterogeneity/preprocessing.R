library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------------------------

saveRDS(readRDS("../../R_databases/HSMC_heterogeneity/HSMC_data_corrected.Rds"),
        file="data/data_corrected.Rds")

#sex comparison
saveRDS(readRDS("../../R_databases/HSMC_heterogeneity/HSMC_sex_data.Rds"),
        file="data/sex_data.Rds")

saveRDS(readRDS("../../R_databases/HSMC_heterogeneity/HSMC_sex_stats.Rds"),
        file="data/sex_stats.Rds")





