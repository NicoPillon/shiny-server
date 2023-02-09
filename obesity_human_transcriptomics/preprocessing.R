library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------------------------
library(readxl)

saveRDS(readRDS("../../R_databases/Obesity_Human/Obesity_Human_data.Rds"),
        file="data/data_raw.Rds")


#list of samples and metadata
saveRDS(readRDS("../../R_databases/Obesity_Human/Obesity_Human_samples.Rds"),
        file="data/metadata.Rds")

# list of datasets
datasets <- read_xlsx("../../R_databases/Obesity_Human/Obesity_Human.xlsx")
datasets <- datasets[1:37, c(1,2,20)]
datasets <- datasets[!duplicated(datasets$GEO),]
saveRDS(datasets, file="data/datasets.Rds")
