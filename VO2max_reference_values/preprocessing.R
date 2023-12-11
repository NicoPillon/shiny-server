library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)
createLink <- function(val) {
  sprintf(paste0('<a href="', URLdecode(val),'" target="_blank">', gsub("(.*org/)|(.*=)", "", val) ,'</a>'))
}

# data
saveRDS(readRDS("../../R_databases/VO2max_reference_values/data/pool_data.Rds"),
        file="data/data.Rds")


studies_included <- read_xlsx("../../R_databases/VO2max_reference_values/000_Studies.xlsx", sheet = 1)[,c(1:6)]
studies_included$DOI  <- sapply(studies_included$DOI, createLink)
saveRDS(studies_included, file = "data/studies_included.Rds")


studies_excluded <- read_xlsx("../../R_databases/VO2max_reference_values/000_Studies.xlsx", sheet = 2)[,c(1,6,7)]
studies_excluded$DOI  <- sapply(studies_excluded$DOI, createLink)
saveRDS(studies_excluded, file = "data/studies_excluded.Rds")



