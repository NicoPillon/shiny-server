library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)

# data
saveRDS(readRDS("../../R_databases/VO2_range/data/pool_data.Rds"),
        file="data/data.Rds")


studies <- read_xlsx("../../R_databases/VO2_range/000_Studies.xlsx")[,c(1:6)]


createLink <- function(val) {
  sprintf(paste0('<a href="', URLdecode(val),'" target="_blank">', gsub("(.*org/)|(.*=)", "", val) ,'</a>'))
}
studies$DOI  <- sapply(studies$DOI, createLink)

saveRDS(studies, file = "data/studies.Rds")
