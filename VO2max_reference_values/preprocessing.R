library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)
library(tidyverse)
createLink <- function(val) {
  sprintf(paste0('<a href="', URLdecode(val),'" target="_blank">', gsub("(.*org/)|(.*=)", "", val) ,'</a>'))
}

# data
VO2_data <- readRDS("../../R_databases/VO2max_reference_values/data_out/dat_adjusted.Rds")

VO2_data <- pivot_longer(VO2_data,
                         cols = c("P0.1", "P0.25", "P0.4", "P0.5", "P0.6", "P0.75", "P0.9"),
                         names_to = "Percentile",
                         values_to = "VO2max")
VO2_data$Percentile <- gsub("P0.5", "Median", VO2_data$Percentile)
VO2_data$Modality <-  str_to_title(VO2_data$Modality)


saveRDS(VO2_data, file = "data/data.Rds")


studies_included <- read.csv("../../R_databases/VO2max_reference_values/data_out/references_included.csv")
saveRDS(studies_included, file = "data/studies_included.Rds")

