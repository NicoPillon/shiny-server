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
VO2_data <- readRDS("../../../R_databases/exercise/VO2max_reference_values/data_out/data_raw.Rds")

VO2_data <- pivot_longer(VO2_data,
                         cols = c("P0.01", "P0.1", "P0.25", "P0.4", "P0.5", "P0.6", "P0.75", "P0.9", "P0.99"),
                         names_to = "Percentile",
                         values_to = "VO2max")
VO2_data$Percentile <- gsub("P0.5", "Median", VO2_data$Percentile)
VO2_data$Modality <-  str_to_title(VO2_data$Modality)

# make factor to place mean in between percentiles
VO2_data$Percentile <- factor(VO2_data$Percentile,
                              levels = c("P0.99", "P0.9", "P0.75", "P0.6",
                                         "Median", 
                                         "P0.4", "P0.25", "P0.1", "P0.01"))

saveRDS(VO2_data, file = "data/data.Rds")

# table of studies included
studies_included <- read_xlsx("../../../R_databases/exercise/VO2max_reference_values/data_out/references_included.xlsx")

# make linkg to pubmed
studies_included$link <- paste("https://pubmed.ncbi.nlm.nih.gov/", studies_included$PMID, sep="")
studies_included$link  <- sapply(studies_included$link, createLink)

# remove columns
colnames(studies_included)
studies_included[,c("Recruitment", "Male, Mean ± Sd", "Female, Mean ± Sd", "Type")] <- NULL

saveRDS(studies_included, file = "data/studies_included.Rds")

