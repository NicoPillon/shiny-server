library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)
library(feather)
library(tidyverse)

#--------------------------------------------------------------------------------------------------
# MATRIX
#--------------------------------------------------------------------------------------------------
datamatrix <- readRDS("../../rawdata/myotube_palmitate/data_out/datamatrix.Rds")
write_feather(datamatrix, "data/datamatrix.feather")

#--------------------------------------------------------------------------------------------------
# METADATA
#--------------------------------------------------------------------------------------------------
gene_list <- readRDS("../../rawdata/myotube_palmitate/data_out/gene_list.Rds")
saveRDS(gene_list, "data/list_gene.Rds")

#--------------------------------------------------------------------------------------------------
# METADATA
#--------------------------------------------------------------------------------------------------
metadata <- readRDS("../../rawdata/myotube_palmitate/data_out/metadata.Rds")
saveRDS(metadata, file="data/metadata.Rds")


#--------------------------------------------------------------------------------------------------
# REFERENCES
#--------------------------------------------------------------------------------------------------
references <- readRDS("../../rawdata/myotube_palmitate/data_out/references.Rds")
saveRDS(references, "data/references.Rds")