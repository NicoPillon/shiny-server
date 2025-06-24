#--------------------------------------------------------------------------------------------------------
#
# Overview
#
#--------------------------------------------------------------------------------------------------------

# Load data and libraries
library(feather)
library(shinycssloaders)
library(tidyverse)

library(arrow)

library(shinyjs)
library(r3dmol)


# for plot
library(highcharter)

# to retrieve gene description
library(rentrez)
library(httr)
library(jsonlite)

gene_list <- readRDS("data/list_gene.Rds")
gene_protein_map <- readRDS("gene_uniprot_map.rds")

