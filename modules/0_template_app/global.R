#--------------------------------------------------------------------------------------------------------
#
# App short description
#
#--------------------------------------------------------------------------------------------------------

# Load data and libraries
library(feather)
library(shinycssloaders)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(cowplot)
library(DT)
library(matrixStats)
library(pheatmap)
library(arrow)

# Load static reference data
gene_list <- readRDS("data/gene_list.Rds")
sample_list <- readRDS("data/sample_list.Rds")
references <- read_feather("data/dataset.feather")