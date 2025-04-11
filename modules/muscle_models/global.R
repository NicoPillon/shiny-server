#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle cells
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
human_genelist <- readRDS("data/gene_names.Rds")
samples_list <- readRDS("data/samples_list.Rds")
human_references <- read_feather("data/dataset.feather")