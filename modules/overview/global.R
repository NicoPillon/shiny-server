#--------------------------------------------------------------------------------------------------------
#
# Overview
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
library(httr)
library(jsonlite)

# function to format p values
p_value_formatter <<- function(p) { # <<- operator forces the function into the global environment.
  sapply(p, function(x) {
    if (x < 0.001) {
      return("italic(p) < 0.001")
    } else {
      return(sprintf("italic(p) == %.3f", x))
    }
  })
}

# Load static reference data
genelist1 <- readRDS("../fiber_types/data/gene_list_transcriptome.Rds")
genelist2 <- readRDS("../human_muscle_aging/data/gene_list.Rds")
genelist3 <- readRDS("../human_muscle_obesity/data/gene_list.Rds")
gene_list <- rbind(genelist1, genelist2, genelist3)
gene_list <- unique(gene_list$SYMBOL)
