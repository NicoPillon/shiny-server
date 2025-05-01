#-----------------------------------------------------------------------
#
# Human Obesity
#
#----------------------------------------------------------------------
# Load libraries
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

#load lists to locate files
gene_to_file <- readRDS("data/list_gene.Rds")
metabolite_to_file <- readRDS("data/list_metabolite.Rds")

target_to_file <- rbind(gene_to_file,
                        metabolite_to_file)
target_list <- c(target_to_file$TARGET)

# metadata
metadata <- readRDS("data/metadata.Rds")

# references
references <- readRDS("data/references.Rds")
