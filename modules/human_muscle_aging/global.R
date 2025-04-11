#-----------------------------------------------------------------------
#
# Human Aging in skeletal muscle
#
#----------------------------------------------------------------------

# Load libraries
library(shinycssloaders)
library(stringr)
library(DT)
library(plyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(feather)
library(ggforce)
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

#load metadata
gene_to_file <- readRDS("data/gene_list.Rds")
gene_list <- gene_to_file$SYMBOL
metadata <- readRDS("data/metadata.Rds")
references <- readRDS("data/references.Rds")
