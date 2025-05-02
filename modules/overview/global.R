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
library(shinyjs)
library(r3dmol)

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
genelist3 <- readRDS("../human_muscle_obesity/data/list_gene.Rds")
gene_list <- rbind(genelist1, genelist2)
gene_list <- unique(gene_list$SYMBOL)



# library(biomaRt)
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# gene_map <- getBM(
#   attributes = c("hgnc_symbol", "uniprot_gn_id"),
#   mart = ensembl
# )
# 
# # Filter out empty IDs
# gene_map <- gene_map[gene_map$uniprot_gn_id != "", ]
# 
# # Save for later use
# saveRDS(gene_map, "gene_uniprot_map.rds")
gene_protein_map <- readRDS("gene_uniprot_map.rds")


#----------------------------------------------------------------------------------------------
# MetaMEx overview
#----------------------------------------------------------------------------------------------
#MetaMEx_human <- readRDS("../../../Project_MetaMEx/limma/limma_results.Rds")
#saveRDS(MetaMEx_human, "MetaMEx_human_limma_results.Rds")
# load all data - Human
MetaMEx_human <- readRDS("MetaMEx_human_limma_results.Rds")
