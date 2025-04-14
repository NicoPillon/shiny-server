#-----------------------------------------------------------------------
#
# Profile of skeletal muscle fibers
#
#----------------------------------------------------------------------

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


# metadata
metadata_proteome <- readRDS("data/metadata_proteome.Rds")[,1:3]
metadata_transcriptome <- readRDS("data/metadata_transcriptome.Rds")

# List of genes
gene_to_file_proteome <- readRDS("data/gene_list_proteome.Rds")
gene_list_proteome <- gene_to_file_proteome$SYMBOL

gene_to_file_transcriptome <- readRDS("data/gene_list_transcriptome.Rds")
gene_list_transcriptome <- gene_to_file_transcriptome$SYMBOL

gene_list_all <- union(gene_list_proteome, gene_list_transcriptome)
