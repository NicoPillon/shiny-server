#--------------------------------------------------------------------------------------------------------
#
# App short description
#
#--------------------------------------------------------------------------------------------------------

# set working directory to source file location
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

# load libraries
library(readxl)
library(feather)
library(tidyverse)

# make a list of genes/proteins/metabolites
gene_list <- paste0("GENE", 1:100)
saveRDS(gene_list, file = "data/gene_list.Rds")


# make a list of samples
sample_list <- tibble(
  sample_id = paste0("Sample_", 1:20),
  condition = rep(c("Control", "Treatment"), each = 10),
  species = rep(c("Human", "Mouse"), length.out = 20)
)
saveRDS(sample_list, file = "data/sample_list.Rds")


# Make table with datasets included in the analysis
references <- data.frame(
  Dataset = paste0("Dataset", 1:10),
  Species = c(rep("Human", 5), rep("Mouse", 5)),
  Plateform = paste0("platform", 1:10)
)
write_feather(references, "data/dataset.feather")

# Export the datamatrix in feather format
n_genes <- length(gene_list)
n_samples <- nrow(sample_list)
expression_matrix <- matrix(
  rnorm(n_genes * n_samples, mean = 8, sd = 1),  # Simulated log2 expression
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(gene_list, sample_list$sample_id)
) %>%
  data.frame()
write_feather(expression_matrix, "data/datamatrix.feather")
