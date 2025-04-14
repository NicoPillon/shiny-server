#--------------------------------------------------------------------------------------------------------
#
# App heavy template for testing load time
#
#--------------------------------------------------------------------------------------------------------

# set working directory to source file location
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

# load libraries
library(readxl)
library(feather)
library(tidyverse)

# make a larger list of genes
gene_list <- paste0("GENE", 1:5000)  # ⬅️ from 100 → 10,000
saveRDS(gene_list, file = "data/gene_list.Rds")

# make a larger list of samples
sample_list <- tibble(
  sample_id = paste0("Sample_", 1:1000),  # ⬅️ from 20 → 1000
  condition = rep(c("Control", "Treatment"), each = 500),
  species = rep(c("Human", "Mouse", "Rat", "Dog", "Pig"), length.out = 1000)
)
saveRDS(sample_list, file = "data/sample_list.Rds")

# Make table with more datasets
references <- tibble(
  Dataset = paste0("Dataset", 1:20),
  Species = rep(c("Human", "Mouse", "Rat"), length.out = 20),
  Platform = paste0("Platform_", sample(1:5, 20, replace = TRUE))
)
write_feather(references, "data/dataset.feather")

# Export the large expression matrix in feather format
n_genes <- length(gene_list)
n_samples <- nrow(sample_list)

set.seed(123)  # for reproducibility
expression_matrix <- matrix(
  rnorm(n_genes * n_samples, mean = 8, sd = 2),  # ~80 million values
  nrow = n_genes,
  ncol = n_samples,
  dimnames = list(gene_list, sample_list$sample_id)
) %>%
  data.frame()

write_feather(expression_matrix, "data/datamatrix.feather")
