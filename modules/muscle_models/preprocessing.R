#==================================================================================================
# Preprocessing Script for Shiny App
#==================================================================================================

# Set working directory to the script's location (RStudio only)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

#==================================================================================================
# Load Required Libraries
#==================================================================================================

library(readxl)     # For reading Excel files (not used here, but possibly needed upstream)
library(arrow)      # For writing Parquet files (used for efficient data access in the app)
library(tidyverse)  # For data manipulation (dplyr, purrr, etc.)

#==================================================================================================
# SECTION 1: PROCESS AND EXPORT EXPRESSION MATRIX
#==================================================================================================

# Load full expression matrix (genes x samples)
datamatrix <- readRDS("../../rawdata/muscle_models/Data_Processed/GENENAME_norm.Rds") %>%
  data.frame()

# exclude HEK and Hela cells
datamatrix <- datamatrix[!grepl("HEK", colnames(datamatrix))]
datamatrix <- datamatrix[!grepl("HeLa", colnames(datamatrix))]

# Determine how many genes to include per chunk to divide the data in ~3 parts
n_genes <- nrow(datamatrix)
chunk_size <- ceiling(n_genes / 6)  # Number of genes per file
splits <- split(1:n_genes, ceiling(seq_along(1:n_genes) / chunk_size))  # Split row indices
file_names <- paste0("datamatrix_", seq_along(splits), ".parquet")      # Output file names

# Process each chunk: attach gene names, reorder columns, save as Parquet
gene_list <- do.call(rbind, lapply(seq_along(splits), function(i) {
  rows <- splits[[i]]
  chunk <- datamatrix[rows, ]
  
  # Add gene symbol as a dedicated column to enable filtering in `arrow::open_dataset()`
  chunk$TARGET <- rownames(chunk)
  chunk <- chunk[, c("TARGET", setdiff(colnames(chunk), "TARGET"))]  # Ensure TARGET is the first column
  
  # Write chunk as a Parquet file (efficient, columnar format)
  write_parquet(chunk, file.path("data", file_names[i]))
  
  # Return mapping: which gene is in which file
  data.frame(
    TARGET = rownames(chunk),
    file = file_names[i],
    stringsAsFactors = FALSE
  )
}))

# Save lookup table to match genes to their file location
saveRDS(gene_list, "data/list_gene.Rds")

#==================================================================================================
# SECTION 2: EXPORT METADATA
#==================================================================================================
metadata <- str_split_fixed(colnames(datamatrix), "_|\\.", 4) %>%
  data.frame()
colnames(metadata) <- c("Platform", "cell_type", "geo_accession")
metadata <- data.frame(sample = colnames(datamatrix),
                       metadata)
# rename human
metadata$cell_type <- gsub("HumanCell", "HumanPrimary", metadata$cell_type)
metadata$cell_type <- gsub("Human", "Human ", metadata$cell_type)
metadata$cell_type <- gsub("Mouse", "Mouse ", metadata$cell_type)
metadata$cell_type <- gsub("Rat", "Rat ", metadata$cell_type)

# list of samples
metadata <- metadata %>%
  select(sample, cell_type, geo_accession) %>%
  mutate(
    species = case_when(
      grepl("Human", sample, ignore.case = TRUE) ~ "Human",
      grepl("Mouse", sample, ignore.case = TRUE) ~ "Mouse",
      grepl("Rat", sample, ignore.case = TRUE) ~ "Rat",
      TRUE ~ "Other"
    ),
    cell_tissue = case_when(
      grepl("Tissue", sample, ignore.case = TRUE) ~ "Tissue",
      TRUE ~ "Cell"
    )
  )

# Save to app directory
saveRDS(metadata, file = "data/metadata.Rds")

#==================================================================================================
# SECTION 3: EXPORT REFERENCE TABLE
#==================================================================================================

# Load reference table containing study-level information (e.g. publication info)
references <- read_xlsx("../../rawdata/muscle_models/Datasets.xlsx")
references <- references[,c("GEO", "Sample", "Species", "Source", "Plateform")]

# Save to app directory
saveRDS(references, "data/references.Rds")












# library(rstudioapi)
# setwd(dirname(getActiveDocumentContext()$path))
# ##########################################################################################################
# #
# # Collect source data for shiny app 
# #
# ##########################################################################################################
# library(readxl)
# library(feather)
# library(tidyverse)
# 
# # list of datasets included
# dataset <- read_xlsx("../../../R_databases/muscle_composition_remodelling/muscle_models/Datasets.xlsx")
# dataset <- dataset[,c("GEO", "Sample", "Species", "Source", "Plateform")]
# write_feather(dataset, "data/dataset.feather")
# 
# # Load the data
# datamatrix <- readRDS("../../../R_databases/muscle_composition_remodelling/muscle_models/Data_Processed/GENENAME_norm.Rds")
# datamatrix <- datamatrix[!grepl("HEK", colnames(datamatrix))]
# datamatrix <- datamatrix[!grepl("HeLa", colnames(datamatrix))]
# write_feather(datamatrix, "data/datamatrix.feather")
# 
# # gene names
# gene_names <- rownames(datamatrix)
# saveRDS(gene_names, "data/gene_names.Rds")
# 
# # make samples
# samples_list <- str_split_fixed(colnames(datamatrix), "_|\\.", 4) %>%
#   data.frame()
# colnames(samples_list) <- c("Platform", "cell_tissue", "geo_accession")
# samples_list <- data.frame(sample = colnames(datamatrix),
#                            samples_list)
# 
# # list of samples
# samples_list <- samples_list %>%
#   select(sample, cell_tissue, geo_accession) %>%
#   mutate(
#     species = case_when(
#       grepl("Human", sample, ignore.case = TRUE) ~ "Human",
#       grepl("Mouse", sample, ignore.case = TRUE) ~ "Mouse",
#       grepl("Rat", sample, ignore.case = TRUE) ~ "Rat",
#       TRUE ~ "Other"
#     ),
#     cell_type = case_when(
#       grepl("Tissue", sample, ignore.case = TRUE) ~ "Tissue",
#       TRUE ~ "Cell"
#     )
#   )
# 
# # Define a colorblind-friendly palette (Okabe-Ito)
# okabe_ito <- c(
#   "#E69F00",     # #E69F00
#   "#56B4E9",    # #56B4E9
#   "#009E73"# #009E73
# )
# 
# # Unique species
# unique_species <- unique(samples_list$species)
# 
# # Assign colors
# samples_list$species_colors <- setNames(okabe_ito[seq_along(unique_species)], unique_species)
# 
# saveRDS(samples_list, "data/samples_list.Rds")
# 
# 
# 
# # Select relevant columns
# muscle <- datamatrix %>%
#   select(matches('HumanCell|MouseC2C12|RatL6|HumanTissue|MouseTissue|RatTissue'))
# 
# # Create a tibble with gene names as a column
# res <- muscle %>%
#   rownames_to_column("Gene") %>%
#   pivot_longer(-Gene, names_to = "Sample", values_to = "y") %>%
#   mutate(
#     x = case_when(
#       str_detect(Sample, "HumanCell")    ~ "A1",
#       str_detect(Sample, "MouseC2C12")   ~ "A2",
#       str_detect(Sample, "RatL6")        ~ "A3",
#       str_detect(Sample, "HumanTissue")  ~ "A4",
#       str_detect(Sample, "MouseTissue")  ~ "A5",
#       str_detect(Sample, "RatTissue")    ~ "A6",
#       TRUE ~ NA_character_
#     )
#   ) %>%
#   select(x, y, Gene)
# 
# # Example: filter for one gene
# res %>% filter(Gene == "SLC2A4")
# 
# 
# 
# 
# datamatrix <- readRDS("../../../R_databases/muscle_composition_remodelling/muscle_models/Data_Processed/GENENAME_norm.Rds")
# 
# muscle <- cbind(datamatrix[grepl('HumanCell', colnames(datamatrix))],
#                 datamatrix[grepl('MouseC2C12', colnames(datamatrix))],
#                 datamatrix[grepl('RatL6', colnames(datamatrix))],
#                 datamatrix[grepl('HumanTissue', colnames(datamatrix))],
#                 datamatrix[grepl('MouseTissue', colnames(datamatrix))],
#                 datamatrix[grepl('RatTissue', colnames(datamatrix))])
# 
# 
# #Make a list of data for ggplot
# res <- muscle
# x      <- c(rep('A1', length(grep('HumanCell',         colnames(muscle)))),  #list of sample types
#             rep('A2', length(grep('MouseC2C12',   colnames(muscle)))),
#             rep('A3', length(grep('RatL6',    colnames(muscle)))),
#             rep('A4', length(grep('HumanTissue',    colnames(muscle)))),
#             rep('A5', length(grep('MouseTissue', colnames(muscle)))),
#             rep('A6', length(grep('RatTissue', colnames(muscle)))))
# datalist <- vector("list", nrow(res))
# names(datalist) <- rownames(res)
# for (i in 1:nrow(res)){
#   data   <- data.frame(x=factor(), y=numeric(), Gene=character(), stringsAsFactors=FALSE) #empty dataframe to collect data
#   y     <- as.numeric(res[i,])            #collect data for gene name i
#   data <- data.frame(x, y, rep(rownames(res[i,]))) #create table with x="sample type", y="data", "gene name"
#   colnames(data) <- c("x","y","Gene")              #rename column names to make it possible to rbind later
#   datalist[[i]] <- data
# }
# 
# datalist[['SLC2A4']] #check example
# saveRDS(datalist, file="data/Muscle_Models_Profiling_data.Rds")
# 
# 
# #stats for all genes
# res <- muscle
# library(matrixStats)
# statslist <- vector("list", nrow(res))
# names(statslist) <- rownames(res)
# 
# for (i in 1:nrow(res)){
#   mean <- cbind(
#     rowMeans(res[i, grepl('HumanCell',   colnames(res))], na.rm=T),
#     rowMeans(res[i, grepl('MouseC2C12',        colnames(res))], na.rm=T),
#     rowMeans(res[i, grepl('RatL6',  colnames(res))], na.rm=T),
#     rowMeans(res[i, grepl('HumanTissue',colnames(res))], na.rm=T),
#     rowMeans(res[i, grepl('MouseTissue',   colnames(res))], na.rm=T),
#     rowMeans(res[i, grepl('RatTissue',colnames(res))], na.rm=T))
#   Sd <- cbind(
#     rowSds(as.matrix(res[i, grepl('HumanCell',   colnames(res))]), na.rm=T),
#     rowSds(as.matrix(res[i, grepl('MouseC2C12',        colnames(res))]), na.rm=T),
#     rowSds(as.matrix(res[i, grepl('RatL6',  colnames(res))]), na.rm=T),
#     rowSds(as.matrix(res[i, grepl('HumanTissue',colnames(res))]), na.rm=T),
#     rowSds(as.matrix(res[i, grepl('MouseTissue',   colnames(res))]), na.rm=T),
#     rowSds(as.matrix(res[i, grepl('RatTissue',colnames(res))]), na.rm=T))
#   nsize <- cbind(
#     rowSums(!is.na(res[i, grepl('HumanCell',   colnames(res))])),
#     rowSums(!is.na(res[i, grepl('MouseC2C12',        colnames(res))])),
#     rowSums(!is.na(res[i, grepl('RatL6',  colnames(res))])),
#     rowSums(!is.na(res[i, grepl('HumanTissue',colnames(res))])),
#     rowSums(!is.na(res[i, grepl('MouseTissue',   colnames(res))])),
#     rowSums(!is.na(res[i, grepl('RatTissue',colnames(res))])))
#   stats <- data.frame(t(mean), t(Sd), t(nsize))
#   stats[,3] <- as.factor(stats[,3])
#   colnames(stats) <- c('Mean', 'Sd', 'n')
#   rownames(stats) <- c("HumanCell", "MouseC2C12", "RatL6", 
#                        "HumanTissue", "MouseTissue", "RatTissue") 
#   statslist[[i]] <- stats
# }
# 
# statslist[['SLC2A4']] #check example
# saveRDS(statslist, file="data/Muscle_Models_Profiling_statslist.Rds")
# 
