library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)
library(feather)
library(tidyverse)

# list of datasets included
dataset <- read_xlsx("../../../R_databases/muscle_composition_remodelling/muscle_models/Datasets.xlsx")
dataset <- dataset[,c("GEO", "Sample", "Species", "Source", "Plateform")]
write_feather(dataset, "data/dataset.feather")

# Load the data
datamatrix <- readRDS("../../../R_databases/muscle_composition_remodelling/muscle_models/Data_Processed/GENENAME_norm.Rds")
datamatrix <- datamatrix[!grepl("HEK", colnames(datamatrix))]
datamatrix <- datamatrix[!grepl("HeLa", colnames(datamatrix))]
write_feather(datamatrix, "data/datamatrix.feather")

# gene names
gene_names <- rownames(datamatrix)
saveRDS(gene_names, "data/gene_names.Rds")

# make samples
samples_list <- str_split_fixed(colnames(datamatrix), "_|\\.", 4) %>%
  data.frame()
colnames(samples_list) <- c("Platform", "cell_tissue", "geo_accession")
samples_list <- data.frame(sample = colnames(datamatrix),
                           samples_list)

# list of samples
samples_list <- samples_list %>%
  select(sample, cell_tissue, geo_accession) %>%
  mutate(
    species = case_when(
      grepl("Human", sample, ignore.case = TRUE) ~ "Human",
      grepl("Mouse", sample, ignore.case = TRUE) ~ "Mouse",
      grepl("Rat", sample, ignore.case = TRUE) ~ "Rat",
      TRUE ~ "Other"
    ),
    cell_type = case_when(
      grepl("Tissue", sample, ignore.case = TRUE) ~ "Tissue",
      TRUE ~ "Cell"
    )
  )

# Define a colorblind-friendly palette (Okabe-Ito)
okabe_ito <- c(
  "#E69F00",     # #E69F00
  "#56B4E9",    # #56B4E9
  "#009E73"# #009E73
)

# Unique species
unique_species <- unique(samples_list$species)

# Assign colors
samples_list$species_colors <- setNames(okabe_ito[seq_along(unique_species)], unique_species)

saveRDS(samples_list, "data/samples_list.Rds")



# Select relevant columns
muscle <- datamatrix %>%
  select(matches('HumanCell|MouseC2C12|RatL6|HumanTissue|MouseTissue|RatTissue'))

# Create a tibble with gene names as a column
res <- muscle %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "y") %>%
  mutate(
    x = case_when(
      str_detect(Sample, "HumanCell")    ~ "A1",
      str_detect(Sample, "MouseC2C12")   ~ "A2",
      str_detect(Sample, "RatL6")        ~ "A3",
      str_detect(Sample, "HumanTissue")  ~ "A4",
      str_detect(Sample, "MouseTissue")  ~ "A5",
      str_detect(Sample, "RatTissue")    ~ "A6",
      TRUE ~ NA_character_
    )
  ) %>%
  select(x, y, Gene)

# Example: filter for one gene
res %>% filter(Gene == "SLC2A4")




datamatrix <- readRDS("../../../R_databases/muscle_composition_remodelling/muscle_models/Data_Processed/GENENAME_norm.Rds")

muscle <- cbind(datamatrix[grepl('HumanCell', colnames(datamatrix))],
                datamatrix[grepl('MouseC2C12', colnames(datamatrix))],
                datamatrix[grepl('RatL6', colnames(datamatrix))],
                datamatrix[grepl('HumanTissue', colnames(datamatrix))],
                datamatrix[grepl('MouseTissue', colnames(datamatrix))],
                datamatrix[grepl('RatTissue', colnames(datamatrix))])


#Make a list of data for ggplot
res <- muscle
x      <- c(rep('A1', length(grep('HumanCell',         colnames(muscle)))),  #list of sample types
            rep('A2', length(grep('MouseC2C12',   colnames(muscle)))),
            rep('A3', length(grep('RatL6',    colnames(muscle)))),
            rep('A4', length(grep('HumanTissue',    colnames(muscle)))),
            rep('A5', length(grep('MouseTissue', colnames(muscle)))),
            rep('A6', length(grep('RatTissue', colnames(muscle)))))
datalist <- vector("list", nrow(res))
names(datalist) <- rownames(res)
for (i in 1:nrow(res)){
  data   <- data.frame(x=factor(), y=numeric(), Gene=character(), stringsAsFactors=FALSE) #empty dataframe to collect data
  y     <- as.numeric(res[i,])            #collect data for gene name i
  data <- data.frame(x, y, rep(rownames(res[i,]))) #create table with x="sample type", y="data", "gene name"
  colnames(data) <- c("x","y","Gene")              #rename column names to make it possible to rbind later
  datalist[[i]] <- data
}

datalist[['SLC2A4']] #check example
saveRDS(datalist, file="data/Muscle_Models_Profiling_data.Rds")


#stats for all genes
res <- muscle
library(matrixStats)
statslist <- vector("list", nrow(res))
names(statslist) <- rownames(res)

for (i in 1:nrow(res)){
  mean <- cbind(
    rowMeans(res[i, grepl('HumanCell',   colnames(res))], na.rm=T),
    rowMeans(res[i, grepl('MouseC2C12',        colnames(res))], na.rm=T),
    rowMeans(res[i, grepl('RatL6',  colnames(res))], na.rm=T),
    rowMeans(res[i, grepl('HumanTissue',colnames(res))], na.rm=T),
    rowMeans(res[i, grepl('MouseTissue',   colnames(res))], na.rm=T),
    rowMeans(res[i, grepl('RatTissue',colnames(res))], na.rm=T))
  Sd <- cbind(
    rowSds(as.matrix(res[i, grepl('HumanCell',   colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[i, grepl('MouseC2C12',        colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[i, grepl('RatL6',  colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[i, grepl('HumanTissue',colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[i, grepl('MouseTissue',   colnames(res))]), na.rm=T),
    rowSds(as.matrix(res[i, grepl('RatTissue',colnames(res))]), na.rm=T))
  nsize <- cbind(
    rowSums(!is.na(res[i, grepl('HumanCell',   colnames(res))])),
    rowSums(!is.na(res[i, grepl('MouseC2C12',        colnames(res))])),
    rowSums(!is.na(res[i, grepl('RatL6',  colnames(res))])),
    rowSums(!is.na(res[i, grepl('HumanTissue',colnames(res))])),
    rowSums(!is.na(res[i, grepl('MouseTissue',   colnames(res))])),
    rowSums(!is.na(res[i, grepl('RatTissue',colnames(res))])))
  stats <- data.frame(t(mean), t(Sd), t(nsize))
  stats[,3] <- as.factor(stats[,3])
  colnames(stats) <- c('Mean', 'Sd', 'n')
  rownames(stats) <- c("HumanCell", "MouseC2C12", "RatL6", 
                       "HumanTissue", "MouseTissue", "RatTissue") 
  statslist[[i]] <- stats
}

statslist[['SLC2A4']] #check example
saveRDS(statslist, file="data/Muscle_Models_Profiling_statslist.Rds")

