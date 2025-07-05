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
library(ggplot2)
library(ggiraph)
library(arrow)
library(highcharter)

#--------------------------------------------------------------------------------------------------
# Collect statistics from all apps
#--------------------------------------------------------------------------------------------------
# Need to spell out the whole web address to get out of the overview iframe.

# fiber_types
fiber2vs1_mRNA <- readRDS("../../rawdata/fiber_types/data_out/transcriptome_stats_IIvsI.Rds")
fiber2vs1_mRNA$SYMBOL <- rownames(fiber2vs1_mRNA)
fiber2vs1_mRNA$experiment <- "Fiber Type 2 vs Type 1 (mRNA)"
fiber2vs1_mRNA$model <- "Muscle Fiber"
fiber2vs1_mRNA$omics <- "Transcriptomics"
fiber2vs1_mRNA$url <- "https://shiny.nicopillon.com/modules/fiber_types.html"

fiber2vs1_protein <- readRDS("../../rawdata/fiber_types/data_out/proteome_stats_IIvsI.Rds")
fiber2vs1_protein$SYMBOL <- rownames(fiber2vs1_protein)
fiber2vs1_protein$experiment <- "Fiber Type 2 vs Type 1 (protein)"
fiber2vs1_protein$model <- "Muscle Fiber"
fiber2vs1_protein$omics <- "Proteomics"
fiber2vs1_protein$url <- "https://shiny.nicopillon.com/modules/fiber_types.html"

# human_muscle_aging
age <- readRDS("../../rawdata/human_muscle_aging/transcriptomics/data_out/AgingDataset_stats.Rds")
age$SYMBOL <- rownames(age)
age$experiment <- "Old vs Young Age (mRNA)"
age$model <- "Muscle Tissue"
age$omics <- "Transcriptomics"
age$url <- "https://shiny.nicopillon.com/modules/human_muscle_aging.html"

# human_muscle_obesity
obesity_mRNA <- readRDS("../../rawdata/human_muscle_obesity/transcriptomics/data_out/stats.Rds")
obesity_mRNA$SYMBOL <- rownames(obesity_mRNA)
obesity_mRNA$experiment <- "Obesity vs Lean Weight (mRNA)"
obesity_mRNA$model <- "Muscle Tissue"
obesity_mRNA$omics <- "Transcriptomics"
obesity_mRNA$url <- "https://shiny.nicopillon.com/modules/human_muscle_obesity.html"

# MetaMEx
exercise.acute_mRNA <- read.csv("../../../Project_MetaMEx/3_MetaMEx_statistics/human_limma_AA.csv", row.names = 1)
exercise.acute_mRNA$SYMBOL <- rownames(exercise.acute_mRNA)
exercise.acute_mRNA$experiment <- "Acute Exercise (mRNA)"
exercise.acute_mRNA$model <- "Muscle Tissue"
exercise.acute_mRNA$omics <- "Transcriptomics"
exercise.acute_mRNA$url <- "https://shiny.nicopillon.com/modules/MetaMEx.html"

exercise.training_mRNA <- read.csv("../../../Project_MetaMEx/3_MetaMEx_statistics/human_limma_TA.csv", row.names = 1)
exercise.training_mRNA$SYMBOL <- rownames(exercise.training_mRNA)
exercise.training_mRNA$experiment <- "Exercise Training (mRNA)"
exercise.training_mRNA$model <- "Muscle Tissue"
exercise.training_mRNA$omics <- "Transcriptomics"
exercise.training_mRNA$url <- "https://shiny.nicopillon.com/modules/MetaMEx.html"

inactivity_mRNA <- read.csv("../../../Project_MetaMEx/3_MetaMEx_statistics/human_limma_IN.csv", row.names = 1)
inactivity_mRNA$SYMBOL <- rownames(inactivity_mRNA)
inactivity_mRNA$experiment <- "Inactivity (mRNA)"
inactivity_mRNA$model <- "Muscle Tissue"
inactivity_mRNA$omics <- "Transcriptomics"
inactivity_mRNA$url <- "https://shiny.nicopillon.com/modules/MetaMEx.html"

# muscle_clamp
clamp <- readRDS("../../rawdata/muscle_clamp/data_out/stats.Rds")
clamp$SYMBOL <- rownames(clamp)
clamp$experiment <- "Hyperinsulinemic Euglycemic Clamp (mRNA)"
clamp$model <- "Muscle Tissue"
clamp$omics <- "Transcriptomics"
clamp$url <- "https://shiny.nicopillon.com/modules/muscle_clamp.html"

# muscle_myocyte_sex
sex_mRNA <- readRDS("../../rawdata/muscle_myocyte_sex/human/data_out/stats.Rds")
colnames(sex_mRNA)[9] <- "SYMBOL"
sex_mRNA$experiment <- "Male vs Female (mRNA)"
sex_mRNA$model <- "Muscle Tissue"
sex_mRNA$omics <- "Transcriptomics"
sex_mRNA$url <- "https://shiny.nicopillon.com/modules/muscle_myocyte_sex.html"

# myotube_EPS
EPS <- readRDS("../../rawdata/myotube_EPS/data_out/stats.Rds")
EPS$SYMBOL <- rownames(EPS)
EPS$experiment <- "Myotubes Exposed to Electrical Pulse Stimulation (mRNA)"
EPS$model <- "Myotube"
EPS$omics <- "Transcriptomics"
EPS$url <- "https://shiny.nicopillon.com/modules/myotube_EPS.html"

# myotube_palmitate
palmitate <- readRDS("../../rawdata/myotube_palmitate/data_out/stats.Rds")
palmitate$SYMBOL <- rownames(palmitate)
palmitate$experiment <- "Myotubes Exposed to Palmitate (mRNA)"
palmitate$model <- "Myotube"
palmitate$omics <- "Transcriptomics"
palmitate$url <- "https://shiny.nicopillon.com/modules/myotube_palmitate.html"


# merge all 
datamatrix <- rbind(
  fiber2vs1_mRNA, fiber2vs1_protein,
  age,
  obesity_mRNA,
  clamp,
  sex_mRNA,
  palmitate,
  EPS,
  exercise.acute_mRNA, exercise.training_mRNA, inactivity_mRNA
)
rownames(datamatrix) <- gsub(" ", "_", paste(datamatrix$SYMBOL, datamatrix$experiment))

# Determine how many genes to include per chunk to divide the data in ~3 parts
n_genes <- nrow(datamatrix)
chunk_size <- ceiling(n_genes / 3)  # Number of genes per file
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
gene_list$SYMBOL <- gsub("_.*", "", gene_list$TARGET)
saveRDS(gene_list, "data/list_gene.Rds")


###################################################################################################

#dat <- selected_row
dat <- datamatrix[datamatrix$SYMBOL == "PDK4", ]  # Keep as data.frame

# Add Significance category
dat <- dat %>%
  mutate(Significance = case_when(
    adj.P.Val < 0.001 ~ "FDR < 0.001",
    adj.P.Val < 0.01  ~ "FDR < 0.01",
    adj.P.Val < 0.05  ~ "FDR < 0.05",
    TRUE              ~ "ns"
  )) %>%
  mutate(Significance = factor(Significance, levels = c("FDR < 0.001", "FDR < 0.01", "FDR < 0.05", "ns")))

# Sort and preserve order
dat <- dat %>%
  arrange(desc(logFC)) %>%
  mutate(experiment = factor(experiment, levels = experiment))

# Define color mapping for significance
signif_colors <- c(
  "FDR < 0.001" = "darkgreen",
  "FDR < 0.01"  = "orange",
  "FDR < 0.05"  = "yellow",
  "ns"          = "gray"
)

# Convert to list of points for highcharter
data_points <- purrr::pmap(dat, function(experiment, logFC, adj.P.Val, url, Significance, model, omics,...) {
  list(
    y = round(logFC, 2),
    name = experiment,
    omics = omics,
    color = signif_colors[[Significance]],
    url = url,
    custom = list(
      model = model,
      logFC = round(logFC, 2),
      FDR = signif(adj.P.Val, 2)
    )
  )
})

hc_theme_custom <- hc_theme(
  chart = list(
    style = list(fontFamily = "Arial")
  ),
  title = list(
    style = list(fontFamily = "Arial")
  )
)

highchart() %>%
  hc_add_theme(hc_theme_custom) %>%
  hc_chart(
    type = "bar",
    plotBorderWidth = 1,
    plotBorderColor = NULL
  ) %>%
  hc_xAxis(
    type = "category",
    gridLineWidth = 1,
    gridLineColor = "#e0e0e0",
    labels = list(style = list(color = "black", 
                               fontSize = "14px",
                               whiteSpace = "nowrap",
                               textOverflow = "none",
                               overflow = "allow"))
  ) %>%
  hc_yAxis(
    title = list(text = "logFC", style = list(color = "black", fontSize = "14px")),
    labels = list(style = list(color = "gray", fontSize = "12px")),
    gridLineWidth = 1,
    gridLineColor = "#e0e0e0",
    plotLines = list(
      list(value = 0, width = 2, color = "black", zIndex = 5)
    )
  ) %>%
  hc_add_series(
    name = "Experiments",
    data = data_points,
    showInLegend = FALSE
  ) %>%
  hc_plotOptions(series = list(
    cursor = "pointer",
    point = list(events = list(
      click = JS("function() { window.open(this.url, '_blank'); }")
    ))
  )) %>%
  hc_tooltip(
    useHTML = TRUE,
    pointFormat = paste0(
      "<b>OMICS:</b> {point.omics}<br>",
      "<b>Model:</b> {point.custom.model}<br>",
      "<b>logFC:</b> {point.custom.logFC}<br>",
      "<b>FDR:</b> {point.custom.FDR}<br>",
      "<i>Click to explore this dataset!</i>"
    )
  )

