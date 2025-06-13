library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
library(readxl)
library(tidyverse)

# data
VO2_data <- read_xlsx("../data_analysis/studies_all.xlsx")
colnames(VO2_data)

# select columns of interest
VO2_data <- VO2_data[,c("first_author", "publication_year", "journal", "country", "sex", "modality",
                        "age_mean", "height_cm", "weight_kg", "bmi", "HDI_group",
                        "mean", "sd", "n_size", "P01","P10","P25","P50","P75","P90","P99")]

# Calculate P40 and P60 by linear interpolation
VO2_data$P40 <- VO2_data$P25 + (VO2_data$P50 - VO2_data$P25) * (40 - 25) / (50 - 25)
VO2_data$P60 <- VO2_data$P50 + (VO2_data$P75 - VO2_data$P50) * (60 - 50) / (75 - 50)

# make nice titles
colnames(VO2_data) <- str_to_title(colnames(VO2_data))

# pivot percentles
VO2_data <- pivot_longer(VO2_data,
                         cols = c("P01","P10","P25","P40","P50","P60","P75","P90","P99"),
                         names_to = "Percentile",
                         values_to = "VO2max")
VO2_data$Percentile <- gsub("P50", "Median", VO2_data$Percentile)

# make factor to place mean in between percentiles
VO2_data$Percentile <- factor(VO2_data$Percentile,
                              levels = c("P99", "P90", "P75", "P60",
                                         "Median", 
                                         "P40", "P25", "P10", "P01"))

saveRDS(VO2_data, file = "data/data.Rds")


ggplot(VO2_data, aes(x = Age_mean, y = VO2max,
                color = Percentile, linetype = Percentile)) +
  theme_bw(16) +
  theme(panel.grid.major = element_line(color = "gray80",
                                        linewidth = 0.5,
                                        linetype = 1),
        legend.title = element_blank()) +
  scale_x_continuous(breaks = seq(0, 100, 10),
                     minor_breaks = seq(0, 100, 1)) +
  scale_y_continuous(breaks = seq(0, 100, 10),
                     minor_breaks = seq(0, 100, 1)) +
  labs(x = "Age (Years)",
       y = expression(dot("V")["O"[2]] ~ "peak (mL/min/kg)")) +
  facet_wrap(.~Sex, scales = "free_y", ncol = 2) +
  geom_smooth(se = FALSE, method = "lm", linewidth = 0.5, alpha = 0.9) +
  scale_colour_manual(values = c("#097910", "#1d7c0f", "#38800d", "#53840c",
                                 "black",
                                 "#8c8d08", "#a99106", "#c99604", "#e69a02")) +
  scale_linetype_manual(values = c(5,4,3,2,1,2,3,4,5))




# Function to create hyperlink from DOI
createLink <- function(doi) {
  if (is.na(doi)) return(NA)
  url <- paste0("https://doi.org/", doi)
  text <- doi
  sprintf('<a href="%s" target="_blank">%s</a>', url, text)
}

# make table with datasets
studies_included <- read_xlsx("../study_included.xlsx")
colnames(studies_included)
studies_included <- studies_included[,c("Authors", "Published Year", "Title", "Journal", "DOI")]

studies_included$DOI <- sapply(studies_included$DOI, createLink)

saveRDS(studies_included, file = "data/studies_included.Rds")

