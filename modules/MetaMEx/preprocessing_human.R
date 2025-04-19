library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------------------------------
library(tidyverse)
library(feather)
library(readxl)
library(limma)

#================================================================================================================
#
# Human
#
#================================================================================================================
version <- "v4.250110"

# Load MetaMEx metadata
human_meta <- readRDS(paste0("../../../Project_MetaMEx/data_out/MetaMEx_human_", version, "_metadata_withMLpredictions.Rds"))
human_genes <- readRDS(paste0("../../../Project_MetaMEx/data_out/MetaMEx_human_", version, "_genelist.Rds"))

# add groups
human_meta$predicted_ageGroup <- "18.34"
human_meta$predicted_ageGroup[human_meta$predicted_age > 34] <- "35.64"
human_meta$predicted_ageGroup[human_meta$predicted_age > 64] <- "65.90"

human_meta$predicted_BMIgroup <- "Lean"
human_meta$predicted_BMIgroup[human_meta$predicted_BMI >24.999] <- "Overweight"
human_meta$predicted_BMIgroup[human_meta$predicted_BMI >29.999] <- "Obesity"


# load MetaMEx data
human_matrix <- readRDS(paste0("../../../Project_MetaMEx/data_out/MetaMEx_human_", version, "_datamatrixCleaned.Rds"))

# #light dataset for optimization
human_genes <- human_genes[human_genes$SYMBOL %in% c("NR4A3", "PPARGC1A", "HES1",
                                                     "COL4A1", "ZNF697", "TBXAS1", "TBXA2R",
                                                     sample(na.omit(human_genes$SYMBOL), 100)),]
human_matrix <- human_matrix[human_genes$SYMBOL,]
all(human_genes$SYMBOL == rownames(human_matrix))

# Get the number of rows in the matrix
total_rows <- nrow(human_matrix)
batch_size <- 1500
num_batches <- ceiling(total_rows / batch_size)

# Split and write matrix in batches
for (i in 1:num_batches) {
  start_row <- (i - 1) * batch_size + 1
  end_row <- min(i * batch_size, total_rows)
  
  batch_matrix <- human_matrix[start_row:end_row, ]
  write_feather(batch_matrix, paste0("data/human_matrix_", i, ".feather"))
}
rm(total_rows, batch_size, num_batches, i, start_row, end_row, batch_matrix)

# Write metadata
saveRDS(human_meta, "annotation/human_metadata.Rds")

# save gene list (SYMBOL)
saveRDS(human_genes, "annotation/human_genes.Rds")

#------------------------------------------------------------------------------------------------------
# If a study is acute is cannot have NA in "Acute_exercise_hours" and only "PRE" in "Training_weeks"
# If a study is training it cannot have NA in "Training_weeks" and only "PRE" in "Acute_exercise_hours"
#------------------------------------------------------------------------------------------------------
# Acute Aerobic
human_meta_AA <- human_meta[!is.na(human_meta$Acute_exercise_hours) &
                                    !grepl("W", human_meta$Training_weeks) &
                                    human_meta$Protocol_group %in% c("Basal", "Aerobic"),]

# Acute HIT
human_meta_AH <- human_meta[!is.na(human_meta$Acute_exercise_hours) &
                                    !grepl("W", human_meta$Training_weeks) &
                                    human_meta$Protocol_group %in% c("Basal", "HIT"),]

# Acute Resistance
human_meta_AR <- human_meta[!is.na(human_meta$Acute_exercise_hours) &
                                    !grepl("W", human_meta$Training_weeks) &
                                    human_meta$Protocol_group %in% c("Basal", "Resistance"),]

# Inactivity
human_meta_IN <- human_meta[human_meta$Study_group %in% c("Inactivity"),]

# Training Aerobic
human_meta_TA <- human_meta[!is.na(human_meta$Training_weeks) &
                                    !grepl("H", human_meta$Acute_exercise_hours) &
                                    human_meta$Protocol_group %in% c("Basal", "Aerobic"),]

# Training Combined
human_meta_TC <- human_meta[!is.na(human_meta$Training_weeks) &
                                    !grepl("H", human_meta$Acute_exercise_hours) &
                                    human_meta$Protocol_group %in% c("Basal", "Combined"),]

# Training HIT
human_meta_TH <- human_meta[!is.na(human_meta$Training_weeks) &
                                    !grepl("H", human_meta$Acute_exercise_hours) &
                                    human_meta$Protocol_group %in% c("Basal", "HIT"),]

# Training Resistance
human_meta_TR <- human_meta[!is.na(human_meta$Training_weeks) &
                                    !grepl("H", human_meta$Acute_exercise_hours) &
                                    human_meta$Protocol_group %in% c("Basal", "Resistance"),]


#================================================================================================================
#
# Human - lists of conditions
#
#================================================================================================================
# Set up the different categories to be selected
unique(human_meta$Health_group) %>%
  sort()

list_categories <- list(
  GEO_AA = unique(human_meta_AA$GEO[order(as.numeric(sub("GSE", "", human_meta_AA$GEO)))]),
  GEO_AH = unique(human_meta_AH$GEO[order(as.numeric(sub("GSE", "", human_meta_AH$GEO)))]),
  GEO_AR = unique(human_meta_AR$GEO[order(as.numeric(sub("GSE", "", human_meta_AR$GEO)))]),
  GEO_IN = unique(human_meta_IN$GEO[order(as.numeric(sub("GSE", "", human_meta_IN$GEO)))]),
  GEO_TA = unique(human_meta_TA$GEO[order(as.numeric(sub("GSE", "", human_meta_TA$GEO)))]),
  GEO_TC = unique(human_meta_TC$GEO[order(as.numeric(sub("GSE", "", human_meta_TC$GEO)))]),
  GEO_TH = unique(human_meta_TH$GEO[order(as.numeric(sub("GSE", "", human_meta_TH$GEO)))]),
  GEO_TR = unique(human_meta_TR$GEO[order(as.numeric(sub("GSE", "", human_meta_TR$GEO)))]),
  
  protocol_choice = c("Acute Aerobic"       = "Acute Aerobic",
                      "Acute HIT"           = "Acute HIT",
                      "Acute Resistance"    = "Acute Resistance",
                      "Inactivity"          = "Inactivity",
                      "Aerobic Training"    = "Training Aerobic",
                      "Resistance Training" = "Training Resistance",
                      "Combined Training"   = "Training Combined",
                      "HIT training"        = "Training HIT"),

  muscle_choice = c("Vastus Lateralis" = "Vastus lateralis",
                    "Biceps Brachii"   = "Biceps brachii",
                    "Soleus"           = "Soleus",
                    "Gastrocnemius"    = "Gastrocnemius"),
  
  sex_choice = c("Male"       = "M",
                 "Female"     = "F"),
  
  age_choice = c("< 35"   = "18.34",
                 "35-64"  = "35.64",
                 "> 65" = "65.90"),
  
  training_choice = c("Sedentary" = "SED",
                      "Active"    = "ACT",
                      "Athlete"   = "ATH"),
  
  obesity_choice = c("Lean" = "Lean",
                     "Overweight" = "Overweight",
                     "Obesity" = "Obesity"),
  
  disease_choice = c("Bone fracture",
                     "Chronic Kidney Disease" = "Chronic Kidney Disease",
                     "Chronic Obstructive Pulmonary Disease" = "COPD",
                     "Dyslipidemia" = "Dyslipidemia",
                     "Frailty" = "Frailty",
                     "Healthy" = "Healthy",
                     "Heart failure" = "Heart failure",
                     "Metabolic Syndrome" = "Metabolic syndrome",
                     "Myopathy" = "Myopathy",
                     "Parkinson's disease" = "Parkinson's disease",
                     "Peripheral Arterial Disease" = "Peripheral arterial disease",
                     "Sarcopenia" = "Sarcopenia",
                     "Type 1 diabetes" = "Type 1 diabetes",
                     "Type 2 diabetes" = "Type 2 diabetes"),
  
  exercise_type_choice = c("Concentric" = "CON",
                           "Eccentric"  = "ECC",
                           "Mixed"      = "MIX"),
  # 
  # acute_biopsy_choice = c("Immediate" = "H00",
  #                         "1 hour"    = "H01",
  #                         "3 hours" = "H03",
  #                         "4 hours" = "H04",
  #                         "5 hours" = "H05",
  #                         "6 hours" = "H06",
  #                         "8 hours" = "H08",
  #                         "18 hours" = "H18",
  #                         "24 hours" = "H24",
  #                         "48 hours" = "H48",
  #                         "96 hours" = "H96"),
  
  # inactivity_duration_choice = c("2 days" = "D02",
  #                                "5 days" = "D05",
  #                                "10 days" = "D10",
  #                                "14 days" = "D14",
  #                                "21 days" = "D21",
  #                                "35 days" = "D35",
  #                                "60 days" = "D60",
  #                                "70 days" = "D70",
  #                                "84 days" = "D84"),
  
  inactivity_protocol_choice = c("Bed Rest" = 'Bed rest',
                                 "Cast" = "Cast",
                                 "Limb Immobilization" = 'Immobilization',
                                 "Suspension" = "ULLS"),
  
  # training_duration_choice = c("1 week" = "W01",
  #                              "2 weeks" = "W02",
  #                              "3 weeks" = "W03",
  #                              "6 weeks" = "W06",
  #                              "8 weeks" = "W08",
  #                              "10 weeks" = "W10",
  #                              "12 weeks" = "W12",
  #                              "14 weeks" = "W14",
  #                              "16 weeks" = "W16",
  #                              "18 weeks" = "W18",
  #                              "20 weeks" = "W20",
  #                              "24 weeks" = "W24",
  #                              "26 weeks" = "W26",
  #                              "32 weeks" = "W32",
  #                              "36 weeks" = "W36",
  #                              "52 weeks" = "W52",
  #                              "60 weeks" = "W60",
  #                              "4 years"  = "W208",
  #                              "8 years"  = "W416",
  #                              "Lifelong" = "W999"),
  
  training_biopsy_choice = c("3 hours" = 2.5,
                             "24 hours" = 24,
                             "48 hours" = 48,
                             "72 hours" = 72,
                             "96 hours" = 96,
                             "N.A." = NA)
)

saveRDS(list_categories, "annotation/human_input_categories.Rds")


#================================================================================================================
#
# Human - Timeline acute exercise
#
#================================================================================================================
acute_metadata <- human_meta[!is.na(human_meta$Acute_exercise_hours),]
acute_metadata <- acute_metadata[!grepl("W", acute_metadata$Training_weeks),]

# keep only Healthy
acute_metadata <- acute_metadata[acute_metadata$Health_group %in% "Healthy",]
acute_metadata <- acute_metadata[acute_metadata$predicted_BMI < 30,]

# make bins by days
acute_metadata$bins <- "PRE"
acute_metadata$bins[acute_metadata$Acute_exercise_hours %in% c("H00", "H01")] <- "0 - 1"
acute_metadata$bins[acute_metadata$Acute_exercise_hours %in% c("H03", "H04")] <- "3 - 4"
acute_metadata$bins[acute_metadata$Acute_exercise_hours %in% c("H05", "H06", "H08")] <- "5 - 8"
acute_metadata$bins[acute_metadata$Acute_exercise_hours %in% c("H18", "H24")] <- "18 - 24"
acute_metadata$bins[acute_metadata$Acute_exercise_hours %in% c("H48", "H96")] <- "48 - 96"
acute_metadata$bins <- factor(acute_metadata$bins, levels = c("PRE", "0 - 1", "3 - 4", "5 - 8", "18 - 24", "48 - 96"))
table(acute_metadata$bins)

# Select the data matrix
res <- human_matrix[,acute_metadata$SampleID]

# Combine metadata with expression data
combined_data <- as.data.frame(t(res))
combined_data$GEO <- acute_metadata$GEO
combined_data$bins <- acute_metadata$bins
combined_data$sex <- acute_metadata$predicted_sex

# Average the data by GEO, bins, and sex
averaged_data <- combined_data %>%
  group_by(GEO, bins, sex) %>%
  summarize(across(everything(), mean), .groups = 'drop')

# Convert the averaged data back to a matrix form for further analysis
averaged_matrix <- as.matrix(t(averaged_data[, -(1:3)]))  # Exclude GEO, bins, and sex columns

# Prepare the metadata to match the averaged matrix
averaged_metadata <- acute_metadata %>%
  group_by(GEO, bins, predicted_sex) %>%
  summarize(
    across(
      everything(),
      ~ if(length(unique(.)) == 1) unique(.) else NA,
      .names = "avg_{col}"
    ),
    .groups = 'drop'
  )

# make linear model, block for confounders
design <- model.matrix(~factor(averaged_metadata$bins)
                       + factor(averaged_metadata$GEO)
                       + factor(averaged_metadata$predicted_sex))
fit <- lmFit(averaged_matrix, design)
fit2 <- eBayes(fit)
summary(decideTests(fit), adjust.method = "bonferroni", p.value = 0.05)[,1:6]

# Statistics - F test
Ftest <- data.frame(topTable(fit2, coef=2:6, adjust.method="bonferroni",n=Inf, confint=T, sort.by="none"))

# Statistics - individual comparisons
TTest.0_1 <- data.frame(topTable(fit2, coef=2, adjust.method="bonferroni",n=Inf, confint=T, sort.by="none"))
hist(TTest.0_1$P.Value)
TTest.3_4 <- data.frame(topTable(fit2, coef=3, adjust.method="bonferroni",n=Inf, confint=T, sort.by="none"))
hist(TTest.3_4$P.Value)
TTest.5_8 <- data.frame(topTable(fit2, coef=4, adjust.method="bonferroni",n=Inf, confint=T, sort.by="none"))
hist(TTest.5_8$P.Value)
TTest.18_24 <- data.frame(topTable(fit2, coef=5, adjust.method="bonferroni",n=Inf, confint=T, sort.by="none"))
hist(TTest.18_24$P.Value)
TTest.48_96 <- data.frame(topTable(fit2, coef=6, adjust.method="bonferroni",n=Inf, confint=T, sort.by="none"))
hist(TTest.48_96$P.Value)

# merge all stats
timeline_acute_stats <- data.frame(SYMBOL = rownames(averaged_matrix),
                                   Ftest.F = Ftest$F,
                                   Ftest.P.Value = Ftest$P.Value,
                                   Ftest.adj.P.Val = Ftest$adj.P.Val,
                                   TTest.0_1$logFC,
                                   TTest.3_4$logFC,
                                   TTest.5_8$logFC,
                                   TTest.18_24$logFC,
                                   TTest.48_96$logFC,
                                   TTest.0_1$P.Value,
                                   TTest.3_4$P.Value,
                                   TTest.5_8$P.Value,
                                   TTest.18_24$P.Value,
                                   TTest.48_96$P.Value,
                                   TTest.0_1$adj.P.Val,
                                   TTest.3_4$adj.P.Val,
                                   TTest.5_8$adj.P.Val,
                                   TTest.18_24$adj.P.Val,
                                   TTest.48_96$adj.P.Val)
timeline_acute_stats[timeline_acute_stats$SYMBOL == "NR4A3",]

#save
write_feather(timeline_acute_stats, "data/human_timelineStats_acute.feather")



#================================================================================================================
#
# Human - Timeline inactivity
#
#================================================================================================================
inactivity_metadata <- human_meta[!is.na(human_meta$Inactivity_days),]
inactivity_metadata <- inactivity_metadata[!grepl("RELOAD", inactivity_metadata$Inactivity_days),]

# keep only Healthy
inactivity_metadata <- inactivity_metadata[inactivity_metadata$Health_group %in% "Healthy",]
inactivity_metadata <- inactivity_metadata[inactivity_metadata$predicted_BMI < 30,]

# make bins by weeks
inactivity_metadata$bins <- "PRE"
inactivity_metadata$bins[inactivity_metadata$Inactivity_days %in% c("D02", "D05")] <- "< 1"
inactivity_metadata$bins[inactivity_metadata$Inactivity_days %in% c("D10", "D14", "D21")] <- "1 - 4"
inactivity_metadata$bins[inactivity_metadata$Inactivity_days %in% c("D35", "D60", "D70", "D84")] <- "> 5"
inactivity_metadata$bins <- factor(inactivity_metadata$bins, levels = c("PRE", "< 1", "1 - 4", "> 5"))
table(inactivity_metadata$bins)

# select matrix
res <- human_matrix[,inactivity_metadata$SampleID]

# Combine metadata with expression data
combined_data <- as.data.frame(t(res))
combined_data$GEO <- inactivity_metadata$GEO
combined_data$bins <- inactivity_metadata$bins
combined_data$sex <- inactivity_metadata$predicted_sex

# Average the data by GEO, bins, and sex
averaged_data <- combined_data %>%
  group_by(GEO, bins, sex) %>%
  summarize(across(everything(), mean), .groups = 'drop')

# Convert the averaged data back to a matrix form for further analysis
averaged_matrix <- as.matrix(t(averaged_data[, -(1:3)]))  # Exclude GEO, bins, and sex columns

# Prepare the metadata to match the averaged matrix
averaged_metadata <- inactivity_metadata %>%
  group_by(GEO, bins, predicted_sex) %>%
  summarize(
    across(
      everything(),
      ~ if(length(unique(.)) == 1) unique(.) else NA,
      .names = "avg_{col}"
    ),
    .groups = 'drop'
  )

# make linear model, block for confounders
design <- model.matrix(~factor(averaged_metadata$bins)
                       + factor(averaged_metadata$GEO)
                       + factor(averaged_metadata$predicted_sex))
fit <- lmFit(averaged_matrix, design)
fit2 <- eBayes(fit)
summary(decideTests(fit), adjust.method = "bonferroni", p.value = 0.05)[,1:6]

# Statistics - F test
Ftest <- data.frame(topTable(fit2, coef=2:4, adjust.method="bonferroni",n=Inf, confint=T, sort.by="none"))

# Statistics - individual comparisons
TTest.0_1 <- data.frame(topTable(fit2, coef=2, adjust.method="bonferroni",n=Inf, confint=T, sort.by="none"))
hist(TTest.0_1$P.Value)
TTest.1_4 <- data.frame(topTable(fit2, coef=3, adjust.method="bonferroni",n=Inf, confint=T, sort.by="none"))
hist(TTest.1_4$P.Value)
TTest.5_12 <- data.frame(topTable(fit2, coef=4, adjust.method="bonferroni",n=Inf, confint=T, sort.by="none"))
hist(TTest.5_12$P.Value)

# merge all stats
timeline_inactivity_stats <- data.frame(SYMBOL = rownames(res),
                                   Ftest.F = Ftest$F,
                                   Ftest.P.Value = Ftest$P.Value,
                                   Ftest.adj.P.Val = Ftest$adj.P.Val,
                                   TTest.0_1$logFC,
                                   TTest.1_4$logFC,
                                   TTest.5_12$logFC,
                                   TTest.0_1$P.Value,
                                   TTest.1_4$P.Value,
                                   TTest.5_12$P.Value,
                                   TTest.0_1$adj.P.Val,
                                   TTest.1_4$adj.P.Val,
                                   TTest.5_12$adj.P.Val)
timeline_inactivity_stats[timeline_inactivity_stats$SYMBOL == "NR4A3",]

#save
write_feather(timeline_inactivity_stats, "data/human_timelineStats_inactivity.feather")



#================================================================================================================
#
# Human - lists of all studies and their characteristics
#
#================================================================================================================
summary_table <- human_meta %>%
  group_by(GEO) %>%
  summarise(
    Male_Female = paste(sum(predicted_sex == "M", na.rm = TRUE), sum(Sex == "F", na.rm = TRUE), sep = " / "),
    Age = signif(mean(predicted_age, na.rm = TRUE), 3),
    BMI = signif(mean(predicted_BMI, na.rm = TRUE), 3),
    Diagnosis = toString(unique(Health_group)),
    Protocols = toString(unique(Protocol_group)),
    Publication = toString(unique(Publication)),
    Platform = toString(unique(Platform)),
    Authors = toString(unique(Authors))
  ) %>%
  mutate(GEO_numeric = as.numeric(sub("GSE", "", GEO))) %>%  # Extract numeric part of GEO
  arrange(GEO_numeric) %>%  # Sort by the numeric part
  dplyr::select(-GEO_numeric)  # Remove the temporary numeric column


#function to make hyperlinks
createLink <- function(val) {
  sprintf(paste0('<a href="', URLdecode(val),'" target="_blank">', gsub("(.*org/)|(.*=)", "", val) ,'</a>'))
}

# Apply function to make links
summary_table$GEO <- paste("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", summary_table$GEO, sep="")
summary_table$GEO  <- sapply(summary_table$GEO, createLink)
summary_table$Publication <- sapply(summary_table$Publication, createLink)
saveRDS(summary_table, "annotation/human_references.Rds")

