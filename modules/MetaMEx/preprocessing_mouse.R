library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------------------------------
library(tidyverse)
library(feather)
library(readxl)
library(limma)

#================================================================================================================
#
# Mouse
#
#================================================================================================================
version <- "v4.240824"

# Load MetaMEx metadata
mouse_meta <- readRDS(paste0("../../../Project_MetaMEx/data_out/", version, "_metadata.Rds"))
mouse_genes <- readRDS(paste0("../MetaMEx_database/mouse/MetaMEx_mouse_", version, "_genelist.Rds"))

# load MetaMEx data
mouse_matrix <- readRDS(paste0("../MetaMEx_database/mouse/MetaMEx_mouse_", version, "_datamatrix.Rds"))

# #light dataset for optimization
mouse_genes <- mouse_genes[mouse_genes$SYMBOL %in% c("Nr4a3", "Ppargc1a", "Hes1",
                                                     "Col4a1", "znf697", "Tbxas1", "Tbxa2r",
                                                     sample(na.omit(mouse_genes$SYMBOL), 500)),]
mouse_matrix <- mouse_matrix[mouse_genes$SYMBOL,]
all(mouse_genes$SYMBOL == rownames(mouse_matrix))

# Get the number of rows in the matrix
total_rows <- nrow(mouse_matrix)
batch_size <- 1500
num_batches <- ceiling(total_rows / batch_size)

# Split and write matrix in batches
for (i in 1:num_batches) {
  start_row <- (i - 1) * batch_size + 1
  end_row <- min(i * batch_size, total_rows)
  
  batch_matrix <- mouse_matrix[start_row:end_row, ]
  write_feather(batch_matrix, paste0("data/mouse_matrix_", i, ".feather"))
}
rm(total_rows, batch_size, num_batches, i, start_row, end_row, batch_matrix)

# Write metadata
saveRDS(mouse_meta, "annotation/mouse_metadata.Rds")

# save gene list (SYMBOL)
saveRDS(mouse_genes, "annotation/mouse_genes.Rds")

#------------------------------------------------------------------------------------------------------
# If a study is acute is cannot have NA in "Acute_exercise_hours" and only "PRE" in "Training_weeks"
# If a study is training it cannot have NA in "Training_weeks" and only "PRE" in "Acute_exercise_hours"
#------------------------------------------------------------------------------------------------------
# Acute Aerobic
mouse_meta_AA <- mouse_meta[!is.na(mouse_meta$Acute_exercise_hours) &
                              !grepl("W", mouse_meta$Training_weeks),]

# Inactivity
mouse_meta_IN <- mouse_meta[mouse_meta$Study_type %in% c("Inactivity"),]

# Training Aerobic
mouse_meta_TA <- mouse_meta[!is.na(mouse_meta$Training_weeks) &
                              !grepl("H", mouse_meta$Acute_exercise_hours),]


#================================================================================================================
#
# Human - lists of conditions
#
#================================================================================================================
# Set up the different categories to be selected
unique(mouse_meta$Protocol_type) %>%
  sort()

list_categories <- list(
  GEO_AA = unique(mouse_meta_AA$GEO[order(as.numeric(sub("GSE", "", mouse_meta_AA$GEO)))]),
  GEO_IN = unique(mouse_meta_IN$GEO[order(as.numeric(sub("GSE", "", mouse_meta_IN$GEO)))]),
  GEO_TA = unique(mouse_meta_TA$GEO[order(as.numeric(sub("GSE", "", mouse_meta_TA$GEO)))]),
  
  protocol_choice = c("Acute"       = "Acute",
                      "Inactivity"  = "Inactivity",
                      "Training"    = "Training"),
  
  muscle_choice = c("Calf muscle"   = "CAF",
                    "EDL"           = "EDL",
                    "Gastrocnemius" = "GAS",
                    "Plantaris"     = "PLA",
                    "Quadriceps"    = "QUA",
                    "Quadriceps + Gastrocnemius" = "QUAGAS",
                    "Soleus"        = "SOL",
                    "Tibialis"      = "TA",
                    "Triceps"       = "TRI",
                    "Vastus Lateralis" = "VAL",
                    "Unkwnown"      = "muscle"),
  
  sex_choice = c("Male"       = "M",
                 "Female"     = "F"),
  
  age_choice = c("5 weeks"  = 5, 
                 "6 weeks"  = 6, 
                 "7 weeks"  = 7, 
                 "8 weeks"  = 8,
                 "10 weeks" = 10, 
                 "11 weeks" = 11, 
                 "12 weeks" = 12, 
                 "13 weeks" = 13, 
                 "14 weeks" = 14, 
                 "16 weeks" = 16, 
                 "19 weeks" = 19, 
                 "20 weeks" = 20, 
                 "6 months" = 24,
                 "7 months" = 28,
                 "18 months" = 72),
  
  disease_choice = c("Control/untreated" = "WT",
                     "alpha7BX2 overexpression" = "a7BX2",
                     "Bcl3-null" = "Bcl3",
                     "knock-in mutation in CaMKII" = "CaMK2gVV",
                     "Ex44 correction" = "correctedEx44",
                     "Ex44 deletion" = "deletedEx44",
                     "DNMT3a-null" = "DNMT3a",
                     "PPARδ agonist (GW1516)" = "GW1516",
                     "HDAC3-null" = "HDAC3",
                     "High fat diet" = "HFD",
                     "IL13-null" = "IL13KO",
                     "MSTN-null" = "MSTNKO",
                     "Nfkb1-null" = "Nfkb1",
                     "Nicotinamide riboside" = "NR",
                     "Nicotinamide riboside + Pterostilbene" = "NRPT",
                     "knock-in mtDNA mutator (PolG)" = "PolG",
                     "Pterostilbene" = "PT",
                     "REDD1-null" = "REDD1",
                     "Streptozotocin" = "STZ"),
  
  acute_biopsy_choice = c("Immediate" = "H00",
                          "1 hour"   = "H01",
                          "2 hours" = "H02",
                          "3 hours" = "H03",
                          "4 hours" = "H04",
                          "8 hours" = "H08",
                          "12 hours" = "H12",
                          "16 hours" = "H16",
                          "20 hours" = "H20"),
  
  inactivity_duration_choice = c("1 day" = "D01",
                                 "2 days" = "D02",
                                 "3 days" = "D03",
                                 "5 days" = "D05",
                                 "6 days" = "D06",
                                 "7 days" = "D07",
                                 "12 days" = "D12",
                                 "14 days" = "D14",
                                 "30 days" = "D30"),
  
  exercise_type_choice = c("Treadmill" = 'RUN',
                           "Treadmill, Eccentric" = 'REC'),
  
  training_protocol_choice = c("Treadmill" = 'RUN',
                               "Treadmill, Eccentric" = 'REC',
                               "Volontary Wheel" = 'WHE',
                               "Progressive weighted wheel running" = 'POW',
                               "Swimming" = 'SWI'),
  
  inactivity_protocol_choice = c("Unilateral hindlimb immobilization" = 'UHI',
                                 "Spaceflight" = 'SPF',
                                 "Hindlimb suspension" = "SUS"),
  
  training_duration_choice = c("3 weeks" = "W03",
                               "4 weeks" = "W04",
                               "5 weeks" = "W05",
                               "8 weeks" = "W08",
                               "10 weeks" = "W10",
                               "12 weeks" = "W12",
                               "24 weeks" = "W24")
)

saveRDS(list_categories, "annotation/mouse_input_categories.Rds")




#================================================================================================================
#
# Mouse - lists of all studies and their characteristics
#
#================================================================================================================
colnames(mouse_meta)
summary_table <- mouse_meta %>%
  group_by(GEO) %>%
  summarise(
    Male_Female = paste(sum(predicted_sex == "M", na.rm = TRUE), sum(Sex == "F", na.rm = TRUE), sep = " / "),
    Age = signif(mean(Age, na.rm = TRUE), 3),
    Strain = toString(unique(Strain)),
    Genotype = toString(unique(Genotype)),
    Protocols = toString(unique(Protocol_type)),
    Publication = toString(unique(Publication)),
    Platform = toString(unique(Platform)),
    Authors = toString(unique(Authors))
  ) %>%
  mutate(GEO_numeric = as.numeric(sub("GSE", "", GEO))) %>%  # Extract numeric part of GEO
  arrange(GEO_numeric) %>%  # Sort by the numeric part
  select(-GEO_numeric)  # Remove the temporary numeric column


#function to make hyperlinks
createLink <- function(val) {
  sprintf(paste0('<a href="', URLdecode(val),'" target="_blank">', gsub("(.*org/)|(.*=)", "", val) ,'</a>'))
}

# Apply function to make links
summary_table$GEO <- paste("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", summary_table$GEO, sep="")
summary_table$GEO  <- sapply(summary_table$GEO, createLink)
summary_table$Publication <- sapply(summary_table$Publication, createLink)
saveRDS(summary_table, "annotation/mouse_references.Rds")

