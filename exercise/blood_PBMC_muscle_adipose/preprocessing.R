library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
##########################################################################################################
#
# Collect source data for shiny app 
#
##########################################################################################################
All_matrix_human <- readRDS("../../R_databases/Blood_PBMC_Muscle_Metabolic_Stress/All_matrix_human.Rds")
All_metadata_human <- readRDS("../../R_databases/Blood_PBMC_Muscle_Metabolic_Stress/All_metadata_human.Rds")

#only keep healthy
All_metadata_human <- All_metadata_human[All_metadata_human$Diagnosis %in% c("healthy", "Healthy"),]
All_matrix_human <- All_matrix_human[,All_metadata_human$sampleID]

# remove GSE83578 - blood was incuated 1h on the bench --> platelet activation
All_metadata_human <- All_metadata_human[!All_metadata_human$GEO %in% "GSE83578",]
All_matrix_human <- All_matrix_human[,All_metadata_human$sampleID]

# slip in 3 datasets
All_matrix_human_part1 <- All_matrix_human[,1:500]
All_matrix_human_part2 <- All_matrix_human[,501:1000]
All_matrix_human_part3 <- All_matrix_human[,1001:ncol(All_matrix_human)]

saveRDS(All_matrix_human_part1, "data/human_matrix_part1.Rds")
saveRDS(All_matrix_human_part2, "data/human_matrix_part2.Rds")
saveRDS(All_matrix_human_part3, "data/human_matrix_part3.Rds")
saveRDS(All_metadata_human, "data/human_metadata.Rds")

All_matrix_mouse <- readRDS("../../R_databases/Blood_PBMC_Muscle_Metabolic_Stress/All_matrix_mouse.Rds")
All_metadata_mouse <- readRDS("../../R_databases/Blood_PBMC_Muscle_Metabolic_Stress/All_metadata_mouse.Rds")
saveRDS(All_matrix_mouse, "data/mouse_matrix.Rds")
saveRDS(All_metadata_mouse, "data/mouse_metadata.Rds")


library(readxl)
path <- "../../R_databases/SingleCellRNAseq/Lovric_exercise/42003_2022_4088_MOESM8_ESM.xlsx"
excel_sheets(path = path) #get tab names

data2_endo <- data.frame(read_xlsx(path, sheet = excel_sheets(path = path)[2]),
                         cell.type = "Endothelial")
data2_lympho <- data.frame(read_xlsx(path, sheet = excel_sheets(path = path)[3]),
                           cell.type = "Lymphocyte")
data2_mesen <- data.frame(read_xlsx(path, sheet = excel_sheets(path = path)[4]),
                          cell.type = "Mesenchymal")
data2_mono <- data.frame(read_xlsx(path, sheet = excel_sheets(path = path)[5]),
                         cell.type = "Monocyte")
data2_myo <- data.frame(read_xlsx(path, sheet = excel_sheets(path = path)[6]),
                        cell.type = "Myocyte")
data2_peri <- data.frame(read_xlsx(path, sheet = excel_sheets(path = path)[7]),
                         cell.type = "Pericyte")

scRNAseq <- rbind(
  data2_endo,
  data2_lympho,
  data2_mesen,
  data2_mono,
  data2_myo,
  data2_peri
)

saveRDS(scRNAseq, "data/scRNAseq.Rds")