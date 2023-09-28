#set working directory to source file location
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

# Load libraries
library(openxlsx)
library(DT)
library(dplyr)
library(readxl)

#scrape data from order list - 2020
order_file_2020 <- read_xlsx("P:/C3_Integrative_Physiology_Group/Orders/Older files/OrderList 2020.xlsx",
                             skip = 6)
order_file_2020 <- order_file_2020[,c(3,4,7,10,13)]
colnames(order_file_2020) <- c("Product", "Product.number", "in.pack", "Ordered.for", "Date")

#scrape data from order list - 2021
order_file_2021 <- read_xlsx("P:/C3_Integrative_Physiology_Group/Orders/Older files/OrderList 2021.xlsx",
                             skip = 6)
order_file_2021 <- order_file_2021[,c(3,4,8,11,13)]
colnames(order_file_2021) <- c("Product", "Product.number", "in.pack", "Ordered.for", "Date")

#scrape data from order list - 2022
order_file_2022 <- read_xlsx("P:/C3_Integrative_Physiology_Group/Orders/Older files/OrderList 2022.xlsx",
                             skip = 6)
order_file_2022 <- order_file_2022[,c(3,4,8,11,13)]
colnames(order_file_2022) <- c("Product", "Product.number", "in.pack", "Ordered.for", "Date")

#scrape data from order list - 2023
order_file_2023 <- read_xlsx("P:/C3_Integrative_Physiology_Group/Orders/OrderList 2023.xlsx",
                             skip = 6)
order_file_2023 <- order_file_2023[,c(3,4,8,11,13)]
colnames(order_file_2023) <- c("Product", "Product.number", "in.pack", "Ordered.for", "Date")

# merge orders
order_file <- rbind(order_file_2020,
                    order_file_2021,
                    order_file_2022,
                    order_file_2023)

# keep only rows with keywords
order_taqman <- order_file[grepl("taqman", order_file$Product, ignore.case = TRUE), ]
order_taqman <- order_taqman[!grepl("mix|kit", order_taqman$Product, ignore.case = TRUE), ]

# fix dates
# replace missing dates (NA) with the preceding value
date_vec <- order_taqman$Date
date_vec
for (i in 2:length(date_vec)) {
  if (is.na(date_vec[i])) {
    date_vec[i] <- date_vec[i - 1]
  }
}
# replace the dates in the dataframe and format
order_taqman$Date <- as.Date(as.numeric(date_vec), origin = "1900-01-01")

# Find species
order_taqman$Species <- gsub(".*uman.*", "Human", order_taqman$Product)
order_taqman$Species <- gsub(".*hsa.*", "Human", order_taqman$Species)

order_taqman$Species <- gsub(".*ouse.*", "Mouse", order_taqman$Species)
order_taqman$Species <- gsub(".*mmu.*", "Mouse", order_taqman$Species)

order_taqman$Species <- gsub(".*cel.*", "C.elegans", order_taqman$Species)

# Find gene names
order_taqman$Gene <- gsub("taqman", "", order_taqman$Product, ignore.case = TRUE)
order_taqman$Gene <- gsub("human", "", order_taqman$Gene, ignore.case = TRUE)
order_taqman$Gene <- gsub("mouse", "", order_taqman$Gene, ignore.case = TRUE)
order_taqman$Gene <- gsub("rat", "", order_taqman$Gene, ignore.case = TRUE)

# fix other small things
order_taqman$in.pack <- gsub("[^0-9.-]", "", order_taqman$in.pack) %>%
  as.numeric()

# rename column names
colnames(order_taqman) <- c("order", "cat.number", "Reactions", "Name", "received.date", "Species", "Gene")

#Gets names of excel files for inventory
taqman_folder <- "P:/C3_Integrative_Physiology_Group/Reagent Lists/TaqMan_Assay_List/"

inventory_file_human <- read.xlsx(paste0(taqman_folder, "1_Human_Taqman_assays_220412.xlsx"), sheet = 1)
inventory_file_human <- inventory_file_human[,c("Box","Pos","Protein", "Gene",
                                                "CAT.#.or.accession.#" ,"Name", "Exp.date")]
colnames(inventory_file_human) <- c("Box", "Position", "Protein", "Gene", "cat.number", "Name", "Exp.date")
inventory_file_human$Species <- "Human"

inventory_file_mouse <- read.xlsx(paste0(taqman_folder, "2_Mouse_Taqman_assays_220414.xlsx"), sheet = 1)
inventory_file_mouse <- inventory_file_mouse[,c("Box#", "Pos", "Protein", "Gene",
                                                "CAT#.or.Accession.#" ,"Name", "Exp.date")]
colnames(inventory_file_mouse) <- c("Box", "Position", "Protein", "Gene", "cat.number", "Name", "Exp.date")
inventory_file_mouse$Species <- "Mouse"

inventory_file_rat <- read.xlsx(paste0(taqman_folder, "3_Rat_Taqman_assays_220419.xlsx"), sheet = 1)
inventory_file_rat <- inventory_file_rat[,c("Box","Position", "Protein", "Gene",
                                            "CAT#.or.Accession.#" ,"Name", "Exp.date")]
colnames(inventory_file_rat) <- c("Box", "Position", "Protein", "Gene", "cat.number", "Name", "Exp.date")
inventory_file_rat$Species <- "Rat"

inventory_file_hamster <- read.xlsx(paste0(taqman_folder, "4_ChineseHamster_Taqman_assays_220414.xlsx"), sheet = 1)
inventory_file_hamster <- inventory_file_hamster[,c("Box", "Pos", "Protein", "Gene",
                                                    "CAT.#.or.accession.#" ,"Name", "Exp.date")]
colnames(inventory_file_hamster) <- c("Box", "Position", "Protein", "Gene", "cat.number", "Name", "Exp.date")
inventory_file_hamster$Species <- "ChineseHamster"

inventory_file <- rbind(
  inventory_file_human,
  inventory_file_mouse,
  inventory_file_rat,
  inventory_file_hamster
)

# common boxes don't need names
inventory_file$Name <- "Common"

# merge protein and gene names
inventory_file$Gene <- paste(inventory_file$Gene, inventory_file$Protein)
inventory_file$Gene <- gsub("NA", "", inventory_file$Gene)
inventory_file$Protein <- NULL

# full inventory
full_inventory <- full_join(inventory_file,
                            order_taqman)

# order columns
full_inventory <- full_inventory[,
                                 c("Gene", "cat.number", "Species", 
                                   "Name", "Box", "Position", 
                                   "received.date", "Reactions", "Exp.date")]

# clean rows without cat.number
full_inventory <- full_inventory[!is.na(full_inventory$cat.number),]

# clean Species
table(full_inventory$Species)
full_inventory$Species[!full_inventory$Species %in% c("Human", "Mouse", "Rat", 
                                                      "ChineseHamster", "C.elegans")] <- NA
species_from_ref <- full_inventory$cat.number
species_from_ref <- gsub("Hs.*", "Human", species_from_ref, ignore.case = T)
species_from_ref <- gsub("Mm.*", "Mouse", species_from_ref, ignore.case = T)
species_from_ref <- gsub("Rn.*", "Rat", species_from_ref, ignore.case = T)
species_from_ref <- gsub("Cg.*", "ChineseHamster", species_from_ref, ignore.case = T)
table(species_from_ref)

# replace NA by the new species
full_inventory$Species <- ifelse(is.na(full_inventory$Species),
                                 species_from_ref,
                                 full_inventory$Species
                                 )
table(full_inventory$Species)

# if the same probe is present in both the common box and the order list, 
# it means that it was transferred in the common boxes during an inventory.
# Only the reference to the common box should be kept.
duplicated_probes <- full_inventory$cat.number[duplicated(full_inventory$cat.number)]
asd <- full_inventory[full_inventory$cat.number %in% duplicated_probes,]

saveRDS(full_inventory, "full_inventory.Rds")

