#set working directory to source file location
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

# Load libraries
library(openxlsx)
library(dplyr)
library(readxl)

#scrape data from order list - 2025
order_file_2025 <- read_xlsx("../../000_IntFys_documents/Order_lists/OrderList 2025.xlsx",
                             skip = 6)
order_file_2025 <- order_file_2025[,c(3,4,7,10,11)]
colnames(order_file_2025) <- c("Product", "Product.number", "in.pack", "Ordered.for", "Date")


# merge orders
order_file <- rbind(order_file_2025)

rm(order_file_2025)

# keep only rows with keywords
order_taqman <- order_file[grepl("taqman", order_file$Product, ignore.case = TRUE), ]
order_taqman <- order_taqman[!grepl("mix|kit|array", order_taqman$Product, ignore.case = TRUE), ]

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

# Find gene names
order_taqman$Gene <- gsub("taqman", "", order_taqman$Product, ignore.case = TRUE)
order_taqman$Gene <- gsub("human", "", order_taqman$Gene, ignore.case = TRUE)
order_taqman$Gene <- gsub("mouse", "", order_taqman$Gene, ignore.case = TRUE)
order_taqman$Gene <- gsub("rat", "", order_taqman$Gene, ignore.case = TRUE)

# fix other small things
order_taqman$in.pack <- gsub("[^0-9.-]", "", order_taqman$in.pack) %>%
  as.numeric()

# rename column names
colnames(order_taqman) <- c("order", "Cat_number", "Reactions", "Name", "received.date", "Gene_symbol")

#Gets names of excel files for inventory
taqman_folder <- "../../000_IntFys_documents/Reagents/Primers_TaqMan/"

inventory_file_human <- read.xlsx(paste0(taqman_folder, "1_Human_Taqman_assays.xlsx"), sheet = 1)
inventory_file_human$Species <- "Human"

inventory_file_mouse <- read.xlsx(paste0(taqman_folder, "2_Mouse_Taqman_assays.xlsx"), sheet = 1)
inventory_file_mouse$Species <- "Mouse"

inventory_file_rat <- read.xlsx(paste0(taqman_folder, "3_Rat_Taqman_assays.xlsx"), sheet = 1)
inventory_file_rat$Species <- "Rat"

inventory_file <- rbind(
  inventory_file_human,
  inventory_file_mouse,
  inventory_file_rat
)

# common boxes don't need names
inventory_file$Name <- NULL

# merge protein and gene names
inventory_file$Gene <- paste(inventory_file$Gene, inventory_file$Protein)
inventory_file$Gene <- gsub("NA", "", inventory_file$Gene)
inventory_file$Protein <- NULL

# full inventory - left join because everything in order list should be in the inventory
full_inventory <- full_join(order_taqman,
                            inventory_file,
                            by = "Cat_number")
colnames(full_inventory)

# merge the "gene columns"
full_inventory$Gene_symbol.x[is.na(full_inventory$Gene_symbol.x)] <- ""
full_inventory$Gene_symbol.y[is.na(full_inventory$Gene_symbol.y)] <- ""
full_inventory$Gene <- paste(full_inventory$Gene_symbol.x, full_inventory$Gene_symbol.y)
full_inventory$Gene_symbol.x <- NULL
full_inventory$Gene_symbol.y <- NULL

# order columns
full_inventory <- full_inventory[,
                                 c("Gene", "Alternate_name", "Cat_number", "Species", 
                                   "Name", "Box", "Pos", "Exp_date", "received.date")]

# clean rows without cat.number
full_inventory <- full_inventory[!is.na(full_inventory$Cat_number),]

# clean Species
table(full_inventory$Species, useNA = "ifany")

species_from_ref <- full_inventory$Cat_number
species_from_ref <- gsub("Hs.*", "Human", species_from_ref, ignore.case = T)
species_from_ref <- gsub("Mm.*", "Mouse", species_from_ref, ignore.case = T)
species_from_ref <- gsub("Rn.*", "Rat", species_from_ref, ignore.case = T)

# replace NA by the new species
full_inventory$Species <- ifelse(is.na(full_inventory$Species),
                                 species_from_ref,
                                 full_inventory$Species
                                 )
table(full_inventory$Species, useNA = "ifany")

# Remove names of people who left the lab OR replace with name of person who took over the probes
# Probes from those people should have been transfered into the common boxes
full_inventory$Name <- gsub(".*for all.*", "Common", full_inventory$Name)
full_inventory$Name <- gsub(" .*", "", full_inventory$Name)
full_inventory$Name <- gsub("/.*", "", full_inventory$Name)
table(full_inventory$Name)

old_members <- c("Ahmed",
                 "Elena",
                 "Flavia",
                 "Ilke",
                 "Jon",
                 "Logan",
                 "Lucile",
                 "Lukas",
                 "Rosamaria")
full_inventory <- full_inventory[!full_inventory$Name %in% old_members,]
table(full_inventory$Name)

#--------------------------------------------------------------------------------
#some more cleaning of gene names
#--------------------------------------------------------------------------------
full_inventory$Gene <- gsub("assay ", "", full_inventory$Gene, ignore.case = TRUE)
full_inventory$Gene <- gsub("Primer ", "", full_inventory$Gene, ignore.case = TRUE)
full_inventory$Gene <- gsub("miRNA:", "", full_inventory$Gene, ignore.case = TRUE)
full_inventory$Gene <- gsub(".*miRNA Assay, ", "", full_inventory$Gene, ignore.case = TRUE)
full_inventory$Gene <- gsub("\\(FAM-MGB\\)", "", full_inventory$Gene, ignore.case = TRUE)

# Removing leading spaces
full_inventory$Gene <- sub("^\\s+", "", full_inventory$Gene)

# remove duplicates
full_inventory$Gene <- sapply(full_inventory$Gene, function(s) {
  words <- strsplit(s, " ")[[1]]  # Split the string into words
  unique_words <- unique(words)   # Remove duplicate words
  paste(unique_words, collapse = " ")  # Recombine into a single string
})

# save file for Shiny
saveRDS(full_inventory, "full_inventory.Rds")

#save date of last update
last_update <- format(Sys.Date(), "%B %d, %Y")
saveRDS(last_update, "last_update.Rds")


###########################################################################
#
# Best way to organize boxes
#
##########################################################################
# Function to find the optimal splits
find_splits <- function(cumsums, target, num_boxes) {
  split_indices <- numeric(num_boxes - 1)
  for (i in 1:(num_boxes - 1)) {
    split_indices[i] <- which(cumsums >= i * target)[1]
  }
  
  boxes <- vector("list", num_boxes)
  start_index <- 1
  for (j in 1:(num_boxes - 1)) {
    boxes[[j]] <- names(cumsums)[start_index:(split_indices[j])]
    start_index <- split_indices[j] + 1
  }
  boxes[[num_boxes]] <- names(cumsums)[start_index:length(cumsums)]
  return(boxes)
}

#------------------------------------------------------------------------
# Human
#------------------------------------------------------------------------
# extract the first letter of each mouse gene
human_genes <- full_inventory[full_inventory$Species == "Human",]
human_genes$Gene <- gsub(" ", "", human_genes$Gene)
human_genes$Gene <- toupper(human_genes$Gene)
human_letters <- substr(human_genes$Gene, 1, 1)
human_letters_counts <- table(human_letters)
barplot(human_letters_counts)

# Sort the table by the names (letters)
sorted_counts <- human_letters_counts[order(names(human_letters_counts))]

# Calculate cumulative sums
cumulative_sums <- cumsum(sorted_counts)

# Total sum of genes
total_sum <- sum(sorted_counts)

# Target sum for each box
target_per_box <- total_sum / 15

# Calculate the splits
boxes <- find_splits(cumulative_sums, target_per_box, 15)

# Print the optimal box distributions
print(boxes)

# Count the number of genes in each box
box_counts <- sapply(boxes, function(letters) sum(sorted_counts[names(sorted_counts) %in% letters]))

# Print summary of tubes per box
barplot(box_counts)



#------------------------------------------------------------------------
# Mouse (3 boxes)
#------------------------------------------------------------------------
# extract the first letter of each mouse gene
mouse_genes <- full_inventory[full_inventory$Species == "Mouse",]
mouse_genes$Gene <- gsub(" ", "", mouse_genes$Gene)
mouse_genes$Gene <- toupper(mouse_genes$Gene)
mouse_letters <- substr(mouse_genes$Gene, 1, 1)
mouse_letters_counts <- table(mouse_letters)
barplot(mouse_letters_counts)

# Sort the table by the names (letters)
sorted_counts <- mouse_letters_counts[order(names(mouse_letters_counts))]

# Calculate cumulative sums
cumulative_sums <- cumsum(sorted_counts)

# Total sum of genes
total_sum <- sum(sorted_counts)

# Target sum for each box
target_per_box <- total_sum / 3

# Calculate the splits
boxes <- find_splits(cumulative_sums, target_per_box, 3)

# Print the optimal box distributions
print(boxes)

