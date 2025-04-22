# Load necessary libraries
library(biomaRt)
library(rstudioapi)
library(openxlsx)
library(dplyr)
library(readxl)

# Set working directory to source file location
setwd(dirname(getActiveDocumentContext()$path))

###################################################################################
# How to organize boxes?
###################################################################################
# Connect to the Ensembl database and retrieve human gene information
# use only protein coding genes
human = useMart("ensembl", dataset="hsapiens_gene_ensembl")
human_genes = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "biotype",
                    values = "protein_coding",
                    mart = human,
                    uniqueRows = TRUE)

# Remove rows with missing gene names and convert names to uppercase
human_genes <- na.omit(human_genes)
human_genes <- human_genes[!human_genes$external_gene_name %in% c(""),]
human_genes$external_gene_name <- toupper(human_genes$external_gene_name)

# Extract the first letter of each gene name and remove genes starting with numbers
human_letters <- substr(human_genes$external_gene_name, 1, 1)
human_letters <- human_letters[!human_letters %in% as.character(1:9)]

# Count the frequency of each letter
human_letters_counts <- table(human_letters)
barplot(human_letters_counts)

# Sort the table by letter names
sorted_counts <- human_letters_counts[order(names(human_letters_counts))]

# Calculate cumulative sums and total sum of genes
cumulative_sums <- cumsum(sorted_counts)
total_sum <- sum(sorted_counts)

# Function to calculate optimal number of boxes within a given range
calculate_optimal_boxes <- function(total_sum, sorted_counts, min_boxes, max_boxes) {
  num_letters <- length(sorted_counts)
  optimal_boxes <- min_boxes
  min_max_difference <- Inf
  
  for (num_boxes in min_boxes:max_boxes) {
    target_per_box <- total_sum / num_boxes
    box_sums <- rep(0, num_boxes)
    current_box <- 1
    
    for (i in seq_along(sorted_counts)) {
      if (box_sums[current_box] + sorted_counts[i] > target_per_box && current_box < num_boxes) {
        current_box <- current_box + 1
      }
      box_sums[current_box] <- box_sums[current_box] + sorted_counts[i]
    }
    
    max_difference <- max(abs(box_sums - target_per_box))
    
    if (max_difference < min_max_difference) {
      min_max_difference <- max_difference
      optimal_boxes <- num_boxes
    }
  }
  
  return(optimal_boxes)
}

# Calculate optimal number of boxes within the range of 6 to 10
optimal_boxes <- calculate_optimal_boxes(total_sum, sorted_counts, 6, 10)
optimal_boxes

# Function to assign letters to boxes in alphabetical order
assign_boxes_in_order <- function(counts, num_boxes) {
  letters <- names(counts)
  num_letters <- length(letters)
  box_size <- ceiling(num_letters / num_boxes)
  
  boxes <- vector("list", num_boxes)
  start_index <- 1
  
  for (i in 1:num_boxes) {
    end_index <- min(start_index + box_size - 1, num_letters)
    boxes[[i]] <- letters[start_index:end_index]
    start_index <- end_index + 1
    if (start_index > num_letters) break
  }
  
  return(boxes)
}

# Assign letters to boxes using optimal number of boxes
boxes <- assign_boxes_in_order(sorted_counts, optimal_boxes)

# Calculate the number of genes in each box and plot
genes_per_box <- sapply(boxes, function(letters) sum(sorted_counts[letters]))
names(genes_per_box) <- paste("Box", 1:length(genes_per_box))
barplot(genes_per_box)

# Print the optimal number of boxes, the boxes, and number of genes per box
print(paste("Optimal number of boxes:", optimal_boxes))
print(boxes)
print(genes_per_box)
