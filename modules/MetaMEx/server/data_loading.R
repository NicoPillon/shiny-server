#===========================================================================================
# Lists of the different studies, categories and gene names
#===========================================================================================
# Human
human_meta <- readRDS("annotation/human_metadata.Rds")
human_genes <- readRDS("annotation/human_genes.Rds")
human_input <- readRDS("annotation/human_input_categories.Rds")
human_references <- readRDS("annotation/human_references.Rds") # Load the table describing the studies

# Mouse
mouse_meta <- readRDS("annotation/mouse_metadata.Rds")
mouse_genes <- readRDS("annotation/mouse_genes.Rds")
mouse_input <- readRDS("annotation/mouse_input_categories.Rds")
mouse_references <- readRDS("annotation/mouse_references.Rds") # Load the table describing the studies

#===========================================================================================
# Function to transform human gene into mouse genename
#===========================================================================================
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#===========================================================================================
# URLs for sharing buttons
#===========================================================================================
url_twitter  <- "https://twitter.com/intent/tweet?text=MetaMEx:%20Meta-Analysis%20of%20skeletal%20muscle%20response%20to%20inactivity%20and%20exercise.%20@NicoPillon%20@JuleenRZierath&url=http://www.metamex.eu"
url_linkedin <- "https://www.linkedin.com/shareArticle?mini=true&url=http://www.metamex.eu&title=MetaMEx:%20Meta-Analysis%20of%20skeletal%20muscle%20response%20to%20inactivity%20and%20exercise.%20&summary=MetaMEx&source=LinkedIn"
url_facebook <- "https://www.facebook.com/sharer.php?u=https://www.metamex.eu&title=MetaMEx:%20Meta-Analysis%20of%20skeletal%20muscle%20response%20to%20inactivity%20and%20exercise.%20"

