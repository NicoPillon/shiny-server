#===========================================================================================
# Graphical parameters
#===========================================================================================
# forrest plot
theme_forestPlot <- fpTxtGp()
theme_forestPlot$ticks$cex <- 1.1 #tick labels
theme_forestPlot$xlab$cex <- 1.1
theme_forestPlot$label$cex <- 1.1
theme_forestPlot$summary$cex <- 1.2

# ggplots
theme_ggplot <- theme(plot.title  = element_text(face="bold", color="black", size=14, angle=0),
                      axis.text.x = element_text(color="black", size=12, angle=0, hjust = 0.5),
                      axis.text.y = element_text(color="black", size=12, angle=0, hjust = 1),
                      axis.title  = element_text(face="bold", color="black", size=14, angle=0),
                      legend.text   = element_text(color="black", size=14, angle=0),
                      legend.title  = element_blank(),
                      legend.position="none")


#===========================================================================================
# Generalized function to perform t-tests for a single GEO considering all relevant parameters
#===========================================================================================
fct_mouse_ttests <- function(data, geo, condition_col) {
  # Filter data for the selected GEO
  geo_data <- data %>% filter(GEO == geo)
  
  # Get the unique combinations of all relevant parameters
  group_combinations <- unique(geo_data %>% 
                                 select(predicted_sex, Age, Genotype, Tissue))
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Loop through each combination of the relevant parameters
  for (i in 1:nrow(group_combinations)) {
    sex <- group_combinations$predicted_sex[i]
    age <- group_combinations$Age[i]
    genotype <- group_combinations$Genotype[i]
    tissue <- group_combinations$Tissue[i]
    
    # Filter data for the current combination of parameters
    group_data <- geo_data %>% 
      filter(predicted_sex == sex, 
             Age == age, 
             Genotype == genotype, 
             Tissue == tissue)
    
    # Extract the control group data
    pre_data <- group_data %>% filter(.data[[condition_col]] == "CTRL") %>% pull(gene_data)
    
    # Identify unique groups other than the control group
    post_conditions <- group_data %>% 
      filter(.data[[condition_col]] != "CTRL") %>% 
      pull(.data[[condition_col]]) %>% 
      unique()
    
    # Loop through each post condition and perform the t-test
    for (condition in post_conditions) {
      # Extract data for the current condition
      post_data <- group_data %>% filter(.data[[condition_col]] == condition) %>% pull(gene_data)
      
      # Check if there are enough observations in both control and post condition groups
      if (length(pre_data) < 2 || all(is.na(pre_data)) || length(post_data) < 2 || all(is.na(post_data))) {
        # Store NA results if data is insufficient
        results_list[[paste(sex, age, genotype, tissue, condition, sep = "_")]] <- data.frame(
          GEO = geo,
          predicted_sex = sex,
          Age = age,
          Genotype = genotype,
          Tissue = tissue,
          Condition = condition,
          logFC = NA,
          CI.L = NA,
          CI.R = NA,
          p_value = NA,
          POST_n = length(post_data)  # Include the sample size for the post condition
        )
      } else {
        # Perform the t-test comparing POST to CTRL
        t_test_result <- t.test(post_data, pre_data)
        
        # Store the results in a list
        results_list[[paste(sex, age, genotype, tissue, condition, sep = "_")]] <- data.frame(
          GEO = geo,
          predicted_sex = sex,
          Age = age,
          Genotype = genotype,
          Tissue = tissue,
          Condition = condition,
          logFC = mean(post_data, na.rm = TRUE) - mean(pre_data, na.rm = TRUE),
          CI.L = t_test_result$conf.int[1],
          CI.R = t_test_result$conf.int[2],
          p_value = t_test_result$p.value,
          POST_n = length(post_data)  # Include the sample size for the post condition
        )
      }
    }
  }
  
  # Combine the list of results into a single data frame
  if (length(results_list) > 0) {
    results_df <- do.call(rbind, results_list) %>%
      data.frame(row.names = NULL) %>% 
      arrange(desc(logFC))
  } else {
    results_df <- data.frame()  # Return an empty data frame if no results
  }
  
  # remove data with n = 1
  results_df <- results_df[results_df$POST_n > 1,]
  
  return(results_df)
}



#===========================================================================================
# Function to perform meta-analysis from logFC and CI
#===========================================================================================
fct_mouse_performMetaAnalysis <- function(paired_data) {
  paired_data <- paired_data
  meta <- rma(yi=paired_data$logFC,
              sei=(paired_data$CI.R - paired_data$CI.L) / 3.92,
              method = "REML", measure = "MD", 
              data = paired_data, 
              weighted = TRUE, 
              weights = paired_data$POST_n)
  
  fdr <- p.adjust(meta$pval, method = 'bonferroni', n = nrow(paired_data))
  
  return(list(meta = meta, fdr = fdr, paired_data = paired_data))
}

