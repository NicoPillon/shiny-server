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
# Function to calculate t-tests - Acute studies
#===========================================================================================
# Function to perform t-tests for a single GEO considering all relevant parameters
fct_human_ttests_acute <- function(data, geo) {
  # Filter data for the selected GEO
  geo_data <- data %>% filter(GEO == geo)
  
  # Get the unique combinations of all relevant parameters
  group_combinations <- unique(geo_data %>% 
                                 select(predicted_sex, predicted_ageGroup, predicted_fitness, predicted_BMIgroup, Tissue, Health_group, 
                                        Procotol_ConcentricEccentricMixed))
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Loop through each combination of the relevant parameters
  for (i in 1:nrow(group_combinations)) {
    sex <- group_combinations$predicted_sex[i]
    age <- group_combinations$predicted_ageGroup[i]
    fitness <- group_combinations$predicted_fitness[i]
    bmi <- group_combinations$predicted_BMIgroup[i]
    tissue <- group_combinations$Tissue[i]
    health <- group_combinations$Health_group[i]
    protocol <- group_combinations$Procotol_ConcentricEccentricMixed[i]
    
    # Filter data for the current combination of parameters
    group_data <- geo_data %>% 
      filter(predicted_sex == sex, predicted_ageGroup == age, predicted_fitness == fitness, 
             predicted_BMIgroup == bmi, Tissue == tissue, Health_group == health, 
             Procotol_ConcentricEccentricMixed == protocol)
    
    # Extract the PRE data
    pre_data <- group_data %>% filter(Acute_exercise_hours == "PRE") %>% pull(gene_data)
    
    # Identify unique groups other than "PRE"
    post_conditions <- group_data %>% 
      filter(Acute_exercise_hours != "PRE") %>% 
      pull(Acute_exercise_hours) %>% 
      unique()
    
    # Loop through each post condition and perform the t-test
    for (condition in post_conditions) {
      # Extract data for the current condition
      post_data <- group_data %>% filter(Acute_exercise_hours == condition) %>% pull(gene_data)
      
      # Check if there are enough observations in both PRE and POST condition groups
      if (length(pre_data) < 2 || all(is.na(pre_data)) || length(post_data) < 2 || all(is.na(post_data))) {
        # Store NA results if data is insufficient
        results_list[[paste(sex, age, fitness, bmi, tissue, health, protocol, condition, sep = "_")]] <- data.frame(
          GEO = geo,
          predicted_sex = sex,
          predicted_ageGroup = age,
          predicted_fitness = fitness,
          predicted_BMIgroup = bmi,
          Tissue = tissue,
          Health_group = health,
          Procotol_ConcentricEccentricMixed = protocol,
          Acute_exercise_hours = condition,
          logFC = NA,
          CI.L = NA,
          CI.R = NA,
          p_value = NA,
          POST_n = length(post_data)  # Include the sample size for the post condition
        )
      } else {
        # Perform the t-test comparing POST to PRE
        t_test_result <- t.test(post_data, pre_data)
        
        # Store the results in a list
        results_list[[paste(sex, age, fitness, bmi, tissue, health, protocol, condition, sep = "_")]] <- data.frame(
          GEO = geo,
          predicted_sex = sex,
          predicted_ageGroup = age,
          predicted_fitness = fitness,
          predicted_BMIgroup = bmi,
          Tissue = tissue,
          Health_group = health,
          Procotol_ConcentricEccentricMixed = protocol,
          Acute_exercise_hours = condition,
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
# Function to calculate Mean, Sd for all pairs - Inactivity studies (PRE / D00)
#===========================================================================================
# Function to perform t-tests for a single GEO considering all relevant parameters
fct_human_ttests_inactivity <- function(data, geo) {
  # Filter data for the selected GEO
  geo_data <- data %>% filter(GEO == geo)
  
  # Get the unique combinations of all relevant parameters
  group_combinations <- unique(geo_data %>% 
                                 select(predicted_sex, predicted_ageGroup, predicted_fitness, predicted_BMIgroup, Tissue, Health_group, 
                                        Protocol_group))
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Loop through each combination of the relevant parameters
  for (i in 1:nrow(group_combinations)) {
    sex <- group_combinations$predicted_sex[i]
    age <- group_combinations$predicted_ageGroup[i]
    fitness <- group_combinations$predicted_fitness[i]
    bmi <- group_combinations$predicted_BMIgroup[i]
    tissue <- group_combinations$Tissue[i]
    health <- group_combinations$Health_group[i]
    protocol <- group_combinations$Protocol_group[i]
    
    # Filter data for the current combination of parameters
    group_data <- geo_data %>% 
      filter(predicted_sex == sex, predicted_ageGroup == age, predicted_fitness == fitness, 
             predicted_BMIgroup == bmi, Tissue == tissue, Health_group == health, 
             Protocol_group == protocol)
    
    # Extract the PRE data
    pre_data <- group_data %>% filter(Inactivity_days == "PRE") %>% pull(gene_data)
    
    # Identify unique groups other than "PRE"
    post_conditions <- group_data %>% 
      filter(Inactivity_days != "PRE") %>% 
      pull(Inactivity_days) %>% 
      unique()
    
    # Loop through each post condition and perform the t-test
    for (condition in post_conditions) {
      # Extract data for the current condition
      post_data <- group_data %>% filter(Inactivity_days == condition) %>% pull(gene_data)
      
      # Check if there are enough observations in both PRE and POST condition groups
      if (length(pre_data) < 2 || all(is.na(pre_data)) || length(post_data) < 2 || all(is.na(post_data))) {
        # Store NA results if data is insufficient
        results_list[[paste(sex, age, fitness, bmi, tissue, health, protocol, condition, sep = "_")]] <- data.frame(
          GEO = geo,
          predicted_sex = sex,
          predicted_ageGroup = age,
          predicted_fitness = fitness,
          predicted_BMIgroup = bmi,
          Tissue = tissue,
          Health_group = health,
          Protocol_group = protocol,
          Inactivity_days = condition,
          logFC = NA,
          CI.L = NA,
          CI.R = NA,
          p_value = NA,
          POST_n = length(post_data)  # Include the sample size for the post condition
        )
      } else {
        # Perform the t-test comparing POST to PRE
        t_test_result <- t.test(post_data, pre_data)
        
        # Store the results in a list
        results_list[[paste(sex, age, fitness, bmi, tissue, health, protocol, condition, sep = "_")]] <- data.frame(
          GEO = geo,
          predicted_sex = sex,
          predicted_ageGroup = age,
          predicted_fitness = fitness,
          predicted_BMIgroup = bmi,
          Tissue = tissue,
          Health_group = health,
          Protocol_group = protocol,
          Inactivity_days = condition,
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
# Function to calculate Mean, Sd for all pairs - Training studies (PRE / W00)
#===========================================================================================
# Function to perform t-tests for a single GEO considering all relevant parameters
fct_human_ttests_training <- function(data, geo) {
  # Filter data for the selected GEO
  geo_data <- data %>% filter(GEO == geo)
  
  # Get the unique combinations of all relevant parameters
  group_combinations <- unique(geo_data %>% 
                                 select(predicted_sex, predicted_ageGroup, predicted_fitness, predicted_BMIgroup, Tissue, Health_group, 
                                        predicted_TimeFromLastExerciseBout))
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Loop through each combination of the relevant parameters
  for (i in 1:nrow(group_combinations)) {
    sex <- group_combinations$predicted_sex[i]
    age <- group_combinations$predicted_ageGroup[i]
    fitness <- group_combinations$predicted_fitness[i]
    bmi <- group_combinations$predicted_BMIgroup[i]
    tissue <- group_combinations$Tissue[i]
    health <- group_combinations$Health_group[i]
    biopsy_time <- group_combinations$predicted_TimeFromLastExerciseBout[i]
    
    # Filter data for the current combination of parameters
    group_data <- geo_data %>% 
      filter(predicted_sex == sex, 
             predicted_ageGroup == age, 
             predicted_fitness == fitness, 
             predicted_BMIgroup == bmi, 
             Tissue == tissue, 
             Health_group == health, 
             predicted_TimeFromLastExerciseBout == biopsy_time)
    
    # Extract the PRE data
    pre_data <- group_data %>% filter(Training_weeks == "PRE") %>% pull(gene_data)
    
    # Identify unique groups other than "PRE"
    post_conditions <- group_data %>% 
      filter(Training_weeks != "PRE") %>% 
      pull(Training_weeks) %>% 
      unique()
    
    # Loop through each post condition and perform the t-test
    for (condition in post_conditions) {
      # Extract data for the current condition
      post_data <- group_data %>% filter(Training_weeks == condition) %>% pull(gene_data)
      
      # Check if there are enough observations in both PRE and POST condition groups
      if (length(pre_data) < 2 || all(is.na(pre_data)) || length(post_data) < 2 || all(is.na(post_data))) {
        # Store NA results if data is insufficient
        results_list[[paste(sex, age, fitness, bmi, tissue, health, biopsy_time, condition, sep = "_")]] <- data.frame(
          GEO = geo,
          predicted_sex = sex,
          predicted_ageGroup = age,
          predicted_fitness = fitness,
          predicted_BMIgroup = bmi,
          Tissue = tissue,
          Health_group = health,
          Training_weeks = condition,
          predicted_TimeFromLastExerciseBout = biopsy_time,
          logFC = NA,
          CI.L = NA,
          CI.R = NA,
          p_value = NA,
          POST_n = length(post_data)  # Include the sample size for the post condition
        )
      } else {
        # Perform the t-test comparing POST to PRE
        t_test_result <- t.test(post_data, pre_data)
        
        # Store the results in a list
        results_list[[paste(sex, age, fitness, bmi, tissue, health, biopsy_time, condition, sep = "_")]] <- data.frame(
          GEO = geo,
          predicted_sex = sex,
          predicted_ageGroup = age,
          predicted_fitness = fitness,
          predicted_BMIgroup = bmi,
          Tissue = tissue,
          Health_group = health,
          Training_weeks = condition,
          predicted_TimeFromLastExerciseBout = biopsy_time,
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
fct_human_performMetaAnalysis <- function(paired_data) {
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




#===========================================================================================
# Function to generate forest plots - HUMAN ACUTE
#===========================================================================================
fct_human_generateForestPlot_acute <- function(meta, fdr, paired_data, color) {

  tabledata <- data.frame(mean = c(NA, meta$data$logFC, NA, meta$beta), 
                          lower = c(NA, meta$data$CI.L, NA, meta$ci.lb),
                          upper = c(NA, meta$data$CI.R, NA, meta$ci.ub))
  
  tabletext <- cbind(c('Study', meta$data$GEO, NA, ""),
                     c('Protocol', as.character(meta$data$Procotol_ConcentricEccentricMixed), NA, ""),
                     c('Sex', as.character(meta$data$predicted_sex), NA, ""),
                     c('Age', as.character(meta$data$predicted_ageGroup), NA, ""),
                     c('BMI', as.character(meta$data$predicted_BMIgroup), NA, ""),
                     c('Muscle', as.character(meta$data$Tissue), NA, ""),
                     c('Health', as.character(meta$data$Health_group), NA, ""),
                     c('Fitness', as.character(meta$data$predicted_fitness), NA, ""),
                     c('Timepoint', as.character(meta$data$Acute_exercise_hours), NA, ""),
                     c("logFC", format(round(meta$data$logFC, digits = 2)), NA, format(round(meta$beta, digits = 2))),
                     c("P.Val", format(meta$data$p_value, scientific = TRUE, digits = 2), NA, format(fdr, scientific = TRUE, digits = 2)),
                     c("n", meta$data$POST_n, NA, sum(meta$data$POST_n)))
  
  # fix age
  tabletext[2:nrow(tabletext), 3] <- gsub("Age", "", tabletext[2:nrow(tabletext), 3])
  tabletext[2:nrow(tabletext), 3] <- gsub("\\.", "-", tabletext[2:nrow(tabletext), 3])
  
  forestplot(tabletext,
             tabledata, 
             new_page = TRUE,
             align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "r", "r", "r"),
             colgap = unit(3, "mm"),
             is.summary = c(TRUE, rep(FALSE, (nrow(meta$data))), TRUE),
             xlog = FALSE, txt_gp = theme_forestPlot, xlab = "logFC",
             col = fpColors(box = color, line = "black", summary = color),
             boxsize = c(NA, meta$data$POST_n / 20, NA, 1.5)
  )
}



#===========================================================================================
# Function to generate forest plots - HUMAN TRAINING
#===========================================================================================
fct_human_generateForestPlot_training <- function(meta, fdr, paired_data, color) {
  
  tabledata <- data.frame(mean = c(NA, meta$data$logFC, NA, meta$beta), 
                          lower = c(NA, meta$data$CI.L, NA, meta$ci.lb),
                          upper = c(NA, meta$data$CI.R, NA, meta$ci.ub))
  
  tabletext <- cbind(c('Study', meta$data$GEO, NA, ""),
                     c('Sex', as.character(meta$data$predicted_sex), NA, ""),
                     c('Age', as.character(meta$data$predicted_ageGroup), NA, ""),
                     c('BMI', as.character(meta$data$predicted_BMIgroup), NA, ""),
                     c('Muscle', as.character(meta$data$Tissue), NA, ""),
                     c('Health', as.character(meta$data$Health_group), NA, ""),
                     c('Fitness', as.character(meta$data$predicted_fitness), NA, ""),
                     c('Duration', as.character(meta$data$Training_weeks), NA, ""),
                     c("logFC", format(round(meta$data$logFC, digits = 2)), NA, format(round(meta$beta, digits = 2))),
                     c("P.Val", format(meta$data$p_value, scientific = TRUE, digits = 2), NA, format(fdr, scientific = TRUE, digits = 2)),
                     c("n", meta$data$POST_n, NA, sum(meta$data$POST_n)))
  
  # fix age
  tabletext[2:nrow(tabletext), 3] <- gsub("Age", "", tabletext[2:nrow(tabletext), 3])
  tabletext[2:nrow(tabletext), 3] <- gsub("\\.", "-", tabletext[2:nrow(tabletext), 3])
  
  forestplot(tabletext,
             tabledata, 
             new_page = TRUE,
             align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "r", "r", "r"),
             colgap = unit(3, "mm"),
             is.summary = c(TRUE, rep(FALSE, (nrow(meta$data))), TRUE),
             xlog = FALSE, txt_gp = theme_forestPlot, xlab = "logFC",
             col = fpColors(box = color, line = "black", summary = color),
             boxsize = c(NA, meta$data$POST_n / 20, NA, 1.5)
  )
}

