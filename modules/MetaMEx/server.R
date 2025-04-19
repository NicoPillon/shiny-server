
####################################################################################################################
# Define server logic ##############################################################################################
server <- function(input, output, session) {
  
  
  #=======================================================================================
  #
  # Load data
  # small annotation datasets are loading in the function files
  #=======================================================================================
  withProgress(message = 'Loading data', value = 1, max=13, {
    
    #--------------------------------------------------------------------------------------------------------------
    # Data matrix for human
    #--------------------------------------------------------------------------------------------------------------
    incProgress(1, detail="Human database")

    # Initialize an empty list to store each batch
    matrix_list <- list()
    
    # Loop through all the batch files
    i <- 1
    repeat {
      file_path <- paste0("data/human_matrix_", i, ".feather")
      
      # Check if the file exists
      if (!file.exists(file_path)) break
      
      # Read the batch and add it to the list
      batch_matrix <- read_feather(file_path)
      matrix_list[[i]] <- batch_matrix
      
      # Increment progress bar
      #incProgress(1 / 15)
      
      i <- i + 1
    }
    
    # Merge all batches back into a single matrix
    human_matrix <- do.call(rbind, matrix_list) %>%
      data.frame(row.names = human_genes$SYMBOL)
    
    #--------------------------------------------------------------------------------------------------------------
    # Data matrix for mouse
    #--------------------------------------------------------------------------------------------------------------
    incProgress(1, detail="Mouse database")
    
    # Initialize an empty list to store each batch
    matrix_list <- list()
    
    # Loop through all the batch files
    i <- 1
    repeat {
      file_path <- paste0("data/mouse_matrix_", i, ".feather")
      
      # Check if the file exists
      if (!file.exists(file_path)) break
      
      # Read the batch and add it to the list
      batch_matrix <- read_feather(file_path)
      matrix_list[[i]] <- batch_matrix
      
      # Increment progress bar
      #incProgress(1 / 15)
      
      i <- i + 1
    }
    
    # Merge all batches back into a single matrix
    mouse_matrix <- do.call(rbind, matrix_list) %>%
      data.frame(row.names = mouse_genes$SYMBOL)
    
    
    #--------------------------------------------------------------------------------------------------------------
    # Data for timeline
    #--------------------------------------------------------------------------------------------------------------
    incProgress(1, detail="Timeline analysis")
    human_timelineStats_acute <- data.frame(read_feather("data/human_timelineStats_acute.feather"), row.names=1)
    human_timelineStats_inactivity <- data.frame(read_feather("data/human_timelineStats_inactivity.feather"), row.names=1)
    
    
  })
  
  
  #=======================================================================================
  #
  # Set up reactivity of home page
  #
  #=======================================================================================
  updateSelectizeInput(session, 'genename_home',
                       choices=human_genes$SYMBOL,
                       server=TRUE,
                       selected='NR4A3',
                       options=NULL)

  # Synchronize genenames according to home page
  observeEvent(input$genename_home, { updateSelectizeInput(session, 'genename_metaanalysis_human',
                                                           choices=human_genes$SYMBOL,
                                                           server=TRUE,
                                                           selected=input$genename_home, options=NULL) })
  observeEvent(input$genename_home, { updateSelectizeInput(session, 'genename_metaanalysis_mouse',
                                                           choices=mouse_genes$SYMBOL,
                                                           server=TRUE,
                                                           selected=firstup(input$genename_home), options=NULL) })
  observeEvent(input$genename_metaanalysis_human, { updateSelectizeInput(session, 'genename_metaanalysis_mouse',
                                                           choices=mouse_genes$SYMBOL,
                                                           server=TRUE,
                                                           selected=firstup(input$genename_metaanalysis_human), options=NULL) })

  # Observe the clicks on buttons on home page and transfer to tab
  observeEvent(input$jumpToHumanOverview, {updateTabsetPanel(session, "inTabset", selected="panelAppHuman")
    updateTabsetPanel(session, "inTabsetMeta", selected="panelAppOverviewHuman")})

  observeEvent(input$jumpToMouseOverview, {updateTabsetPanel(session, "inTabset", selected="panelAppMouse")
    updateTabsetPanel(session, "inTabsetMeta", selected="panelAppOverviewMouse")})

  observeEvent(input$jumpToHelp, {updateTabsetPanel(session, "inTabset", selected="Tutorial")
    updateTabsetPanel(session, "inTabsetMeta", selected="Tutorial")})

  
  #=======================================================================================
  #Make all checkboxes selected by default - necessary for the select all(none) button to work
  #=======================================================================================
  observe({ updateCheckboxGroupInput(session, 'human_muscle', 
                                     choices = human_input[['muscle_choice']], 
                                     selected = if (input$human_muscle_allnone) human_input[['muscle_choice']])})
  observe({ updateCheckboxGroupInput(session, 'human_disease',             
                                     choices = human_input[['disease_choice']], 
                                     selected = if (input$human_disease_allnone) human_input[['disease_choice']])})
  observe({ updateCheckboxGroupInput(session, 'human_exercise_type',        
                                     choices = human_input[['exercise_type_choice']], 
                                     selected = if (input$human_exercise_type_allnone) human_input[['exercise_type_choice']])})
  observe({ updateCheckboxGroupInput(session, 'human_inactivity_protocol', 
                                     choices=human_input[['inactivity_protocol_choice']],
                                     selected = if (input$human_inactivity_protocol_allnone) human_input[['inactivity_protocol_choice']])})
  
  
  observe({ updateCheckboxGroupInput(session, 'mouse_muscle',
                                     choices = mouse_input[['muscle_choice']],
                                     selected = if (input$mouse_muscle_allnone) mouse_input[['muscle_choice']])})
  observe({ updateCheckboxGroupInput(session, 'mouse_disease',
                                     choices = mouse_input[['disease_choice']],
                                     selected = if (input$mouse_disease_allnone) mouse_input[['disease_choice']])})
  observe({ updateCheckboxGroupInput(session, 'mouse_exercise_type',
                                     choices = mouse_input[['exercise_type_choice']],
                                     selected = if (input$mouse_exercise_type_allnone) mouse_input[['exercise_type_choice']])})
  observe({ updateCheckboxGroupInput(session, 'mouse_inactivity_protocol',
                                     choices=mouse_input[['inactivity_protocol_choice']],
                                     selected = if (input$mouse_inactivity_protocol_allnone) mouse_input[['inactivity_protocol_choice']])})

  
  #=======================================================================================
  #
  # Make annotation tables for legend 
  #
  #=======================================================================================
  output$reftable_human <- DT::renderDataTable(escape = FALSE, rownames = FALSE, { 
    human_references 
    
    # Modify the column name
    colnames(human_references)[2] <- "Male / Female"
    
    # Render the datatable with column alignment
    DT::datatable(
      human_references, 
      escape = FALSE, 
      rownames = FALSE,
      options = list(
        columnDefs = list(
          list(className = 'dt-center', targets = 1),  # Align second column to the center
          list(className = 'dt-center', targets = 2),  # Align third column to the right
          list(className = 'dt-center', targets = 3)  # Align third column to the right
          # Add more lines for additional columns if needed
        )
      )
    )
    } )
  
  output$reftable_mouse <- DT::renderDataTable(escape = FALSE, rownames = FALSE, { 
    mouse_references 
    
    # Modify the column name
    colnames(mouse_references)[2] <- "Male / Female"
    
    # Render the datatable with column alignment
    DT::datatable(
      mouse_references, 
      escape = FALSE, 
      rownames = FALSE,
      options = list(
        columnDefs = list(
          list(className = 'dt-center', targets = 1),  # Align second column to the center
          list(className = 'dt-center', targets = 2)  # Align third column to the right
          # Add more lines for additional columns if needed
        )
      )
    )
    } )
  
  
  #=======================================================================================
  #
  # Forest plot - Human Acute Aerobic
  #
  #=======================================================================================
  human_metaAnalysisData_AA <- reactive({
    tryCatch({  
      # select acute aerobic studies
      selected_meta <- human_meta[!is.na(human_meta$Acute_exercise_hours),]
      selected_meta <- selected_meta[!grepl("W", selected_meta$Training_weeks),]
      selected_meta <- selected_meta[selected_meta$Protocol_group %in% c("Basal", "Aerobic"),]
      
      # select corresponding matrix
      selected_matrix <- human_matrix[,selected_meta$SampleID]
      
      # selected gene
      selected_gene <- "NR4A3"
      selected_gene <- input$genename_metaanalysis_human
      
      # merge data
      selected_data <- data.frame(selected_meta[,c("GEO", 
                                                   "predicted_sex", 
                                                   "predicted_ageGroup", 
                                                   "predicted_fitness", 
                                                   "predicted_BMIgroup", 
                                                   "Tissue", 
                                                   "Health_group",
                                                   "Acute_exercise_hours", 
                                                   "Procotol_ConcentricEccentricMixed")],
                                  gene_data = as.numeric(selected_matrix[selected_gene,]))
      
      # Apply the function to each GEO in the dataset
      all_geos <- unique(selected_data$GEO)
      paired_data <- do.call(rbind, lapply(all_geos, function(geo) fct_human_ttests_acute(selected_data, geo)))
      all_geos[!all_geos %in% unique(paired_data$GEO)]
      
      # order data
      paired_data <- paired_data[order(paired_data$logFC, decreasing = TRUE), ]
      
      # change acute time for numeric
      paired_data$Acute_exercise_hours <- as.numeric(gsub("H", "", paired_data$Acute_exercise_hours))
      
      # Filter population of interest
      paired_data <- dplyr::filter(paired_data,
                                   GEO %in% input$human_studies_AA,
                                   predicted_sex %in% input$human_sex,
                                   predicted_ageGroup %in% input$human_age,
                                   predicted_fitness %in% input$human_fitness,
                                   predicted_BMIgroup %in% input$human_weight,
                                   Tissue %in% input$human_muscle,
                                   Health_group %in% input$human_disease,
                                   Acute_exercise_hours >= input$human_acute_biopsy[1] & Acute_exercise_hours <= input$human_acute_biopsy[2],
                                   Procotol_ConcentricEccentricMixed %in% input$human_exercise_type)
      
      # Calculate meta-analysis
      meta_analysis_results <- fct_human_performMetaAnalysis(paired_data)
      
    }, error=function(e) NULL)
  }) 
  
  output$plot_human_AA <- renderPlot({ 
    
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculating", value = 1)
    
    # call data
    meta_analysis_results <- human_metaAnalysisData_AA()
    
    # Handle empty data case: this validates that paired_data is not null and has rows.
    shiny::validate(
      need(!is.null(meta_analysis_results$paired_data) && nrow(meta_analysis_results$paired_data) > 0, 
           "No studies found - try different selection criteria")
    )
    
    # Generate forest plot
    fct_human_generateForestPlot_acute(meta = meta_analysis_results$meta, 
                                   fdr = meta_analysis_results$fdr, 
                                   paired_data = paired_data,
                                   color = "#D55E00")
  })
  
  output$human_forestPlotDynamic_AA <- renderUI({ 
    # call data
    meta_analysis_results <- human_metaAnalysisData_AA()
    
    # Determine height based on whether meta_analysis_results is NULL
    if (is.null(meta_analysis_results) || is.null(meta_analysis_results$paired_data)) {
      # If meta_analysis_results is NULL, set height to 50px
      height <- 50
    } else {
      # Calculate height based on the number of rows in paired_data
      height <- 150 + 17 * nrow(meta_analysis_results$paired_data)
    }
    
    plotOutput("plot_human_AA", height = paste0(height, "px"))
  })
  
  
  #=======================================================================================
  #
  # Forest plot - Human Acute HIT
  #
  #=======================================================================================
  human_metaAnalysisData_AH <- reactive({
    tryCatch({  
      # select acute aerobic studies
      selected_meta <- human_meta[!is.na(human_meta$Acute_exercise_hours),]
      selected_meta <- selected_meta[!grepl("W", selected_meta$Training_weeks),]
      selected_meta <- selected_meta[selected_meta$Protocol_group %in% c("Basal", "HIT"),]
      
      # select corresponding matrix
      selected_matrix <- human_matrix[,selected_meta$SampleID]
      
      # selected gene
      selected_gene <- "NR4A3"
      selected_gene <- input$genename_metaanalysis_human
      
      # merge data
      selected_data <- data.frame(selected_meta[,c("GEO", 
                                                   "predicted_sex", "predicted_ageGroup", "predicted_fitness", "predicted_BMIgroup", "Tissue", "Health_group",
                                                   "Acute_exercise_hours", "Procotol_ConcentricEccentricMixed")],
                                  gene_data = as.numeric(selected_matrix[selected_gene,]))
      
      # Apply the function to each GEO in the dataset
      all_geos <- unique(selected_data$GEO)
      paired_data <- do.call(rbind, lapply(all_geos, function(geo) fct_human_ttests_acute(selected_data, geo)))
      all_geos[!all_geos %in% unique(paired_data$GEO)]
      
      # order data
      paired_data <- paired_data[order(paired_data$logFC, decreasing = TRUE), ]
      
      # change acute time for numeric
      paired_data$Acute_exercise_hours <- as.numeric(gsub("H", "", paired_data$Acute_exercise_hours))
      
      # Filter population of interest
      paired_data <- dplyr::filter(paired_data,
                                   GEO %in% input$human_studies_AH,
                                   predicted_sex %in% input$human_sex,
                                   predicted_ageGroup %in% input$human_age,
                                   predicted_fitness %in% input$human_fitness,
                                   predicted_BMIgroup %in% input$human_weight,
                                   Tissue %in% input$human_muscle,
                                   Health_group %in% input$human_disease,
                                   Acute_exercise_hours >= input$human_acute_biopsy[1] & Acute_exercise_hours <= input$human_acute_biopsy[2],
                                   Procotol_ConcentricEccentricMixed %in% input$human_exercise_type)
      
      # Calculate meta-analysis
      meta_analysis_results <- fct_human_performMetaAnalysis(paired_data)
      
    }, error=function(e) NULL)
  }) 
  
  
  output$plot_human_AH <- renderPlot({ 
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculating", value = 1)
    
    # call data
    meta_analysis_results <- human_metaAnalysisData_AH()
    
    # Handle empty data case: this validates that paired_data is not null and has rows.
    shiny::validate(
      need(!is.null(meta_analysis_results$paired_data) && nrow(meta_analysis_results$paired_data) > 0, 
           "No studies found - try different selection criteria")
    )

    # Generate forest plot
    fct_human_generateForestPlot_acute(meta = meta_analysis_results$meta, 
                             fdr = meta_analysis_results$fdr, 
                             paired_data = meta_analysis_results$paired_data,
                             color = "#D55E00")
  })
  

  output$human_forestPlotDynamic_AH <- renderUI({ 
    # call data
    meta_analysis_results <- human_metaAnalysisData_AH()
    
    # Determine height based on whether meta_analysis_results is NULL
    if (is.null(meta_analysis_results) || is.null(meta_analysis_results$paired_data)) {
      # If meta_analysis_results is NULL, set height to 50px
      height <- 50
    } else {
      # Calculate height based on the number of rows in paired_data
      height <- 150 + 17 * nrow(meta_analysis_results$paired_data)
    }
    
    plotOutput("plot_human_AH", height = paste0(height, "px"))
  })
  
  
  
  #=======================================================================================
  #
  # Forest plot - Human Acute Resistance
  #
  #=======================================================================================
  human_metaAnalysisData_AR <- reactive({
    tryCatch({  
      # select studies
      selected_meta <- human_meta[!is.na(human_meta$Acute_exercise_hours),]
      selected_meta <- selected_meta[!grepl("W", selected_meta$Training_weeks),]
      selected_meta <- selected_meta[selected_meta$Protocol_group %in% c("Basal", "Resistance"),]
      
      # select corresponding matrix
      selected_matrix <- human_matrix[,selected_meta$SampleID]
      
      # selected gene
      selected_gene <- "NR4A3"
      selected_gene <- input$genename_metaanalysis_human
      
      # merge data
      selected_data <- data.frame(selected_meta[,c("GEO", 
                                                   "predicted_sex", "predicted_ageGroup", "predicted_fitness", "predicted_BMIgroup", "Tissue", "Health_group",
                                                   "Acute_exercise_hours", "Procotol_ConcentricEccentricMixed")],
                                  gene_data = as.numeric(selected_matrix[selected_gene,]))
      
      # Apply the function to each GEO in the dataset
      all_geos <- unique(selected_data$GEO)
      paired_data <- do.call(rbind, lapply(all_geos, function(geo) fct_human_ttests_acute(selected_data, geo)))
      paired_data <- paired_data[order(paired_data$logFC, decreasing = TRUE), ]
      all_geos[!all_geos %in% unique(paired_data$GEO)]
      
      # change acute time for numeric
      paired_data$Acute_exercise_hours <- as.numeric(gsub("H", "", paired_data$Acute_exercise_hours))
      
      # Filter population of interest
      paired_data <- dplyr::filter(paired_data,
                                   GEO %in% input$human_studies_AR,
                                   predicted_sex %in% input$human_sex,
                                   predicted_ageGroup %in% input$human_age,
                                   predicted_fitness %in% input$human_fitness,
                                   predicted_BMIgroup %in% input$human_weight,
                                   Tissue %in% input$human_muscle,
                                   Health_group %in% input$human_disease,
                                   Acute_exercise_hours >= input$human_acute_biopsy[1] & Acute_exercise_hours <= input$human_acute_biopsy[2],
                                   Procotol_ConcentricEccentricMixed %in% input$human_exercise_type)
      
      # Calculate meta-analysis
      meta_analysis_results <- fct_human_performMetaAnalysis(paired_data)
      
    }, error=function(e) NULL)
  }) 
  
  
  output$plot_human_AR <- renderPlot({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculating", value = 1)
    
    # call data
    meta_analysis_results <- human_metaAnalysisData_AR()
    
    # Handle empty data case: this validates that paired_data is not null and has rows.
    shiny::validate(
      need(!is.null(meta_analysis_results$paired_data) && nrow(meta_analysis_results$paired_data) > 0, 
           "No studies found - try different selection criteria")
    )

    # Generate forest plot
    fct_human_generateForestPlot_acute(meta = meta_analysis_results$meta, 
                             fdr = meta_analysis_results$fdr, 
                             paired_data = meta_analysis_results$paired_data,
                             color = "#D55E00")
  })
  
  
  output$human_forestPlotDynamic_AR <- renderUI({ 
    # call data
    meta_analysis_results <- human_metaAnalysisData_AR()
    
    # Determine height based on whether meta_analysis_results is NULL
    if (is.null(meta_analysis_results) || is.null(meta_analysis_results$paired_data)) {
      # If meta_analysis_results is NULL, set height to 50px
      height <- 50
    } else {
      # Calculate height based on the number of rows in paired_data
      height <- 150 + 17 * nrow(meta_analysis_results$paired_data)
    }
    
    plotOutput("plot_human_AR", height = paste0(height, "px"))
  })
  
  
  
  #=======================================================================================
  #
  # Forest plot - Human Inactivity
  #
  #=======================================================================================
  human_metaAnalysisData_IN <- reactive({
    tryCatch({  
      # select studies
      selected_meta <- human_meta[!is.na(human_meta$Inactivity_days),]
      selected_meta <- selected_meta[!grepl("RELOAD", selected_meta$Inactivity_days),]
      
      # select corresponding matrix
      selected_matrix <- human_matrix[,selected_meta$SampleID]
      
      # selected gene
      selected_gene <- "NR4A3"
      selected_gene <- input$genename_metaanalysis_human
      
      # merge data
      selected_data <- data.frame(selected_meta[,c("GEO", 
                                                   "predicted_sex", "predicted_ageGroup", 
                                                   "predicted_fitness", "predicted_BMIgroup", 
                                                   "Tissue", "Health_group",
                                                   "Inactivity_days", "Protocol_group")],
                                  gene_data = as.numeric(selected_matrix[selected_gene,]))
      
      # Apply the function to each GEO in the dataset
      all_geos <- unique(selected_data$GEO)
      paired_data <- do.call(rbind, lapply(all_geos, function(geo) fct_human_ttests_inactivity(selected_data, geo)))
      paired_data <- paired_data[order(paired_data$logFC, decreasing = TRUE), ]
      all_geos[!all_geos %in% unique(paired_data$GEO)]
      
      # change acute time for numeric
      paired_data$Inactivity_days <- as.numeric(gsub("D", "", paired_data$Inactivity_days))
      
      # Filter population of interest
      paired_data <- dplyr::filter(paired_data,
                                   GEO %in% input$human_studies_IN,
                                   predicted_sex %in% input$human_sex,
                                   predicted_ageGroup %in% input$human_age,
                                   predicted_fitness %in% input$human_fitness,
                                   predicted_BMIgroup %in% input$human_weight,
                                   Tissue %in% input$human_muscle,
                                   Health_group %in% input$human_disease,
                                   Protocol_group %in% input$human_inactivity_protocol,
                                   Inactivity_days >= input$human_inactivity_duration[1] & Inactivity_days <= input$human_inactivity_duration[2]
                                   )
      
      # Calculate meta-analysis
      meta_analysis_results <- fct_human_performMetaAnalysis(paired_data)
      
    }, error=function(e) NULL)
  }) 
  
  output$plot_human_IN <- renderPlot({ 
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculating", value = 1)
    
    # call data
    meta_analysis_results <- human_metaAnalysisData_IN()
    
    # Handle empty data case: this validates that paired_data is not null and has rows.
    shiny::validate(
      need(!is.null(meta_analysis_results$paired_data) && nrow(meta_analysis_results$paired_data) > 0, 
           "No studies found - try different selection criteria")
    )
    
    # extract data
    meta <- meta_analysis_results$meta
    fdr <- meta_analysis_results$fdr
    
    # make tables
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
                       c('Protocol', as.character(meta$data$Protocol_group), NA, ""),
                       c('Timepoint', as.character(meta$data$Inactivity_days), NA, ""),
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
               col = fpColors(box = "#CC79A7", line = "black", summary = "#CC79A7"),
               boxsize = c(NA, meta$data$POST_n / 20, NA, 1.5)
    )
  })
  
  
  output$human_forestPlotDynamic_IN <- renderUI({ 
    # call data
    meta_analysis_results <- human_metaAnalysisData_IN()
    
    # Determine height based on whether meta_analysis_results is NULL
    if (is.null(meta_analysis_results) || is.null(meta_analysis_results$paired_data)) {
      # If meta_analysis_results is NULL, set height to 50px
      height <- 50
    } else {
      # Calculate height based on the number of rows in paired_data
      height <- 150 + 17 * nrow(meta_analysis_results$paired_data)
    }
    
    plotOutput("plot_human_IN", height = paste0(height, "px"))
  })
  
  
  
  #=======================================================================================
  #
  # Forest plot - Human Aerobic Training
  #
  #=======================================================================================
  human_metaAnalysisData_TA <- reactive({
    tryCatch({  
      # select training aerobic studies
      selected_meta <- human_meta[!is.na(human_meta$Training_weeks),]
      selected_meta <- selected_meta[!grepl("H0", selected_meta$Acute_exercise_hours),]
      selected_meta <- selected_meta[selected_meta$Protocol_group %in% c("Basal", "Aerobic"),]
      
      # select corresponding matrix
      selected_matrix <- human_matrix[,selected_meta$SampleID]
      
      # selected gene
      selected_gene <- "NR4A3"
      selected_gene <- input$genename_metaanalysis_human
      
      # merge data
      selected_data <- data.frame(selected_meta[,c("GEO", 
                                                   "predicted_sex", "predicted_ageGroup", "predicted_fitness", "predicted_BMIgroup", "Tissue", "Health_group",
                                                   "Training_weeks", "predicted_TimeFromLastExerciseBout")],
                                  gene_data = as.numeric(selected_matrix[selected_gene,]))
      
      # Apply the function to each GEO in the dataset
      all_geos <- unique(selected_data$GEO)
      paired_data <- do.call(rbind, lapply(all_geos, function(geo) fct_human_ttests_training(selected_data, geo)))
      paired_data <- paired_data[order(paired_data$logFC, decreasing = TRUE), ]
      all_geos[!all_geos %in% unique(paired_data$GEO)]
      
      # change acute time for numeric
      paired_data$Training_weeks <- as.numeric(gsub("W", "", paired_data$Training_weeks))
      
      # Filter population of interest
      paired_data <- dplyr::filter(paired_data,
                                   GEO %in% input$human_studies_TA,
                                   predicted_sex %in% input$human_sex,
                                   predicted_ageGroup %in% input$human_age,
                                   predicted_fitness %in% input$human_fitness,
                                   predicted_BMIgroup %in% input$human_weight,
                                   Tissue %in% input$human_muscle,
                                   Health_group %in% input$human_disease,
                                   Training_weeks >= input$human_training_duration[1] & Training_weeks <= input$human_training_duration[2],
                                   predicted_TimeFromLastExerciseBout  >= input$human_training_biopsy[1] & predicted_TimeFromLastExerciseBout <= input$human_training_biopsy[2]
                                   )
      
      # Calculate meta-analysis
      meta_analysis_results <- fct_human_performMetaAnalysis(paired_data)
      
    }, error=function(e) NULL)
  }) 
  
  output$plot_human_TA <- renderPlot({ 
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculating", value = 1)
    
    # call data
    meta_analysis_results <- human_metaAnalysisData_TA()
    
    # Handle empty data case: this validates that paired_data is not null and has rows.
    shiny::validate(
      need(!is.null(meta_analysis_results$paired_data) && nrow(meta_analysis_results$paired_data) > 0, 
           "No studies found - try different selection criteria")
    )
    
    # Generate forest plot
    fct_human_generateForestPlot_training(meta = meta_analysis_results$meta, 
                             fdr = meta_analysis_results$fdr, 
                             paired_data = meta_analysis_results$paired_data,
                             color = "#E69F00")
  })
  
  output$human_forestPlotDynamic_TA <- renderUI({ 
    # call data
    meta_analysis_results <- human_metaAnalysisData_TA()
    
    # Determine height based on whether meta_analysis_results is NULL
    if (is.null(meta_analysis_results) || is.null(meta_analysis_results$paired_data)) {
      # If meta_analysis_results is NULL, set height to 50px
      height <- 50
    } else {
      # Calculate height based on the number of rows in paired_data
      height <- 150 + 17 * nrow(meta_analysis_results$paired_data)
    }
    
    plotOutput("plot_human_TA", height = paste0(height, "px"))
  })
  
  
  #=======================================================================================
  #
  # Forest plot - Human Combined Training
  #
  #=======================================================================================
  human_metaAnalysisData_TC <- reactive({
    tryCatch({  
      # select acute aerobic studies
      selected_meta <- human_meta[!is.na(human_meta$Training_weeks),]
      selected_meta <- selected_meta[!grepl("H0", selected_meta$Acute_exercise_hours),]
      selected_meta <- selected_meta[selected_meta$Protocol_group %in% c("Basal", "Combined"),]
      
      # select corresponding matrix
      selected_matrix <- human_matrix[,selected_meta$SampleID]
      
      # selected gene
      selected_gene <- "NR4A3"
      selected_gene <- input$genename_metaanalysis_human
      
      # merge data
      selected_data <- data.frame(selected_meta[,c("GEO", 
                                                   "predicted_sex", "predicted_ageGroup", "predicted_fitness", "predicted_BMIgroup", "Tissue", "Health_group",
                                                   "Training_weeks", "predicted_TimeFromLastExerciseBout")],
                                  gene_data = as.numeric(selected_matrix[selected_gene,]))
      
      # Apply the function to each GEO in the dataset
      all_geos <- unique(selected_data$GEO)
      paired_data <- do.call(rbind, lapply(all_geos, function(geo) fct_human_ttests_training(selected_data, geo)))
      paired_data <- paired_data[order(paired_data$logFC, decreasing = TRUE), ]
      all_geos[!all_geos %in% unique(paired_data$GEO)]
      
      # change acute time for numeric
      paired_data$Training_weeks <- as.numeric(gsub("W", "", paired_data$Training_weeks))
      
      # Filter population of interest
      paired_data <- dplyr::filter(paired_data,
                                   GEO %in% input$human_studies_TC,
                                   predicted_sex %in% input$human_sex,
                                   predicted_ageGroup %in% input$human_age,
                                   predicted_fitness %in% input$human_fitness,
                                   predicted_BMIgroup %in% input$human_weight,
                                   Tissue %in% input$human_muscle,
                                   Health_group %in% input$human_disease,
                                   Training_weeks >= input$human_training_duration[1] & Training_weeks <= input$human_training_duration[2]
                                   #predicted_TimeFromLastExerciseBout %in% input$human_training_biopsy
      )
      
      # Calculate meta-analysis
      meta_analysis_results <- fct_human_performMetaAnalysis(paired_data)
      
    }, error=function(e) NULL)
  }) 
  
  output$plot_human_TC <- renderPlot({ 
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculating", value = 1)
    
    # call data
    meta_analysis_results <- human_metaAnalysisData_TC()
    
    # Handle empty data case: this validates that paired_data is not null and has rows.
    shiny::validate(
      need(!is.null(meta_analysis_results$paired_data) && nrow(meta_analysis_results$paired_data) > 0, 
           "No studies found - try different selection criteria")
    )
    
    # Generate forest plot
    fct_human_generateForestPlot_training(meta = meta_analysis_results$meta, 
                                      fdr = meta_analysis_results$fdr, 
                                      paired_data = meta_analysis_results$paired_data,
                                      color = "#E69F00")
  })
  
  output$human_forestPlotDynamic_TC <- renderUI({ 
    # call data
    meta_analysis_results <- human_metaAnalysisData_TC()
    
    # Determine height based on whether meta_analysis_results is NULL
    if (is.null(meta_analysis_results) || is.null(meta_analysis_results$paired_data)) {
      # If meta_analysis_results is NULL, set height to 50px
      height <- 50
    } else {
      # Calculate height based on the number of rows in paired_data
      height <- 150 + 17 * nrow(meta_analysis_results$paired_data)
    }
    
    plotOutput("plot_human_TC", height = paste0(height, "px"))
  })
  
  #=======================================================================================
  #
  # Forest plot - Human HIT Training
  #
  #=======================================================================================
  human_metaAnalysisData_TH <- reactive({
    tryCatch({  
      # select acute aerobic studies
      selected_meta <- human_meta[!is.na(human_meta$Training_weeks),]
      selected_meta <- selected_meta[!grepl("H0", selected_meta$Acute_exercise_hours),]
      selected_meta <- selected_meta[selected_meta$Protocol_group %in% c("Basal", "HIT"),]
      
      # select corresponding matrix
      selected_matrix <- human_matrix[,selected_meta$SampleID]
      
      # selected gene
      selected_gene <- "NR4A3"
      selected_gene <- input$genename_metaanalysis_human
      
      # merge data
      selected_data <- data.frame(selected_meta[,c("GEO", 
                                                   "predicted_sex", "predicted_ageGroup", "predicted_fitness", "predicted_BMIgroup", "Tissue", "Health_group",
                                                   "Training_weeks", "predicted_TimeFromLastExerciseBout")],
                                  gene_data = as.numeric(selected_matrix[selected_gene,]))
      
      # Apply the function to each GEO in the dataset
      all_geos <- unique(selected_data$GEO)
      paired_data <- do.call(rbind, lapply(all_geos, function(geo) fct_human_ttests_training(selected_data, geo)))
      paired_data <- paired_data[order(paired_data$logFC, decreasing = TRUE), ]
      all_geos[!all_geos %in% unique(paired_data$GEO)]
      
      # change acute time for numeric
      paired_data$Training_weeks <- as.numeric(gsub("W", "", paired_data$Training_weeks))
      
      # Filter population of interest
      paired_data <- dplyr::filter(paired_data,
                                   GEO %in% input$human_studies_TH,
                                   predicted_sex %in% input$human_sex,
                                   predicted_ageGroup %in% input$human_age,
                                   predicted_fitness %in% input$human_fitness,
                                   predicted_BMIgroup %in% input$human_weight,
                                   Tissue %in% input$human_muscle,
                                   Health_group %in% input$human_disease,
                                   Training_weeks >= input$human_training_duration[1] & Training_weeks <= input$human_training_duration[2]
                                   #predicted_TimeFromLastExerciseBout %in% input$human_training_biopsy
      )
      
      # Calculate meta-analysis
      meta_analysis_results <- fct_human_performMetaAnalysis(paired_data)
      
    }, error=function(e) NULL)
  }) 
  
  output$plot_human_TH <- renderPlot({ 
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculating", value = 1)
    
    # call data
    meta_analysis_results <- human_metaAnalysisData_TH()
    
    # Handle empty data case: this validates that paired_data is not null and has rows.
    shiny::validate(
      need(!is.null(meta_analysis_results$paired_data) && nrow(meta_analysis_results$paired_data) > 0, 
           "No studies found - try different selection criteria")
    )
    
    # Generate forest plot
    fct_human_generateForestPlot_training(meta = meta_analysis_results$meta, 
                                      fdr = meta_analysis_results$fdr, 
                                      paired_data = meta_analysis_results$paired_data,
                                      color = "#E69F00")
  })
  
  output$human_forestPlotDynamic_TH <- renderUI({ 
    # call data
    meta_analysis_results <- human_metaAnalysisData_TH()
    
    # Determine height based on whether meta_analysis_results is NULL
    if (is.null(meta_analysis_results) || is.null(meta_analysis_results$paired_data)) {
      # If meta_analysis_results is NULL, set height to 50px
      height <- 50
    } else {
      # Calculate height based on the number of rows in paired_data
      height <- 150 + 17 * nrow(meta_analysis_results$paired_data)
    }
    
    plotOutput("plot_human_TH", height = paste0(height, "px"))
  })
  
  
  #=======================================================================================
  #
  # Forest plot - Human Resistance Training
  #
  #=======================================================================================
  human_metaAnalysisData_TR <- reactive({
    tryCatch({  
      # select acute aerobic studies
      selected_meta <- human_meta[!is.na(human_meta$Training_weeks),]
      selected_meta <- selected_meta[!grepl("H0", selected_meta$Acute_exercise_hours),]
      selected_meta <- selected_meta[selected_meta$Protocol_group %in% c("Basal", "Resistance"),]
      
      # select corresponding matrix
      selected_matrix <- human_matrix[,selected_meta$SampleID]
      
      # selected gene
      selected_gene <- "NR4A3"
      selected_gene <- input$genename_metaanalysis_human
      
      # merge data
      selected_data <- data.frame(selected_meta[,c("GEO", 
                                                   "predicted_sex", "predicted_ageGroup", "predicted_fitness", "predicted_BMIgroup", "Tissue", "Health_group",
                                                   "Training_weeks", "predicted_TimeFromLastExerciseBout")],
                                  gene_data = as.numeric(selected_matrix[selected_gene,]))
      
      # Apply the function to each GEO in the dataset
      all_geos <- unique(selected_data$GEO)
      paired_data <- do.call(rbind, lapply(all_geos, function(geo) fct_human_ttests_training(selected_data, geo)))
      paired_data <- paired_data[order(paired_data$logFC, decreasing = TRUE), ]
      all_geos[!all_geos %in% unique(paired_data$GEO)]
      
      # change acute time for numeric
      paired_data$Training_weeks <- as.numeric(gsub("W", "", paired_data$Training_weeks))
      
      # Filter population of interest
      paired_data <- dplyr::filter(paired_data,
                                   GEO %in% input$human_studies_TR,
                                   predicted_sex %in% input$human_sex,
                                   predicted_ageGroup %in% input$human_age,
                                   predicted_fitness %in% input$human_fitness,
                                   predicted_BMIgroup %in% input$human_weight,
                                   Tissue %in% input$human_muscle,
                                   Health_group %in% input$human_disease,
                                   Training_weeks >= input$human_training_duration[1] & Training_weeks <= input$human_training_duration[2]
                                   #predicted_TimeFromLastExerciseBout %in% input$human_training_biopsy
      )
      
      # Calculate meta-analysis
      meta_analysis_results <- fct_human_performMetaAnalysis(paired_data)
      
    }, error=function(e) NULL)
  }) 
  
  output$plot_human_TR <- renderPlot({ 
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculating", value = 1)
    
    # call data
    meta_analysis_results <- human_metaAnalysisData_TR()
    
    # Handle empty data case: this validates that paired_data is not null and has rows.
    shiny::validate(
      need(!is.null(meta_analysis_results$paired_data) && nrow(meta_analysis_results$paired_data) > 0, 
           "No studies found - try different selection criteria")
    )
    
    # Generate forest plot
    fct_human_generateForestPlot_training(meta = meta_analysis_results$meta, 
                                      fdr = meta_analysis_results$fdr, 
                                      paired_data = meta_analysis_results$paired_data,
                                      color = "#E69F00")
  })
  
  output$human_forestPlotDynamic_TR <- renderUI({ 
    # call data
    meta_analysis_results <- human_metaAnalysisData_TR()
    
    # Determine height based on whether meta_analysis_results is NULL
    if (is.null(meta_analysis_results) || is.null(meta_analysis_results$paired_data)) {
      # If meta_analysis_results is NULL, set height to 50px
      height <- 50
    } else {
      # Calculate height based on the number of rows in paired_data
      height <- 150 + 17 * nrow(meta_analysis_results$paired_data)
    }
    
    plotOutput("plot_human_TR", height = paste0(height, "px"))
  })

  
  #=======================================================================================
  #
  # Overview plot - Human
  #
  #=======================================================================================
  human_overviewData <- reactive({
    
    withProgress(message = 'Calculating meta-analyses', value = 1, max=11, {

    # Helper function to create a data frame or NA values if the input is NULL
    create_data_frame_or_na <- function(data, study_name, color) {
      if (is.null(data)) {
        return(data.frame(
          Study = study_name,
          beta = NA,
          ci.lb = NA,
          ci.ub = NA,
          fdr = NA,
          color = color
        ))
      } else {
        return(data.frame(
          Study = study_name,
          beta = data$meta$beta,
          ci.lb = data$meta$ci.lb,
          ci.ub = data$meta$ci.ub,
          fdr = data$fdr,
          color = color
        ))
      }
    }
    
    # Collect data from the different reactives, returning NA values if any dataset is NULL
    incProgress(1, detail="Acute aerobic exercise")
    human_metaAnalysisData_AA <- create_data_frame_or_na(human_metaAnalysisData_AA(), "Acute\nAerobic", "#E69F00")
    incProgress(1, detail="Acute HIT exercice")
    human_metaAnalysisData_AH <- create_data_frame_or_na(human_metaAnalysisData_AH(), "Acute\nHIT", "#E69F00")
    incProgress(1, detail="Acute resistance exercice")
    human_metaAnalysisData_AR <- create_data_frame_or_na(human_metaAnalysisData_AR(), "Acute\nResistance", "#E69F00")
    incProgress(1, detail="Inactivity")
    human_metaAnalysisData_IN <- create_data_frame_or_na(human_metaAnalysisData_IN(), "Inactivity", "#CC79A7")
    incProgress(1, detail="Aerobic training")
    human_metaAnalysisData_TA <- create_data_frame_or_na(human_metaAnalysisData_TA(), "Aerobic\nTraining", "#0072B2")
    incProgress(1, detail="Combined training")
    human_metaAnalysisData_TC <- create_data_frame_or_na(human_metaAnalysisData_TC(), "Combined\nTraining", "#0072B2")
    incProgress(1, detail="HIT training")
    human_metaAnalysisData_TH <- create_data_frame_or_na(human_metaAnalysisData_TH(), "HIT\nTraining", "#0072B2")
    incProgress(1, detail="Resistance training")
    human_metaAnalysisData_TR <- create_data_frame_or_na(human_metaAnalysisData_TR(), "Resistance\nTraining", "#0072B2")
    
    # Join everything in one table - full_join should now work as expected
    incProgress(1, detail="Making figure")
    overview_data <- full_join(human_metaAnalysisData_AA, human_metaAnalysisData_AH)
    overview_data <- full_join(overview_data, human_metaAnalysisData_AR)
    overview_data <- full_join(overview_data, human_metaAnalysisData_IN)
    overview_data <- full_join(overview_data, human_metaAnalysisData_TA)
    overview_data <- full_join(overview_data, human_metaAnalysisData_TC)
    overview_data <- full_join(overview_data, human_metaAnalysisData_TH)
    overview_data <- full_join(overview_data, human_metaAnalysisData_TR)

    # order dataset
    overview_data$Study <- factor(overview_data$Study, levels = c(
      "Acute\nAerobic", "Acute\nResistance", "Acute\nHIT", "Inactivity", "Aerobic\nTraining", "Resistance\nTraining", "Combined\nTraining", "HIT\nTraining"
    ))
    
    # round fdr
    overview_data$fdr <- format(overview_data$fdr, scientific = TRUE, digits = 2)
    return(overview_data)
    })
  })
    
    
  output$human_overviewPlot <- renderPlot({
    
    overview_data <- human_overviewData()
    selected_gene <- input$genename_metaanalysis_human
    
    # check that data is available
    shiny::validate(need(sum(is.na(overview_data$beta)) != length(overview_data$beta),
                  "No studies found - try different selection criteria"))
    
    # prepare position for stats
    overview_data$y.position <- overview_data$ci.ub + 0.5
    overview_data$y.position <- ifelse(is.na(overview_data$y.position), -1, overview_data$y.position)
    
    # Create the plot using base ggplot2
    ggplot(overview_data, aes(x = Study, y = beta, fill = Study)) +
      geom_col(position = position_dodge(), width = 0.8, 
               size = 0.2, color = "black") +                       # black border on bars
      geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
                    width = 0.2, position = position_dodge(.9)) +
      geom_hline(yintercept = 0) +
      labs(
        x = NULL,
        y = paste0(selected_gene, "\nMeta-analysis score (logFC)")
      ) +
      theme_bw() + theme_ggplot +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
      scale_fill_viridis_d() +
      geom_text(aes(y = y.position, label = fdr), size = 3)
    
    })
  
  
  
  #=======================================================================================
  #
  # Timeline plot - Human acute exercise
  #
  #=======================================================================================
  output$human_timelinePlot_acute <- renderPlot({
    #show progress bar
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculating", value = 1)
    
    gene_of_interest <- "NR4A3"
    gene_of_interest <- input$genename_metaanalysis_human
    
    acute_metadata <- human_meta[!is.na(human_meta$Acute_exercise_hours),]
    acute_metadata <- acute_metadata[!grepl("W", acute_metadata$Training_weeks),]
    acute_matrix <- human_matrix[,acute_metadata$SampleID] 
    
    mydata <- data.frame(acute_metadata, 
                         gene_data = as.numeric(acute_matrix[gene_of_interest,]))
    
    # normalize to median of PRE
    mydata$gene_data <- mydata$gene_data - median(mydata[mydata$Acute_exercise_hours == "PRE",]$gene_data, na.rm = T)
    unique(mydata$Acute_exercise_hours)
    mydata$bins <- "PRE"
    mydata$bins[mydata$Acute_exercise_hours %in% c("H00", "H01")] <- "0 - 1"
    mydata$bins[mydata$Acute_exercise_hours %in% c("H03", "H04")] <- "3 - 4"
    mydata$bins[mydata$Acute_exercise_hours %in% c("H05", "H06", "H08")] <- "5 - 8"
    mydata$bins[mydata$Acute_exercise_hours %in% c("H18", "H24")] <- "18 - 24"
    mydata$bins[mydata$Acute_exercise_hours %in% c("H48", "H96")] <- "48 - 96"
    mydata$bins <- factor(mydata$bins, levels = c("PRE", "0 - 1", "3 - 4", "5 - 8", "18 - 24", "48 - 96"))
    
    # plot
    ggplot(mydata, aes(x=bins, y=gene_data, fill=bins)) + 
      geom_hline(yintercept=0, color="gray50", linetype="dashed") +
      geom_boxplot() +
      labs(x="Time after exercise (h)",
           y=paste(gene_of_interest, ", log2(FC)", sep="")) +
      theme_bw() + theme_ggplot +
      scale_fill_viridis_d()
  })
  
  output$human_timelineTable_acute <- renderTable(rownames = TRUE, align='c', { 
    res <- human_timelineStats_acute["NR4A3",]
    res <- human_timelineStats_acute[input$genename_metaanalysis_human,]
    
    res <- rbind(as.numeric(res[,1:3]),
                 as.numeric(res[grepl("TTest.0_1", colnames(res))]),
                 as.numeric(res[grepl("TTest.3_4", colnames(res))]),
                 as.numeric(res[grepl("TTest.5_8", colnames(res))]),
                 as.numeric(res[grepl("TTest.18_24", colnames(res))]),
                 as.numeric(res[grepl("TTest.48_96", colnames(res))]))
    
    res <- data.frame(res,
                      row.names = c("F-test", "0-1h vs Pre", "2-4h vs Pre", "5-8h vs Pre", "18-24h vs Pre", "48-96h vs Pre"))
    colnames(res) <- c('estimate', 'p.value', 'adj.P.Val')
    
    # format values
    res$estimate <- round(res$estimate, 2)
    res$p.value <- scientific(res$p.value, 1)
    res$adj.P.Val <- scientific(res$adj.P.Val, 1)

    return(res)
  })
  
  
  #=======================================================================================
  #
  # Timeline plot - Human inactivity
  #
  #=======================================================================================
  output$human_timelinePlot_inactivity      <- renderPlot({
    #show progress bar
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculating", value = 1)
    
    gene_of_interest <- "NR4A3"
    gene_of_interest <- input$genename_metaanalysis_human
    
    metadata <- human_meta[!is.na(human_meta$Inactivity_days),]
    metadata <- metadata[!grepl("RELO", metadata$Inactivity_days),]
    matrix <- human_matrix[,metadata$SampleID] 
    
    mydata <- data.frame(metadata, 
                         gene_data = as.numeric(matrix[gene_of_interest,]))
    
    # normalize to median of PRE
    mydata$gene_data <- mydata$gene_data - median(mydata[mydata$Inactivity_days == "PRE",]$gene_data, na.rm = T)
    mydata$weeks <- "PRE"
    mydata$weeks[mydata$Inactivity_days %in% c("D02", "D05", "D10")] <- "< 2"
    mydata$weeks[mydata$Inactivity_days %in% c("D14", "D21", "D35")] <- "2 - 5"
    mydata$weeks[mydata$Inactivity_days %in% c("D60", "D70", "D84")] <- "> 5"
    mydata$weeks <- factor(mydata$weeks, levels = c("PRE", "< 2", "2 - 5", "> 5"))
    
    # plot
    ggplot(mydata, aes(x=weeks, y=gene_data, fill=weeks)) + 
      theme_bw() + theme_ggplot +
      geom_hline(yintercept=0, color="gray50", linetype="dashed") +
      geom_boxplot() +
      labs(x="Inactivity duration (weeks)",
           y=paste(gene_of_interest, ", log2(fold-change)", sep="")) +
      scale_fill_viridis_d()
  })
  
  output$human_timelineTable_inactivity <- renderTable(rownames = TRUE, align='c', { 
    res <- human_timelineStats_inactivity["NR4A3",]
    res <- human_timelineStats_inactivity[input$genename_metaanalysis_human,]
    
    res <- rbind(as.numeric(res[,1:3]),
                 as.numeric(res[grepl("TTest.0_1", colnames(res))]),
                 as.numeric(res[grepl("TTest.1_4", colnames(res))]),
                 as.numeric(res[grepl("TTest.5_12", colnames(res))]))
    
    res <- data.frame(res,
                      row.names = c("F-test", "<1 weeks vs Pre", "2-5 weeks vs Pre", ">5 weeks vs Pre"))
    colnames(res) <- c('estimate', 'p.value', 'adj.P.Val')
    
    # format values
    res$estimate <- round(res$estimate, 2)
    res$p.value <- scientific(res$p.value, 1)
    res$adj.P.Val <- scientific(res$adj.P.Val, 1)
    
    return(res)
  })
  
  
  #=======================================================================================
  #
  # Co-expression plot - Human 
  #
  #=======================================================================================
  
  #subset data for correlations
  human_corrData <- eventReactive(input$updateCorrHuman, {
    
    withProgress(message = 'Calculating correlations:', value = 1, max=10, {
      
    incProgress(1, detail="gathering selected data")
    
    corr_meta <- human_meta
    
    # Extract GEO and data_category columns from corr_meta
    data_category <- paste(corr_meta$sample_type, corr_meta$Protocol_group)
    data_category <- gsub("Basal.*", "Basal", data_category)
    data_category <- gsub("Inactivity.*", "Inactivity", data_category)
    corr_meta$data_category <- data_category
    table(corr_meta$data_category)
    
    # filter data based on input
    corr_meta <- dplyr::filter(corr_meta,
                               data_category %in% c("Basal", input$human_corr_protocol))
    
    # get the matrix based on selected metadata
    corr_matrix <- human_matrix[,corr_meta$SampleID]
    all(colnames(corr_matrix) == corr_meta$SampleID)
    
    # Function to normalize a single gene row based on the GEO grouping
    normalize_gene <- function(gene_data, GEO, data_category) {
      data <- data.frame(Gene = gene_data, GEO = GEO, data_category = data_category)
      data <- data %>%
        group_by(GEO) %>%
        mutate(Gene = Gene - median(Gene[data_category %in% "Basal"], na.rm = TRUE)) %>%
        ungroup()
      return(data$Gene)
    }
    
    # Apply the normalization across all rows in corr_matrix
    incProgress(2, detail="normalizing data")
    normalized_corr_matrix <- apply(corr_matrix, 1, function(gene_data) {
      normalize_gene(gene_data, corr_meta$GEO, corr_meta$data_category)
    })
    
    # Convert the result back to the original matrix structure
    normalized_corr_matrix <- data.frame(t(normalized_corr_matrix))
    rownames(normalized_corr_matrix) <- rownames(corr_matrix)
    colnames(normalized_corr_matrix) <- colnames(corr_matrix)
    
    # remove basal
    corr_meta <- corr_meta[!corr_meta$data_category %in% "Basal",]
    normalized_corr_matrix <- normalized_corr_matrix[,corr_meta$SampleID]
    
    #select data for gene of interest
    genename <- "NR4A3"
    genename <- input$genename_metaanalysis_human
    geneofinterest <- as.numeric(normalized_corr_matrix[genename,])
    
    # Define a function to calculate Spearman correlation or p-value
    calculate_spearman <- function(data, type = c("estimate", "p.value")) {
      type <- match.arg(type)
      lapply(seq(1, nrow(data), by = 1000), function(start_idx) {
        end_idx <- min(start_idx + 999, nrow(data))
        apply(data[start_idx:end_idx, ], 1, function(x) {
          tryCatch({
            if (type == "estimate") {
              cor.test(x, geneofinterest, method = "spearman", exact = F)$estimate
            } else {
              cor.test(x, geneofinterest, method = "spearman", exact = F)$p.value
            }
          }, error = function(e) {
            NA
          })
        })
      }) %>%
        unlist()
    }
    
    # Calculate Spearman correlation coefficients and p-values
    incProgress(1, detail = "calculating spearman statistics")
    Spearman.r <- calculate_spearman(normalized_corr_matrix, type = "estimate")
    Spearman.p <- calculate_spearman(normalized_corr_matrix, type = "p.value")
    
    # Make table with statistics
    Spearman.adj.P.Val <- p.adjust(Spearman.p, method = "bonferroni")

    # make a table    
    coeff <- data.frame(Spearman.r, Spearman.adj.P.Val)
    colnames(coeff) <- c("Spearman.r", "adj.P.Val")
    coeff <- coeff[order(coeff$adj.P.Val), ]
    
    # round numbers
    coeff$Spearman.r <- round(coeff$Spearman.r, digits = 3)
    coeff$adj.P.Val <- format(coeff$adj.P.Val, scientific = T, digits = 2)
    
    # eclude auto-correlation
    coeff <- coeff[!rownames(coeff) %in% genename,]
    
    # save in a list
    incProgress(1, detail = "preparing tables")
    human_corrData <- list(data = normalized_corr_matrix,
                            metadata = corr_meta,
                            stats = coeff)
    return(human_corrData)
    })
  })

  
  #make table output
  output$human_corrTable <- DT::renderDataTable(escape = FALSE, 
                                                 rownames = T, 
                                                 selection = "single", {
                                                   shiny::validate(need(!is.null(human_corrData()),   
                                                                 "Start by selecting a gene in the list of official gene names"))
                                                   
                                                   human_corrData <- human_corrData()
                                                   human_corrData <- human_corrData$stats
                                                 })
  
  # Make plot output
  output$human_corrPlot  <- renderPlot({ 
    shiny::validate(need(!is.null(human_corrData()),     " "))
    shiny::validate(need(input$human_corrTable_rows_selected!="",  "Click on a gene in the table to display the correlation")) 
    
    human_corrData <- human_corrData()

    corr_meta <- human_corrData$metadata
    corr_matrix <- human_corrData$data
    
    #collect names of the 2 genes to correlate
    Gene1name <- "NR4A3"
    Gene1name <- input$genename_metaanalysis_human
    Gene2name <- "PPARGC1A"
    Gene2name <- rownames(human_corrData$stats[input$human_corrTable_rows_selected,])
    
    #find data from the 2 genes and merge
    data  <- data.frame(corr_meta, 
                        Gene1 = as.numeric(corr_matrix[Gene1name,]),
                        Gene2 = as.numeric(corr_matrix[Gene2name,]))

    # plot    
    active <- ggplot(data, aes(x=Gene2, y=Gene1)) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_smooth(method=lm, se=F, fullrange=TRUE, size=0.75, color="black") +
      geom_point(aes(x=Gene2, y=Gene1, 
                     color=data_category,
                     shape=data_category),
                 size=1) +
      labs(x=paste(Gene2name, ", log2(relative to basal)", sep=""),
           y=paste(Gene1name, ", log2(relative to basal)", sep="")) +
      theme_bw() + theme_ggplot + theme(legend.position="right") + 
      scale_shape_manual(values=rep(c(15,0,16,1,17,2),10)) +
      scale_color_viridis_d()
  active
    return(active) 
  })
  
  # Make correlation description
  output$human_corrDescription <- renderText({
    req(input$human_corrTable_rows_selected)
    symbol <- input$human_corrTable_rows_selected
    
    # Step 1: Query gene symbol
    res <- GET(paste0("https://mygene.info/v3/query?q=", symbol, "&species=human"))
    data <- fromJSON(content(res, "text", encoding = "UTF-8"))
    
    # Safety check: ensure we got hits
    if (is.null(data$hits) || length(data$hits) == 0) {
      return("No gene information found.")
    }
    
    # Try to extract entrezgene ID from first hit
    gene_id <- data$hits[[1]][1]
    if (is.null(gene_id)) {
      return("Gene ID not found for this symbol.")
    }
    
    # Step 2: Get full gene info
    res2 <- GET(paste0("https://mygene.info/v3/gene/", gene_id))
    info <- fromJSON(content(res2, "text", encoding = "UTF-8"))
    
    # Return summary, or fallback
    if (!is.null(info$summary)) {
      paste("Description:", info$summary)
    } else {
      "No summary available for this gene."
    }
  })
  

  
  output$human_corrGeneCardsLink <- renderUI({
    shiny::validate(need(!is.null(human_corrData()),     " "))
    shiny::validate(need(input$human_corrTable_rows_selected!="",  " ")) 
    
    human_corrData <- human_corrData()
    
    #find gene selected in the table
    symbol_from_table <- "PPARGC1A"
    symbol_from_table <- rownames(human_corrData$stats[input$human_corrTable_rows_selected,])
    
    #Make link to genecard
    GeneCards <- paste("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", symbol_from_table, sep="")
    GeneCards <- sprintf(paste0('<a href="', paste("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", symbol_from_table, sep=""),
                                '" target="_blank">',
                                'Learn more about ', symbol_from_table, ' on GeneCards' ,'</a>'))
    GeneCards <- HTML(GeneCards)
    return(GeneCards)
  })
  

 
  #=======================================================================================
  #
  # Forest plots - Mouse acute
  #
  #=======================================================================================
  
  mouse_metaAnalysisData_AA <- reactive({
    tryCatch({  
      # select acute aerobic studies
      selected_meta <- mouse_meta[!is.na(mouse_meta$Acute_exercise_hours),]
      selected_meta <- selected_meta[!grepl("W", selected_meta$Training_weeks),]
      
      # select corresponding matrix
      selected_matrix <- mouse_matrix[,selected_meta$SampleID]
      
      # selected gene
      selected_gene <- "Nr4a3"
      selected_gene <- input$genename_metaanalysis_mouse
      
      # merge data
      colnames(selected_meta)
      selected_data <- data.frame(selected_meta[,c("GEO", 
                                                   "predicted_sex", "Age", "Genotype", "Tissue", "Treatment",
                                                   "Acute_exercise_hours")],
                                  gene_data = as.numeric(selected_matrix[selected_gene,]))
      
      # Apply the function to each GEO in the dataset
      all_geos <- unique(selected_data$GEO)
      paired_data <- do.call(rbind, lapply(all_geos, function(geo) fct_mouse_ttests(selected_data, geo, "Acute_exercise_hours")))
      paired_data <- paired_data[order(paired_data$logFC, decreasing = TRUE), ]
      
      # change acute time for numeric
      paired_data$Condition <- as.numeric(gsub("H", "", paired_data$Condition))
      
      # Filter population of interest
      paired_data <- dplyr::filter(paired_data,
                                   GEO %in% input$mouse_studies_AA,
                                   predicted_sex %in% input$mouse_sex,
                                   Age >= input$mouse_age[1] & Age <= input$mouse_age[2],
                                   #Protocol
                                   #Tissue %in% input$mouse_muscle,
                                   #Genotype %in% input$disease_choice,
                                   Condition >= input$mouse_acute_biopsy[1] & Condition <= input$mouse_acute_biopsy[2]
                                   )

      # Calculate meta-analysis
      meta_analysis_results <- fct_mouse_performMetaAnalysis(paired_data)
      
    }, error=function(e) NULL)
  }) 
  
  output$mouse_forestPlot_AA <- renderPlot({ 
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculating", value = 1)
    
    # call data
    meta_analysis_results <- mouse_metaAnalysisData_AA()
    
    # Handle empty data case: this validates that paired_data is not null and has rows.
    shiny::validate(
      need(!is.null(meta_analysis_results$paired_data) && nrow(meta_analysis_results$paired_data) > 0, 
           "No studies found - try different selection criteria")
    )
    
    # extract data
    meta <- meta_analysis_results$meta
    fdr <- meta_analysis_results$fdr
    
    # make tables
    tabledata <- data.frame(mean = c(NA, meta$data$logFC, NA, meta$beta), 
                            lower = c(NA, meta$data$CI.L, NA, meta$ci.lb),
                            upper = c(NA, meta$data$CI.R, NA, meta$ci.ub))
    
    tabletext <- cbind(c('Study', meta$data$GEO, NA, ""),
                       c('Sex', as.character(meta$data$predicted_sex), NA, ""),
                       c('Age', as.character(meta$data$Age), NA, ""),
                       c('Genotype', as.character(meta$data$Genotype), NA, ""),
                       c('Timepoint', as.character(meta$data$Condition), NA, ""),
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
               col = fpColors(box = "#E69F00", line = "black", summary = "#E69F00"),
               boxsize = c(NA, meta$data$POST_n / 20, NA, 1.5)
    )
  })
  
  output$mouse_forestPlotDynamic_AA <- renderUI({ 
    # call data
    meta_analysis_results <- mouse_metaAnalysisData_AA()
    
    # Determine height based on whether meta_analysis_results is NULL
    if (is.null(meta_analysis_results) || is.null(meta_analysis_results$paired_data)) {
      # If meta_analysis_results is NULL, set height to 50px
      height <- 50
    } else {
      # Calculate height based on the number of rows in paired_data
      height <- 150 + 17 * nrow(meta_analysis_results$paired_data)
    }
    
    plotOutput("mouse_forestPlot_AA", height = paste0(height, "px"))
  })
  

  #=======================================================================================
  #
  # Forest plots - Mouse inactivity
  #
  #=======================================================================================
  
  mouse_metaAnalysisData_IN <- reactive({
    tryCatch({  
      # select acute aerobic studies
      selected_meta <- mouse_meta[!is.na(mouse_meta$Inactivity_days),]
      selected_meta <- selected_meta[!selected_meta$Inactivity_days %in% c("RELOAD1", "RELOAD3.5H"),]
      
      # select corresponding matrix
      selected_matrix <- mouse_matrix[,selected_meta$SampleID]
      
      # selected gene
      selected_gene <- "Nr4a3"
      selected_gene <- input$genename_metaanalysis_mouse
      
      # merge data
      colnames(selected_meta)
      selected_data <- data.frame(selected_meta[,c("GEO", 
                                                   "predicted_sex", "Age", "Genotype", "Tissue", "Treatment",
                                                   "Inactivity_days")],
                                  gene_data = as.numeric(selected_matrix[selected_gene,]))
      
      # Apply the function to each GEO in the dataset
      all_geos <- unique(selected_data$GEO)
      paired_data <- do.call(rbind, lapply(all_geos, function(geo) fct_mouse_ttests(selected_data, geo, "Inactivity_days")))
      paired_data <- paired_data[order(paired_data$logFC, decreasing = TRUE), ]
      
      # change acute time for numeric
      paired_data$Condition <- as.numeric(gsub("D", "", paired_data$Condition))
      
      # Filter population of interest
      paired_data <- dplyr::filter(paired_data,
                                   GEO %in% input$mouse_studies_IN,
                                   predicted_sex %in% input$mouse_sex,
                                   Age >= input$mouse_age[1] & Age <= input$mouse_age[2],
                                   #Protocol
                                   #Tissue %in% input$mouse_muscle,
                                   #Genotype %in% input$disease_choice,
                                   Condition >= input$mouse_inactivity_duration[1] & Condition <= input$mouse_inactivity_duration[2]
      )
      
      # Calculate meta-analysis
      meta_analysis_results <- fct_mouse_performMetaAnalysis(paired_data)
      
    }, error=function(e) NULL)
  }) 
  
  output$mouse_forestPlot_IN <- renderPlot({ 
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculating", value = 1)
    
    # call data
    meta_analysis_results <- mouse_metaAnalysisData_IN()
    
    # Handle empty data case: this validates that paired_data is not null and has rows.
    shiny::validate(
      need(!is.null(meta_analysis_results$paired_data) && nrow(meta_analysis_results$paired_data) > 0, 
           "No studies found - try different selection criteria")
    )
    
    # extract data
    meta <- meta_analysis_results$meta
    fdr <- meta_analysis_results$fdr
    
    # make tables
    tabledata <- data.frame(mean = c(NA, meta$data$logFC, NA, meta$beta), 
                            lower = c(NA, meta$data$CI.L, NA, meta$ci.lb),
                            upper = c(NA, meta$data$CI.R, NA, meta$ci.ub))
    
    tabletext <- cbind(c('Study', meta$data$GEO, NA, ""),
                       c('Sex', as.character(meta$data$predicted_sex), NA, ""),
                       c('Age', as.character(meta$data$Age), NA, ""),
                       c('Genotype', as.character(meta$data$Genotype), NA, ""),
                       c('Duration', as.character(meta$data$Condition), NA, ""),
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
               col = fpColors(box = "#CC79A7", line = "black", summary = "#CC79A7"),
               boxsize = c(NA, meta$data$POST_n / 20, NA, 1.5)
    )
  })
  
  output$mouse_forestPlotDynamic_IN <- renderUI({ 
    # call data
    meta_analysis_results <- mouse_metaAnalysisData_IN()
    
    # Determine height based on whether meta_analysis_results is NULL
    if (is.null(meta_analysis_results) || is.null(meta_analysis_results$paired_data)) {
      # If meta_analysis_results is NULL, set height to 50px
      height <- 50
    } else {
      # Calculate height based on the number of rows in paired_data
      height <- 150 + 17 * nrow(meta_analysis_results$paired_data)
    }
    
    plotOutput("mouse_forestPlot_IN", height = paste0(height, "px"))
  })
  
  
  #=======================================================================================
  #
  # Forest plots - Mouse training
  #
  #=======================================================================================
  
  mouse_metaAnalysisData_TA <- reactive({
    tryCatch({  
      # select acute aerobic studies
      selected_meta <- mouse_meta[!is.na(mouse_meta$Training_weeks),]
      selected_meta <- selected_meta[!grepl("H", selected_meta$Acute_exercise_hours),]
      
      # select corresponding matrix
      selected_matrix <- mouse_matrix[,selected_meta$SampleID]
      
      # selected gene
      selected_gene <- "Nr4a3"
      selected_gene <- input$genename_metaanalysis_mouse
      
      # merge data
      colnames(selected_meta)
      selected_data <- data.frame(selected_meta[,c("GEO", 
                                                   "predicted_sex", "Age", "Genotype", "Tissue", "Treatment",
                                                   "Training_weeks")],
                                  gene_data = as.numeric(selected_matrix[selected_gene,]))
      
      # Apply the function to each GEO in the dataset
      all_geos <- unique(selected_data$GEO)
      paired_data <- do.call(rbind, lapply(all_geos, function(geo) fct_mouse_ttests(selected_data, geo, "Training_weeks")))
      paired_data <- paired_data[order(paired_data$logFC, decreasing = TRUE), ]
      
      # change acute time for numeric
      paired_data$Condition <- as.numeric(gsub("W", "", paired_data$Condition))
      
      # Filter population of interest
      paired_data <- dplyr::filter(paired_data,
                                   GEO %in% input$mouse_studies_TA,
                                   predicted_sex %in% input$mouse_sex,
                                   Age >= input$mouse_age[1] & Age <= input$mouse_age[2],
                                   #Protocol
                                   #Tissue %in% input$mouse_muscle,
                                   #Genotype %in% input$disease_choice,
                                   Condition >= input$mouse_training_duration[1] & Condition <= input$mouse_training_duration[2]
      )
      
      # Calculate meta-analysis
      meta_analysis_results <- fct_mouse_performMetaAnalysis(paired_data)
      
    }, error=function(e) NULL)
  }) 
  
  output$mouse_forestPlot_TA <- renderPlot({ 
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Calculating", value = 1)
    
    # call data
    meta_analysis_results <- mouse_metaAnalysisData_TA()
    
    # Handle empty data case: this validates that paired_data is not null and has rows.
    shiny::validate(
      need(!is.null(meta_analysis_results$paired_data) && nrow(meta_analysis_results$paired_data) > 0, 
           "No studies found - try different selection criteria")
    )
    
    # extract data
    meta <- meta_analysis_results$meta
    fdr <- meta_analysis_results$fdr
    
    # make tables
    tabledata <- data.frame(mean = c(NA, meta$data$logFC, NA, meta$beta), 
                            lower = c(NA, meta$data$CI.L, NA, meta$ci.lb),
                            upper = c(NA, meta$data$CI.R, NA, meta$ci.ub))
    
    tabletext <- cbind(c('Study', meta$data$GEO, NA, ""),
                       c('Sex', as.character(meta$data$predicted_sex), NA, ""),
                       c('Age', as.character(meta$data$Age), NA, ""),
                       c('Genotype', as.character(meta$data$Genotype), NA, ""),
                       c('Timepoint', as.character(meta$data$Condition), NA, ""),
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
               col = fpColors(box = "#0072B2", line = "black", summary = "#0072B2"),
               boxsize = c(NA, meta$data$POST_n / 20, NA, 1.5)
    )
  })
  
  output$mouse_forestPlotDynamic_TA <- renderUI({ 
    # call data
    meta_analysis_results <- mouse_metaAnalysisData_TA()
    
    # Determine height based on whether meta_analysis_results is NULL
    if (is.null(meta_analysis_results) || is.null(meta_analysis_results$paired_data)) {
      # If meta_analysis_results is NULL, set height to 50px
      height <- 50
    } else {
      # Calculate height based on the number of rows in paired_data
      height <- 150 + 17 * nrow(meta_analysis_results$paired_data)
    }
    
    plotOutput("mouse_forestPlot_TA", height = paste0(height, "px"))
  })
  

  #=======================================================================================
  #
  # Overview plot - Mouse
  #
  #=======================================================================================
  mouse_overviewData <- reactive({
    
    withProgress(message = 'Calculating meta-analyses', value = 1, max=11, {
      
      # Helper function to create a data frame or NA values if the input is NULL
      create_data_frame_or_na <- function(data, study_name, color) {
        if (is.null(data)) {
          return(data.frame(
            Study = study_name,
            beta = NA,
            ci.lb = NA,
            ci.ub = NA,
            fdr = NA,
            color = color
          ))
        } else {
          return(data.frame(
            Study = study_name,
            beta = data$meta$beta,
            ci.lb = data$meta$ci.lb,
            ci.ub = data$meta$ci.ub,
            fdr = data$fdr,
            color = color
          ))
        }
      }
      
      # Collect data from the different reactives, returning NA values if any dataset is NULL
      incProgress(1, detail="Acute exercise")
      mouse_metaAnalysisData_AA <- create_data_frame_or_na(mouse_metaAnalysisData_AA(), "Acute\nExercise", "#E69F00")
      incProgress(1, detail="Inactivity")
      mouse_metaAnalysisData_IN <- create_data_frame_or_na(mouse_metaAnalysisData_IN(), "Inactivity", "#CC79A7")
      incProgress(1, detail="Exercise training")
      mouse_metaAnalysisData_TA <- create_data_frame_or_na(mouse_metaAnalysisData_TA(), "Exercise\nTraining", "#0072B2")

      # Join everything in one table - full_join should now work as expected
      incProgress(1, detail="Making figure")
      overview_data <- full_join(mouse_metaAnalysisData_AA, mouse_metaAnalysisData_IN)
      overview_data <- full_join(overview_data, mouse_metaAnalysisData_TA)

      # order dataset
      overview_data$Study <- factor(overview_data$Study, levels = c(
        "Acute\nExercise", "Inactivity", "Exercise\nTraining"
      ))
      
      # round fdr
      overview_data$fdr <- format(overview_data$fdr, scientific = TRUE, digits = 2)
      return(overview_data)
    })
  })
  
  
  output$mouse_overviewPlot <- renderPlot({
    
    overview_data <- mouse_overviewData()
    selected_gene <- input$genename_metaanalysis_human
    
    # check that data is available
    shiny::validate(need(sum(is.na(overview_data$beta)) != length(overview_data$beta),
                  "No studies found - try different selection criteria"))
    
    # Assuming `overview_data` and `selected_gene` are defined
    overview_data$y.position <- overview_data$ci.ub + 0.5
    overview_data$y.position <- ifelse(is.na(overview_data$y.position), -1, overview_data$y.position)
    
    # Create the plot using base ggplot2
    ggplot(overview_data, aes(x = Study, y = beta, fill = Study)) +
      geom_col(position = position_dodge(), width = 0.8, 
               size = 0.2, color = "black") +                       # black border on bars
      geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), 
                    width = 0.2, position = position_dodge(.9)) +
      geom_hline(yintercept = 0) +
      labs(
        x = NULL,
        y = paste0(selected_gene, "\nMeta-analysis score (logFC)")
      ) +
      theme_bw() + theme_ggplot + 
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
      scale_fill_viridis_d() +
      geom_text(aes(y = y.position, label = fdr), size = 3)
    
  })
 
  
  #=======================================================================================
  #
  # Download buttons
  #
  #=======================================================================================
  output$human_downloadOverview <- downloadHandler(
    filename = function() { "MetaMEx_v4.2409_human.csv" },
    content = function(file) {
      human_overviewData <- human_overviewData()
      human_overviewData$Study <- gsub("\n", " ", human_overviewData$Study)
      write.csv(human_overviewData, file, row.names = F)
    })
  

  output$mouse_downloadOverview <- downloadHandler(
    filename = function() { "MetaMEx_v4.2409_mouse.csv" },
    content = function(file) {
      mouse_overviewData <- human_overviewData()
      mouse_overviewData$Study <- gsub("\n", " ", mouse_overviewData$Study)
      write.csv(mouse_overviewData, file, row.names = F)
    })



  #-------------------------------
  # HUMAN DOWNLOADS
  #-------------------------------
  
  output$download_human_metadata <- downloadHandler(
    filename = function() { "MetaMEx_Human_Metadata.csv" },
    content = function(file) {
      write.csv(human_meta, file, row.names = FALSE)
    }
  )
  
  output$download_human_genelist <- downloadHandler(
    filename = function() { "MetaMEx_Human_GeneList.csv" },
    content = function(file) {
      write.csv(human_genes, file, row.names = FALSE)
    }
  )
  
  output$download_human_matrix <- downloadHandler(
    filename = function() { "MetaMEx_Human_GeneMatrix.csv" },
    content = function(file) {
      write.csv(human_matrix, file, row.names = TRUE)
    }
  )
  
  #-------------------------------
  # MOUSE DOWNLOADS
  #-------------------------------
  
  output$download_mouse_metadata <- downloadHandler(
    filename = function() { "MetaMEx_Mouse_Metadata.csv" },
    content = function(file) {
      write.csv(mouse_meta, file, row.names = FALSE)
    }
  )
  
  output$download_mouse_genelist <- downloadHandler(
    filename = function() { "MetaMEx_Mouse_GeneList.csv" },
    content = function(file) {
      write.csv(mouse_genes, file, row.names = FALSE)
    }
  )
  
  output$download_mouse_matrix <- downloadHandler(
    filename = function() { "MetaMEx_Mouse_GeneMatrix.csv" },
    content = function(file) {
      write.csv(mouse_matrix, file, row.names = TRUE)
    }
  )
  
  
  #---------------------------------------------------------------------
  #hide loading page
  #Sys.sleep(2)
  shinyjs::hide("loading_page", anim = F, animType = "fade")
  shinyjs::show("main_content")
  
  #---------------------------------------------------------------------
  # Citations
  output$frame <- renderUI({
    my_test <- tags$iframe(src="https://app.dimensions.ai/discover/publication?and_subset_publication_citations=pub.1124285483",
                           height=600, width='100%')
    print(my_test)
    my_test
  })
}

