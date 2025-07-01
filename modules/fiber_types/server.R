#-----------------------------------------------------------------------
#
# Profile of skeletal muscle fibers
#
#----------------------------------------------------------------------

server <- function(input, output, session) {
  
  #-----------------------------------------------
  # Resize the iframe containing the app
  # Sends a custom message to the parent HTML to adjust height dynamically
  session$onFlushed(function() {
    session$sendCustomMessage("resizeFrame", list())
  }, once = FALSE)
  
  #-----------------------------------------------
  # Code to collect target name from loading page    
  observe({
    query <- parseQueryString(session$clientData$url_search)
    
    if (!is.null(query$target)) {
      selected_target <- toupper(query$target)  # normalize to uppercase
      updateSelectizeInput(session, "inputTarget",
                           choices = gene_to_file$TARGET,
                           selected = selected_target,
                           server = TRUE)
    } else {
      updateSelectizeInput(session, "inputTarget",
                           choices = gene_to_file$TARGET,
                           selected = c("LDHA", "LDHB", "MYH7", "MYH1", "MYH2"),  # default genes
                           server = TRUE)
    }
  })
  
  #-----------------------------------------------
  # Reset button functionality: resets all UI inputs to their default values
  observeEvent(input$resetInputs, {
    updateSelectizeInput(session, "inputTarget", selected = character(0))
  })
  
  #-----------------------------------------------
  observeEvent(input$fibers, {
    fibers <- sort(input$fibers)
    
    # Only allow comparisons if exactly two fiber types are selected
    if (length(fibers) == 2) {
      # Create the appropriate comparison string
      comp <- paste(fibers, collapse = " vs ")
      
      # Update the selectInput with only this comparison and preselect it
      updateSelectInput(session, "inputComparison",
                        choices = comp,
                        selected = comp)
    } else {
      # If not exactly two selected, allow all options
      updateSelectInput(session, "inputComparison",
                        choices = c("Type I vs Type IIA", "Type I vs Type IIX", "Type IIA vs Type IIX"),
                        selected = "Type I vs Type IIA")
    }
  })
  
  #-----------------------------------------------------------------
  # REACTIVE: load only selected gene(s) from dataset
  selectedTargetData <- reactive({
    req(input$inputTarget, input$inputOmics)
    targetname <- c("LDHA", "LDHB", "MYH7", "MYH1", "MYH2")
    omics <- "Proteome"
    
    targetname <- toupper(unlist(strsplit(input$inputTarget, "[,;\\s]+")))
    omics <- input$inputOmics
    
    # Match gene to file and row info - gene and OMICS
    file_row <- gene_to_file %>%
      filter(TARGET %in% targetname & OMICS %in% omics)
    
    # Ensure there's something to analyze
    shiny::validate(
      need(nrow(file_row) > 0, "No data available for the selected filters. Please adjust your selections.")
    )
    
    # Split by file
    split_rows <- split(file_row, file_row$file)
    
    # Load selected genes from each file using lazy loading via Arrow
    selected_list <- lapply(split_rows, function(rows) {
      path <- file.path("data", unique(rows$file))  # Path to parquet dataset
      ds <- arrow::open_dataset(path, format = "parquet")  # Lazy load dataset
      selected <- ds %>%
        filter(TARGET %in% rows$TARGET) %>%
        collect() %>%
        as.data.frame()
      
      # Reorder rows to match input gene order
      rownames(selected) <- selected$TARGET
      selected <- selected[match(rows$TARGET, rownames(selected)), , drop = FALSE]
      selected <- selected[, !(colnames(selected) %in% "TARGET"), drop = FALSE]
      return(selected)
    })
    
    # Combine rows across files
    selected_row <- do.call(rbind, selected_list)
    rownames(selected_row) <- file_row$TARGET
    
    # Handle genes not found: fill with NA rows
    missing_targets <- setdiff(targetname, rownames(selected_row))
    if (length(missing_targets) > 0) {
      missing_df <- matrix(NA, nrow = length(missing_targets), ncol = ncol(selected_row))
      rownames(missing_df) <- missing_targets
      colnames(missing_df) <- colnames(selected_row)
      selected_row <- rbind(selected_row, missing_df)
    }
    
    # Ensure output preserves original input order
    df <- selected_row[targetname, , drop = FALSE]
    df <- data.frame(sample_id = colnames(df), t(df))
    
    # Merge with metadata
    dat <- right_join(metadata, df) 
    
    # Convert to long format for stats
    dat <- pivot_longer(dat, cols = c(5:ncol(dat)),
                        values_to = "y",
                        names_to = "Target")
    
    # Apply same filters as for plot
    dat <- dplyr::filter(dat,
                         fiber_type %in% input$fibers)
  })
  
  
  #-----------------------------------------------------------------
  # Boxplot
  boxplotTarget <- reactive({
    
    dat <- selectedTargetData() %>%
      filter(is.finite(y))

    dat$fiber_type <- factor(dat$fiber_type,
                             levels = c("Type I", "Type IIA", "Type IIX", "Mixed"))
    
    # Plot
    ggplot(dat, aes(x = Target, y = y, fill = fiber_type)) +
      geom_boxplot(position = position_dodge(0.8), outlier.size = 0) +
      #geom_sina(size = 0.5, position = position_dodge(0.8), alpha = 0.5) +
      theme_bw(base_size = 16) +
      labs(x = NULL, 
           y = "Relative expression, log2",
           fill = NULL) +
      scale_y_continuous(expand = c(0, 4)) +
      scale_fill_manual(values = c("Type I" = "#8B0000", 
                                   "Type IIA" = "#F5DEB3", 
                                   "Type IIX"= "#D3D3D3", 
                                   "Mixed" = "#A0522D")) 
    
    
  })
  
  output$TargetPlot <- renderPlot({
    boxplotTarget()
  })
  
  #-----------------------------------------------------------------
  # Reactive text description of the current filters (used in UI summary)
  filterSummary <- reactive({
    
    dat <- selectedTargetData()
    
    if (nrow(dat) == 0) {
      return("No data available for the selected filters.")
    }
    
    # Unique fiber types in the filtered dataset
    fiber_types <- sort(unique(dat$fiber_type))
    
    # Format readable list of fiber types
    format_list <- function(x) {
      x <- sort(unique(x))
      n <- length(x)
      if (n == 0) return("")
      if (n == 1) return(x)
      if (n == 2) return(paste(x, collapse = " and "))
      paste(paste(x[-n], collapse = ", "), "and", x[n])
    }
    
    fiber_text <- format_list(fiber_types)
    
    # Early exit if less than 2 types
    if (length(fiber_types) < 2) {
      return(glue::glue(
        "Only {fiber_text} muscle fibers are selected. At least two fiber types are required for statistical comparisons."
      ))
    }
    
    # Flags for comparison logic
    has_I   <- "Type I"   %in% fiber_types
    has_IIA <- "Type IIA" %in% fiber_types
    has_IIX <- "Type IIX" %in% fiber_types
    
    # Generate sentence about what comparisons are possible
    if (all(c("Type I", "Type IIA", "Type IIX") %in% fiber_types)) {
      comparison_sentence <- "Statistics are available for comparing Type I to Type IIA, Type I to Type IIX, and Type IIA to Type IIX fibers."
    } else {
      comparison_text <- c()
      if (has_I & has_IIA) comparison_text <- c(comparison_text, "Type I to Type IIA")
      if (has_I & has_IIX) comparison_text <- c(comparison_text, "Type I to Type IIX")
      if (has_IIA & has_IIX) comparison_text <- c(comparison_text, "Type IIA to Type IIX")
      
      if (length(comparison_text) == 0) {
        comparison_sentence <- "No valid pairwise comparisons can be performed with the current fiber selection."
      } else {
        comparison_sentence <- paste(
          "Statistics are available for comparing",
          paste(comparison_text, collapse = " or "),
          "fibers."
        )
      }
    }
    
    # Final summary
    glue::glue(
      "Based on your selection, the plot and statistics reflect data from {fiber_text} muscle fibers. {comparison_sentence}"
    )
  })
  
  output$filterSummaryText <- renderText({
    filterSummary()
  })
  
  
  #-----------------------------------------------------------------
  # Compute statistics (Wilcoxon test + summary)
  statisticsData_IvsIIA <- reactive({
    # Validate at least two fibers are selected
    validate(
      need(length(input$fibers) >= 2,
           "Please select at least two fiber types to calculate statistics.")
    )
    
    dat <- selectedTargetData()

    # Group and compute statistics per gene - Type I vs Type IIA
    stats_result_IvsIIA <- dat %>%
      filter(fiber_type %in% c("Type I", "Type IIA")) 

    # Make table
    stats_result_IvsIIA <- stats_result_IvsIIA %>%
      group_by(Target) %>%
      summarise(
        mean_I = round(mean(y[fiber_type == "Type I"], na.rm = TRUE), 2),
        sd_I = round(sd(y[fiber_type == "Type I"], na.rm = TRUE), 2),
        n_I = sum(fiber_type == "Type I" & !is.na(y)),
        mean_IIA = round(mean(y[fiber_type == "Type IIA"], na.rm = TRUE), 2),
        sd_IIA = round(sd(y[fiber_type == "Type IIA"], na.rm = TRUE), 2),
        n_IIA = sum(fiber_type == "Type IIA" & !is.na(y)),
        logFoldChange = mean_IIA - mean_I,
        FoldChange = round(2^logFoldChange, 2),
        p_value = tryCatch(wilcox.test(y ~ fiber_type, data = pick(everything()))$p.value, error = function(e) NA),
        FDR = p.adjust(as.numeric(p_value), method = "bonferroni", n = nrow(gene_to_file)),
        .groups = 'drop'
      ) %>%
      mutate(
        summary_I = sprintf("%.1f ± %.1f, n = %d", mean_I, sd_I, n_I),
        summary_IIA = sprintf("%.1f ± %.1f, n = %d", mean_IIA, sd_IIA, n_IIA),
        Significance = case_when(
          FDR < 0.001 ~ "***",
          FDR < 0.01  ~ "**",
          FDR < 0.05  ~ "*",
          TRUE ~ "ns"
        ),
        p_value = format(p_value, scientific = TRUE, digits = 2),
        FDR = format(FDR, scientific = TRUE, digits = 2)
      ) %>%
      select(Target, logFoldChange, FoldChange, p_value, FDR, Significance, summary_I, summary_IIA) %>%
      mutate(across(everything(), ~ ifelse(is.na(.), "NA", .))) %>%
      as.data.frame() %>%
      tibble::column_to_rownames("Target")

    # Reformat table for display
    stats_result_IvsIIA <- data.frame(Statistics = colnames(stats_result_IvsIIA),
                                      t(stats_result_IvsIIA))

    # rename rows
    stats_result_IvsIIA$Statistics <- gsub("summary_IIA", "Type IIA (mean ± sd, n)", stats_result_IvsIIA$Statistics)
    stats_result_IvsIIA$Statistics <- gsub("summary_I", "Type I (mean ± sd, n)", stats_result_IvsIIA$Statistics)
    stats_result_IvsIIA$Statistics <- gsub("logFoldChange", "log2(fold-change)", stats_result_IvsIIA$Statistics)
    stats_result_IvsIIA$Statistics <- gsub("FoldChange", "Fold-change", stats_result_IvsIIA$Statistics)
    stats_result_IvsIIA$Statistics <- gsub("FDR", "FDR (Bonferroni)", stats_result_IvsIIA$Statistics)
    
    return(stats_result_IvsIIA)
  })
  
  
  #-----------------------------------------------------------------
  # Compute statistics (Wilcoxon test + summary)
  statisticsData_IvsIIX <- reactive({
    # Validate at least two fibers are selected
    validate(
      need(length(input$fibers) >= 2,
           "Please select at least two fiber types to calculate statistics.")
    )
    
    dat <- selectedTargetData()
    
    # Group and compute statistics per gene - Type I vs Type IIX
    stats_result_IvsIIX <- dat %>%
      filter(fiber_type %in% c("Type I", "Type IIX")) 

    # Validate that both fiber types are present
    validate(
      need(length(unique(stats_result_IvsIIX$fiber_type)) == 2,
           "Statistics impossible to calculate with the selected criteria: one of the fiber types (Type I or Type IIX) is missing.")
    )
    
    # Make table
    stats_result_IvsIIX <- stats_result_IvsIIX %>%
      group_by(Target) %>%
      summarise(
        mean_I = round(mean(y[fiber_type == "Type I"], na.rm = TRUE), 2),
        sd_I = round(sd(y[fiber_type == "Type I"], na.rm = TRUE), 2),
        n_I = sum(fiber_type == "Type I" & !is.na(y)),
        mean_IIX = round(mean(y[fiber_type == "Type IIX"], na.rm = TRUE), 2),
        sd_IIX = round(sd(y[fiber_type == "Type IIX"], na.rm = TRUE), 2),
        n_IIX = sum(fiber_type == "Type IIX" & !is.na(y)),
        logFoldChange = mean_IIX - mean_I,
        FoldChange = round(2^logFoldChange, 2),
        p_value = tryCatch(wilcox.test(y ~ fiber_type, data = pick(everything()))$p.value, error = function(e) NA),
        FDR = p.adjust(as.numeric(p_value), method = "bonferroni", n = nrow(gene_to_file)),
        .groups = 'drop'
      ) %>%
      mutate(
        summary_I = sprintf("%.1f ± %.1f, n = %d", mean_I, sd_I, n_I),
        summary_IIX = sprintf("%.1f ± %.1f, n = %d", mean_IIX, sd_IIX, n_IIX),
        Significance = case_when(
          FDR < 0.001 ~ "***",
          FDR < 0.01  ~ "**",
          FDR < 0.05  ~ "*",
          TRUE ~ "ns"
        ),
        p_value = format(p_value, scientific = TRUE, digits = 2),
        FDR = format(FDR, scientific = TRUE, digits = 2)
      ) %>%
      select(Target, logFoldChange, FoldChange, p_value, FDR, Significance, summary_I, summary_IIX) %>%
      mutate(across(everything(), ~ ifelse(is.na(.), "NA", .))) %>%
      as.data.frame() %>%
      tibble::column_to_rownames("Target")
    
    # Reformat table for display
    stats_result_IvsIIX <- data.frame(Statistics = colnames(stats_result_IvsIIX),
                                      t(stats_result_IvsIIX))
    
    # rename rows
    stats_result_IvsIIX$Statistics <- gsub("summary_IIX", "Type IIX (mean ± sd, n)", stats_result_IvsIIX$Statistics)
    stats_result_IvsIIX$Statistics <- gsub("summary_I", "Type I (mean ± sd, n)", stats_result_IvsIIX$Statistics)
    stats_result_IvsIIX$Statistics <- gsub("logFoldChange", "log2(fold-change)", stats_result_IvsIIX$Statistics)
    stats_result_IvsIIX$Statistics <- gsub("FoldChange", "Fold-change", stats_result_IvsIIX$Statistics)
    stats_result_IvsIIX$Statistics <- gsub("FDR", "FDR (Bonferroni)", stats_result_IvsIIX$Statistics)
    
    return(stats_result_IvsIIX)
  })
  
  
  #-----------------------------------------------------------------
  # Compute statistics (Wilcoxon test + summary)
  statisticsData_IIAvsIIX <- reactive({
    # Validate at least two fibers are selected
    validate(
      need(length(input$fibers) >= 2,
           "Please select at least two fiber types to calculate statistics.")
    )
    
    dat <- selectedTargetData()
    
    # Group and compute statistics per gene - Type IIA vs Type IIX
    stats_result_IIAvsIIX <- dat %>%
      filter(fiber_type %in% c("Type IIA", "Type IIX")) 

    # Validate that both fiber types are present
    validate(
      need(length(unique(stats_result_IIAvsIIX$fiber_type)) == 2,
           "Statistics impossible to calculate with the selected criteria: one of the fiber types is missing.")
    )
    
    # Make table
    stats_result_IIAvsIIX <- stats_result_IIAvsIIX %>%
      group_by(Target) %>%
      summarise(
        mean_IIA = round(mean(y[fiber_type == "Type IIA"], na.rm = TRUE), 2),
        sd_IIA = round(sd(y[fiber_type == "Type IIA"], na.rm = TRUE), 2),
        n_IIA = sum(fiber_type == "Type IIA" & !is.na(y)),
        mean_IIX = round(mean(y[fiber_type == "Type IIX"], na.rm = TRUE), 2),
        sd_IIX = round(sd(y[fiber_type == "Type IIX"], na.rm = TRUE), 2),
        n_IIX = sum(fiber_type == "Type IIX" & !is.na(y)),
        logFoldChange = mean_IIX - mean_IIA,
        FoldChange = round(2^logFoldChange, 2),
        p_value = tryCatch(wilcox.test(y ~ fiber_type, data = cur_data())$p.value, error = function(e) NA),
        FDR = p.adjust(as.numeric(p_value), method = "bonferroni", n = nrow(gene_to_file)),
        .groups = 'drop'
      ) %>%
      mutate(
        summary_IIA = sprintf("%.1f ± %.1f, n = %d", mean_IIA, sd_IIA, n_IIA),
        summary_IIX = sprintf("%.1f ± %.1f, n = %d", mean_IIX, sd_IIX, n_IIX),
        Significance = case_when(
          FDR < 0.001 ~ "***",
          FDR < 0.01  ~ "**",
          FDR < 0.05  ~ "*",
          TRUE ~ "ns"
        ),
        p_value = format(p_value, scientific = TRUE, digits = 2),
        FDR = format(FDR, scientific = TRUE, digits = 2)
      ) %>%
      select(Target, logFoldChange, FoldChange, p_value, FDR, Significance, summary_IIA, summary_IIX) %>%
      mutate(across(everything(), ~ ifelse(is.na(.), "NA", .))) %>%
      as.data.frame() %>%
      tibble::column_to_rownames("Target")
    
    # Reformat table for display
    stats_result_IIAvsIIX <- data.frame(Statistics = colnames(stats_result_IIAvsIIX),
                                        t(stats_result_IIAvsIIX))
    
    # rename rows
    stats_result_IIAvsIIX$Statistics <- gsub("summary_IIX", "Type IIX (mean ± sd, n)", stats_result_IIAvsIIX$Statistics)
    stats_result_IIAvsIIX$Statistics <- gsub("summary_IIA", "Type IIA (mean ± sd, n)", stats_result_IIAvsIIX$Statistics)
    stats_result_IIAvsIIX$Statistics <- gsub("logFoldChange", "log2(fold-change)", stats_result_IIAvsIIX$Statistics)
    stats_result_IIAvsIIX$Statistics <- gsub("FoldChange", "Fold-change", stats_result_IIAvsIIX$Statistics)
    stats_result_IIAvsIIX$Statistics <- gsub("FDR", "FDR (Bonferroni)", stats_result_IIAvsIIX$Statistics)
    
    return(stats_result_IIAvsIIX)
  })
  
  
  # Display significance table only
  output$statisticsTable <- DT::renderDT({
    req(input$inputComparison)
    
    # Dynamically select the correct data
    stats_df <- switch(input$inputComparison,
                       "Type I vs Type IIA" = statisticsData_IvsIIA(),
                       "Type I vs Type IIX" = statisticsData_IvsIIX(),
                       "Type IIA vs Type IIX" = statisticsData_IIAvsIIX())
    
    # If no valid dataset, don't render the table
    req(!is.null(stats_df))
    
    # Rename first column to match selected comparison
    colnames(stats_df)[1] <- input$inputComparison
    
    DT::datatable(
      stats_df,
      escape = FALSE, 
      rownames = FALSE,
      options = list(
        searching = FALSE,
        paging = FALSE,
        info = FALSE,
        ordering = FALSE,
        dom = 't',
        columnDefs = list(
          list(targets = 0, width = '30rem'),  # Set fixed width
          list(targets = 1:(ncol(stats_df)-1), className = 'dt-center')
        )
      )
    )
  })

  
  #-----------------------------------------------------------------
  # Show references table (metadata about datasets)
  output$references <- DT::renderDT({
    DT::datatable(
      references, 
      escape = FALSE, 
      rownames = FALSE,
      options = list(
        columnDefs = list(
          list(className = 'dt-center', targets = 1),
          list(className = 'dt-center', targets = 2)
        ),
        searching = FALSE,
        paging = FALSE,
        info = FALSE,
        dom = 't'
      )
    )
  })
  
  #-----------------------------------------------------------------
  # Download: Boxplot as PNG
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("FiberTypes_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      plot_obj <- boxplotTarget() + labs(caption = "Plot generated on MuscleOmics.org")
      png(filename = file, width = 3200, height = 1800, res = 300)
      print(plot_obj)
      dev.off()
    }
  )
  
  #-----------------------------------------------------------------
  # Download: Raw expression data (selected gene(s))
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("FiberTypes_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- selectedTargetData()
      df <- pivot_wider(df, names_from = "Target", values_from = "y")
      if (!is.null(df)) {
        write.table(df, file, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)
        cat("\n# Data generated on MuscleOmics.org\n", file = file, append = TRUE)
      }
    }
  )
  
  #-----------------------------------------------------------------
  # Download: Statistics table (all comparisons combined)
  output$downloadStats <- downloadHandler(
    filename = function() {
      paste0("FiberTypes_statistics_", Sys.Date(), ".csv")
    },
    content = function(file) {
      # Get the three data frames from their reactive expressions
      df_IvsIIA <- statisticsData_IvsIIA()
      df_IvsIIX <- statisticsData_IvsIIX()
      df_IIAvsIIX <- statisticsData_IIAvsIIX()
      
      # Add a 'Comparison' column to each one
      df_IvsIIA$Comparison <- "Type I vs Type IIA"
      df_IvsIIX$Comparison <- "Type I vs Type IIX"
      df_IIAvsIIX$Comparison <- "Type IIA vs Type IIX"
      
      # Reorder columns to put 'Comparison' first
      df_IvsIIA <- df_IvsIIA[, c(ncol(df_IvsIIA), 1:(ncol(df_IvsIIA) - 1))]
      df_IvsIIX <- df_IvsIIX[, c(ncol(df_IvsIIX), 1:(ncol(df_IvsIIX) - 1))]
      df_IIAvsIIX <- df_IIAvsIIX[, c(ncol(df_IIAvsIIX), 1:(ncol(df_IIAvsIIX) - 1))]
      
      # Combine all into a single data frame
      combined_df <- rbind(df_IvsIIA, df_IvsIIX, df_IIAvsIIX)
      
      # Write to CSV
      write.table(combined_df, file, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)
      
      # Append generation note
      cat("\n# Data generated on MuscleOmics.org\n", file = file, append = TRUE)
    }
  )
  
}
