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
  
  #-----------------------------------------------------------------
  # REACTIVE: load only selected gene(s) from dataset
  selectedTargetData <- reactive({
    req(input$inputTarget, input$inputOmics)
    # targetname <- c("LDHA", "LDHB", "MYH7", "MYH1", "MYH2")
    # omics <- "Transcriptome"
    
    targetname <- toupper(unlist(strsplit(input$inputTarget, "[,;\\s]+")))
    omics <- input$inputOmics
    
    # Match gene to file and row info - gene and OMICS
    file_row <- gene_to_file %>%
      filter(TARGET %in% targetname & OMICS %in% omics)
    
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
    
    dat <- selectedTargetData()

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
      return(" ")
    }
    
    # Format readable list of fiber types
    format_list <- function(x) {
      x <- sort(unique(x))
      n <- length(x)
      if (n == 0) return("")
      if (n == 1) return(x)
      if (n == 2) return(paste(x, collapse = " and "))
      paste(paste(x[-n], collapse = ", "), "and", x[n])
    }
    
    fiber_types <- sort(unique(dat$fiber_type))
    fiber_text <- format_list(fiber_types)
    
    # Flags
    has_I   <- "Type I"   %in% fiber_types
    has_IIA <- "Type IIA" %in% fiber_types
    has_IIX <- "Type IIX" %in% fiber_types
    
    # Custom phrasing for full selection
    if (all(c("Type I", "Type IIA", "Type IIX") %in% fiber_types)) {
      comparison_sentence <- "Statistics are calculated comparing Type I to Type IIA/X or Type IIX to Type IIA fibers."
    } else {
      comparison_text <- c()
      if (has_I & has_IIA) comparison_text <- c(comparison_text, "Type I to Type IIA")
      if (has_I & has_IIX) comparison_text <- c(comparison_text, "Type I to Type IIX")
      if (has_IIA & has_IIX) comparison_text <- c(comparison_text, "Type IIX to Type IIA")
      
      if (length(comparison_text) == 0) {
        comparison_sentence <- "No valid pairwise comparisons available based on selected fiber types."
      } else {
        comparison_sentence <- paste(
          "Statistics are calculated comparing",
          paste(comparison_text, collapse = " or "),
          "fibers."
        )
      }
    }
    
    # Final output
    glue::glue(
      "Based on your selection, the plot and statistics reflect data from {fiber_text} muscle fibers. {comparison_sentence}"
    )
  })
  
  
  
  output$filterSummaryText <- renderText({
    filterSummary()
  })
  
  
  #-----------------------------------------------------------------
  # Compute statistics (Wilcoxon test + summary)
  statisticsData <- reactive({
    dat <- selectedTargetData()
    
    # Ensure there's y to analyze
    shiny::validate(
      need(nrow(dat) > 0, "No data available for the selected filters. Please adjust your selections.")
    )
    
    # Group and compute statistics per gene - Type I vs Type II
    stats_result_IvsII <- dat %>%
      mutate(fiber_type = ifelse(fiber_type %in% c("Type IIA", "Type IIX"), "Type II", fiber_type)) %>% 
      filter(fiber_type %in% c("Type I", "Type II")) %>%
      group_by(Target) %>%
      summarise(
        mean_I = round(mean(y[fiber_type == "Type I"], na.rm = TRUE), 2),
        sd_I = round(sd(y[fiber_type == "Type I"], na.rm = TRUE), 2),
        n_I = sum(fiber_type == "Type I" & !is.na(y)),
        mean_IIAX = round(mean(y[fiber_type == "Type II"], na.rm = TRUE), 2),
        sd_IIAX = round(sd(y[fiber_type == "Type II"], na.rm = TRUE), 2),
        n_IIAX = sum(fiber_type == "Type II" & !is.na(y)),
        logFoldChange = mean_IIAX - mean_I,
        FoldChange = round(2^logFoldChange, 2),
        p_value = tryCatch(wilcox.test(y ~ fiber_type, data = cur_data())$p.value, error = function(e) NA),
        FDR = p.adjust(as.numeric(p_value), method = "bonferroni", n = nrow(gene_to_file)),
        .groups = 'drop'
      ) %>%
      mutate(
        Significance = case_when(
          FDR < 0.001 ~ "***",
          FDR < 0.01  ~ "**",
          FDR < 0.05  ~ "*",
          TRUE ~ "ns"
        ),
        p_value = format(p_value, scientific = TRUE, digits = 2),
        FDR = format(FDR, scientific = TRUE, digits = 2)
      ) %>%
      mutate(across(everything(), ~ ifelse(is.na(.), "NA", .))) %>%
      as.data.frame() %>%
      tibble::column_to_rownames("Target")

    # Reformat table for display
    stats_result_IvsII <- data.frame(Statistics = colnames(stats_result_IvsII),
                                     t(stats_result_IvsII))
    stats_result_IvsII$Comparison <- "Type IIA/X vs Type I"
    
    # Group and compute statistics per gene - Type IIX vs Type IIA
    stats_result_IIXvsIIA <- dat %>%
      filter(fiber_type %in% c("Type IIA", "Type IIX")) %>%
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
        Significance = case_when(
          FDR < 0.001 ~ "***",
          FDR < 0.01  ~ "**",
          FDR < 0.05  ~ "*",
          TRUE ~ "ns"
        ),
        p_value = format(p_value, scientific = TRUE, digits = 2),
        FDR = format(FDR, scientific = TRUE, digits = 2)
      ) %>%
      mutate(across(everything(), ~ ifelse(is.na(.), "NA", .))) %>%
      as.data.frame() %>%
      tibble::column_to_rownames("Target")

    # Reformat table for display
    stats_result_IIXvsIIA <- data.frame(Statistics = colnames(stats_result_IIXvsIIA),
                                        t(stats_result_IIXvsIIA))
    stats_result_IIXvsIIA$Comparison <- "Type IIX vs Type IIA"
    
    # Merge stats
    stats_result <- rbind(stats_result_IvsII, stats_result_IIXvsIIA)
    
    # Remove rows with only NA or only NaN
    stats_result <- stats_result[
      !apply(stats_result[, input$inputTarget, drop = FALSE], 1, function(x) {
        all(trimws(as.character(x)) %in% c("NaN", "NA", "ns", "0", ""))
      }),
    ]
    
    # rename rows
    stats_result$Statistics <- gsub("_", " Type ", stats_result$Statistics)
    stats_result$Statistics <- gsub("AX", " A/X ", stats_result$Statistics)
    stats_result$Statistics <- gsub("logFoldChange", "log2(fold-change)", stats_result$Statistics)
    stats_result$Statistics <- gsub("FoldChange", "Fold-change", stats_result$Statistics)
    stats_result$Statistics <- gsub("FDR", "FDR (Bonferroni)", stats_result$Statistics)
    
    return(stats_result)
  })
  
  # Display significance table only
  output$statistics1 <- DT::renderDT({
    dat <- statisticsData()
    dat <- dat[!grepl("(mean |sd |n )", dat$Statistics), ]
    dat$Statistics <- paste(dat$Comparison, dat$Statistics, sep = ", ")
    dat$Comparison <- NULL
    colnames(dat)[1] <- "Differential Expression Analysis"
    DT::datatable(
      dat,
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
          list(targets = 1:(ncol(dat)-1), className = 'dt-center')
        )
      )
    )
  })
  
  # Display group statistics table
  output$statistics2 <- DT::renderDT({
    dat <- statisticsData()
    dat <- dat[grepl("(mean |sd |n )", dat$Statistics), ]
    dat <- dat[!grepl("(A/X)", dat$Statistics), ]
    dat$Comparison <- NULL
    colnames(dat)[1] <- "Group Summary Statistics"
    DT::datatable(
      dat,
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
          list(targets = 1:(ncol(dat)-1), className = 'dt-center')
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
      paste0("MyotubePalmitate_plot_", Sys.Date(), ".png")
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
      paste0("MyotubePalmitate_data_", Sys.Date(), ".csv")
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
  # Download: Statistics table
  output$downloadStats <- downloadHandler(
    filename = function() {
      paste0("MyotubePalmitate_statistics_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- statisticsData()
      if (!is.null(df)) {
        write.table(df, file, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)
        cat("\n# Data generated on MuscleOmics.org\n", file = file, append = TRUE)
      }
    }
  )
}
