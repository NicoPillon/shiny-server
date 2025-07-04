#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic myotube response to electrical pulse stimulation
#
#--------------------------------------------------------------------------------------------------------

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
      updateSelectizeInput(session, "inputGeneSymbol",
                           choices = gene_list$TARGET,
                           selected = selected_target,
                           server = TRUE)
    } else {
      updateSelectizeInput(session, "inputGeneSymbol",
                           choices = gene_list$TARGET,
                           selected = c("HMOX1", "NFKBIZ", "NR4A3", "PLEKHA4"),  # default genes
                           server = TRUE)
    }
  })
  
  #-----------------------------------------------
  # Reset button functionality: resets all UI inputs to their default values
  observeEvent(input$resetInputs, {
    updateSelectizeInput(session, "inputGeneSymbol", selected = character(0))
    updateCheckboxGroupInput(session, "duration", selected = c("3", "3+Rest3h", "6", "24", "48", "72", "72(1h+3rest)"))
    updateCheckboxGroupInput(session, "pulse_param", selected = c("14V, 5Hz, 2ms", "3V, 1Hz, 6ms", "13V, 66Hz, 2ms",
                                                               "10V, 1Hz, 2ms", "30V, 1Hz, 2ms", "10V, 1Hz, 10ms", 
                                                               "11.5V, 1Hz, 2ms"))
    updateCheckboxGroupInput(session, "cell_type", selected = c("mouse C2C12", "human primary"))
  })
  
  #-----------------------------------------------------------------
  # REACTIVE: Load only selected gene(s) from on-disk parquet files
  selectedGeneData <- reactive({
    req(input$inputGeneSymbol)  # Ensure input is not empty
    genename <- toupper(input$inputGeneSymbol)  # Convert input to uppercase for consistency
    
    # Find matching file(s) and row(s) for each gene
    file_row <- gene_list[gene_list$TARGET %in% genename, ]
    split_rows <- split(file_row, file_row$file)  # Group by file for efficient loading
    
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
    missing_genes <- setdiff(genename, rownames(selected_row))
    if (length(missing_genes) > 0) {
      missing_df <- matrix(NA, nrow = length(missing_genes), ncol = ncol(selected_row))
      rownames(missing_df) <- missing_genes
      colnames(missing_df) <- colnames(selected_row)
      selected_row <- rbind(selected_row, missing_df)
    }
    
    # Ensure output preserves original input order
    df <- selected_row[genename, , drop = FALSE]
    data.frame(df)
    
    # Combine metadata with expression values (transpose to long format)
    dat <- data.frame(metadata, t(df))
    dat <- pivot_longer(dat, cols = c(15:ncol(dat)),
                        values_to = "y",
                        names_to = "Gene")
    
    # Identify pairIDs where POST Time is within the selected time range
    valid_pairs <- dat %>%
      filter(Condition == "EPS",
             time_hours %in% input$duration) %>%
      pull(pairID)
    
    # Filter the full dataset by valid pairIDs and other filters
    dat <- dat %>%
      filter(cell.type %in% input$cell_type,
             Pulse_parameters %in% input$pulse_param,
             pairID %in% valid_pairs)
  })
  
  
  #-----------------------------------------------------------------
  # Generate Boxplot visualization of gene expression
  plotDataBox <- reactive({
    dat <- selectedGeneData()
    
    # Validate that some data is available
    validate(
      need(nrow(dat) > 0, "No data available for the selected filters. Please adjust your selections.")
    )
    
    # Generate boxplot with ggplot2
    ggplot(dat, aes(x=Gene, y=y, fill=Condition)) + 
      geom_boxplot(outlier.size = 0.1, alpha = 0.5, position = position_dodge(0.8))  + 
      geom_sina(size = 1.5, position = position_dodge(0.8), alpha = 0.1) +  # Sina plot for distribution
      theme_minimal(17, base_family = "Arial") + 
      theme(axis.text.x = element_text(face = "bold", size = 14, color = "black"),
            axis.text.y = element_text("black"),
            axis.title.y = element_text(face = "bold", color = "black")) +
      labs(x="", y="mRNA expression (log2)", fill = "Treatment") +
      scale_shape_manual(values=rep(c(15,16,17), 20)) +
      scale_fill_manual(values = c("#5B768E", "#bd1a0e")) +  # Custom treatment colors
      scale_y_continuous(expand = expansion(mult = c(0.05, .15)))
  })
  
  # Render the boxplot
  output$geneBoxplot <- renderPlot({
    plotDataBox()
  })
  
  
  #-----------------------------------------------------------------
  # Reactive text description of the current filters (used in UI summary)
  filterSummary <- reactive({
    
    # Load gene expression data and merge with metadata
    dat <- selectedGeneData()

    # If no rows remain, return fallback
    if (nrow(dat) == 0) {
      return(" ")
    }
    
    # Format a vector as a human-readable list with "or" before the last item
    format_list1 <- function(x) {
      x <- sort(unique(x))  # ensure unique + sorted
      n <- length(x)
      if (n == 0) return("")
      if (n == 1) return(x)
      if (n == 2) return(paste(x, collapse = " and "))
      paste(paste(x[-n], collapse = ", "), "and", x[n])
    }
    
    format_list2 <- function(x) {
      x <- sort(unique(x))  # ensure unique + sorted
      n <- length(x)
      if (n == 0) return("")
      if (n == 1) return(x)
      if (n == 2) return(paste(x, collapse = " or "))
      paste(paste(x[-n], collapse = ", "), "or", x[n])
    }
    
    # Generate a combined cell type description (e.g., "human C2C12")
    cell_types <- format_list1(dat$cell.type)
    
    # Unique concentrations and durations (sorted for clarity)
    pulse_param <- format_list2(dat$Pulse_parameters)
    time <- format_list2(dat$time_hours)
    
    # Construct final sentence
    glue::glue(
      "Based on your selection, the plot and statistics reflect data from {cell_types} myotubes, ",
      "exposed to {pulse_param} EPS for {time} hours."
    )
  })
  
  output$filterSummaryText <- renderText({
    filterSummary()
  })
  
  
  #-----------------------------------------------------------------
  # Compute statistics (Wilcoxon test + summary)
  statisticsData <- reactive({
    dat <- selectedGeneData()

    # Ensure there's data to analyze
    validate(
      need(nrow(dat) > 0, "No data available for the selected filters. Please adjust your selections.")
    )
    
    # Group and compute statistics per gene
    stats_result <- dat %>%
      group_by(Gene) %>%
      summarise(
        mean_control = round(mean(y[Condition == "Control"], na.rm = TRUE), 2),
        sd_control = round(sd(y[Condition == "Control"], na.rm = TRUE), 2),
        n_control = sum(Condition == "Control" & !is.na(y)),
        mean_EPS = round(mean(y[Condition == "EPS"], na.rm = TRUE), 2),
        sd_EPS = round(sd(y[Condition == "EPS"], na.rm = TRUE), 2),
        n_EPS = sum(Condition == "EPS" & !is.na(y)),
        logFoldChange = mean_EPS - mean_control,
        FoldChange = round(2^logFoldChange, 2),
        p_value = tryCatch(wilcox.test(y ~ Condition, data = pick(everything()), exact = FALSE)$p.value, error = function(e) NA),
        FDR = p.adjust(as.numeric(p_value), method = "bonferroni", n = nrow(gene_list)),
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
      as.data.frame() %>%
      tibble::column_to_rownames("Gene")
    
    # Reformat table for display
    stats_result <- data.frame(Statistics = colnames(stats_result),
                               t(stats_result))
    stats_result$Statistics <- gsub("_", " ", stats_result$Statistics)
    stats_result$Statistics <- gsub("logFoldChange", "log2(fold-change)", stats_result$Statistics)
    stats_result$Statistics <- gsub("FoldChange", "Fold-change", stats_result$Statistics)
    stats_result$Statistics <- gsub("FDR", "FDR (Bonferroni)", stats_result$Statistics)
    
    return(stats_result)
  })
  
  # Display significance table only
  output$statistics1 <- DT::renderDT({
    dat <- statisticsData()
    dat <- dat[!dat$Statistics %in% c("mean control", "sd control", "n control", "mean EPS", "sd EPS", "n EPS"),]
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
    dat <- dat[dat$Statistics %in% c("mean control", "sd control", "n control", "mean EPS", "sd EPS", "n EPS"),]
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
      paste0("MyotubeEPS_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      plot_obj <- plotDataBox() + labs(caption = "Plot generated on MuscleOmics.org")
      png(filename = file, width = 3200, height = 1800, res = 300)
      print(plot_obj)
      dev.off()
    }
  )
  
  #-----------------------------------------------------------------
  # Download: Raw expression data (selected gene(s))
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("MyotubeEPS_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- selectedGeneData()
      df <- pivot_wider(df, names_from = "Gene", values_from = "y")
      df <- data.frame(metadata, t(df))
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
      paste0("MyotubeEPS_statistics_", Sys.Date(), ".csv")
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

  
  
  
  