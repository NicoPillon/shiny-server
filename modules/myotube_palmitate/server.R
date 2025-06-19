#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle cells - Server
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
  # Reset button functionality: resets all UI inputs to their default values
  observeEvent(input$resetInputs, {
    updateSelectizeInput(session, "inputGeneSymbol", selected = character(0))
    updateSliderInput(session, "concentration", value = c(100, 500))
    updateSliderInput(session, "duration", value = c(12, 96))
    updateCheckboxGroupInput(session, "cell_type", selected = c("human primary", "human LHCN-M2", "mouse C2C12", "rat primary"))
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
    selected_row <- selected_row[genename, , drop = FALSE]
    data.frame(selected_row)
  })
  
  #-----------------------------------------------------------------
  # Generate Boxplot visualization of gene expression
  plotDataBox <- reactive({
    df <- selectedGeneData()
    
    # Combine metadata with expression values (transpose to long format)
    dat <- data.frame(metadata, t(df))
    dat <- pivot_longer(dat, cols = c(12:ncol(dat)),
                        values_to = "y",
                        names_to = "Gene")
    
    # Filter dataset based on selected slider/checkbox inputs
    dat <- dplyr::filter(dat,
                         concentration.micromolar >= input$concentration[1] & concentration.micromolar <= input$concentration[2],
                         time.hours >= input$duration[1] & time.hours <= input$duration[2],
                         cell.type %in% input$cell_type)
    
    # Validate that some data is available
    validate(
      need(nrow(dat) > 0, "No data available for the selected filters. Please adjust your selections.")
    )
    
    # Generate boxplot with ggplot2
    ggplot(dat, aes(x=Gene, y=y, fill=treatment)) + 
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
    req(input$cell_type, input$concentration, input$duration)
    
    # Load gene expression data and merge with metadata
    df <- selectedGeneData()
    dat <- data.frame(metadata, t(df))
    
    # Reshape to long format for proper filtering
    dat <- pivot_longer(dat, cols = c(12:ncol(dat)),
                        values_to = "data",
                        names_to = "Gene")
    
    # Apply current filters
    dat <- dplyr::filter(dat,
                         concentration.micromolar >= input$concentration[1] & concentration.micromolar <= input$concentration[2],
                         time.hours >= input$duration[1] & time.hours <= input$duration[2],
                         cell.type %in% input$cell_type)
    
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
    conc <- format_list2(dat$concentration.micromolar)
    time <- format_list2(dat$time.hours)
    
    # Construct final sentence
    glue::glue(
      "Based on your selection, the plot and statistics reflect data from {cell_types} myotubes, ",
      "exposed to {conc} µmol/L palmitate for {time} hours."
    )
  })
  
  output$filterSummaryText <- renderText({
    filterSummary()
  })
  
  #-----------------------------------------------------------------
  # Compute statistics (Wilcoxon test + summary)
  statisticsData <- reactive({
    df <- selectedGeneData()
    dat <- data.frame(metadata, t(df))  # Merge with metadata
    
    # Convert to long format for stats
    dat <- pivot_longer(dat, cols = c(12:ncol(dat)),
                        values_to = "data",
                        names_to = "Gene")
    
    # Apply same filters as for plot
    dat <- dplyr::filter(dat,
                         concentration.micromolar >= input$concentration[1] & concentration.micromolar <= input$concentration[2],
                         time.hours >= input$duration[1] & time.hours <= input$duration[2],
                         cell.type %in% input$cell_type)
    
    # Ensure there's data to analyze
    validate(
      need(nrow(dat) > 0, "No data available for the selected filters. Please adjust your selections.")
    )
    
    # Group and compute statistics per gene
    stats_result <- dat %>%
      group_by(Gene) %>%
      summarise(
        mean_control = round(mean(data[treatment == "control"], na.rm = TRUE), 2),
        sd_control = round(sd(data[treatment == "control"], na.rm = TRUE), 2),
        n_control = sum(treatment == "control" & !is.na(data)),
        mean_palmitate = round(mean(data[treatment == "palmitate"], na.rm = TRUE), 2),
        sd_palmitate = round(sd(data[treatment == "palmitate"], na.rm = TRUE), 2),
        n_palmitate = sum(treatment == "palmitate" & !is.na(data)),
        logFoldChange = mean_palmitate - mean_control,
        FoldChange = round(2^logFoldChange, 2),
        p_value = tryCatch(wilcox.test(data ~ treatment, data = cur_data())$p.value, error = function(e) NA),
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
  output$statistics1 <- DT::renderDataTable({
    dat <- statisticsData()
    dat <- dat[!dat$Statistics %in% c("mean control", "sd control", "n control", "mean palmitate", "sd palmitate", "n palmitate"),]
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
  output$statistics2 <- DT::renderDataTable({
    dat <- statisticsData()
    dat <- dat[dat$Statistics %in% c("mean control", "sd control", "n control", "mean palmitate", "sd palmitate", "n palmitate"),]
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
  output$references <- renderDataTable({
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
      paste0("MyotubePalmitate_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- selectedGeneData()
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
