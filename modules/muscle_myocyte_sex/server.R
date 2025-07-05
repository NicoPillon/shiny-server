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
                           selected = c("RPS4Y1"),  # default genes
                           server = TRUE)
    }
  })
  
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
    df <- selected_row[genename, , drop = FALSE]
    data.frame(df)
    
    # Merge with metadata
    dat <- data.frame(metadata, t(df)) 
    
    # Convert to long format for stats
    dat <- pivot_longer(dat, cols = c(17:ncol(dat)),
                        values_to = "y",
                        names_to = "Gene")
    
  })
  
  #-----------------------------------------------------------------
  # Generate Boxplot visualization of gene expression
  plotDataBox <- reactive({
    dat <- selectedGeneData()
    
    # Validate that some data is available
    shiny::validate(
      need(nrow(dat) > 0, "No data available for the selected filters. Please adjust your selections.")
    )
    
    # Generate boxplot with ggplot2
    ggplot(dat, aes(x=cellType, y=y, fill=sex)) + 
      geom_boxplot(outlier.size = 0.1, alpha = 0.5, position = position_dodge(0.8))  + 
      geom_sina(size = 1.5, position = position_dodge(0.8), alpha = 0.1) +  # Sina plot for distribution
      facet_wrap(~ Gene, ncol =3)+
      theme_minimal(17, base_family = "Arial") + 
      theme(axis.text.x = element_text(face = "bold", size = 14, color = "black"),
            axis.text.y = element_text("black"),
            axis.title.y = element_text(face = "bold", color = "black")) +
      labs(x="", y="mRNA expression (log2)", fill = "Sex") +
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
    cell_types <- format_list1(dat$cellType)
    
  })
  
  output$filterSummaryText <- renderText({
    filterSummary()
  })
  

  #-----------------------------------------------------------------
  # Compute statistics (Wilcoxon test + summary)
  statisticsData <- reactive({
    dat <- selectedGeneData()
    
    # Ensure there's y to analyze
    shiny::validate(
      need(nrow(dat) > 0, "No data available for the selected filters. Please adjust your selections.")
    )
    
    # Group and compute statistics per gene
    stats_result <- dat %>%
      group_by(Gene) %>%
      summarise(
        ## Skeletal muscle
        Female_skeletal_muscle = paste(round(mean(y[sex == "Female" & cellType == "Skeletal muscle"], na.rm = TRUE), 2),
                                            " ± ",
                                                          round(sd(y[sex == "Female" & cellType == "Skeletal muscle"], na.rm = TRUE), 2),
                                            " (n:",
                                            sum(sex == "Female" & !is.na(y) & cellType == "Skeletal muscle"),
                                            ")",
          sep =""),
        Male_skeletal_muscle = paste(round(mean(y[sex == "Male" & cellType == "Skeletal muscle"], na.rm = TRUE), 2),
                                                          " ± ",
                                                          round(sd(y[sex == "Male" & cellType == "Skeletal muscle"], na.rm = TRUE), 2),
                                                          " (n:",
                                                          sum(sex == "Male" & !is.na(y) & cellType == "Skeletal muscle"),
                                                          ")",
                                                          sep =""),
        
        logFoldChange_skeletal_muscle = round(mean(y[sex == "Female" & cellType == "Skeletal muscle"], na.rm = TRUE) -
                                                mean(y[sex == "Male"  & cellType == "Skeletal muscle"], na.rm = TRUE), 2),
        FoldChange_skeletal_muscle = round(2^logFoldChange_skeletal_muscle, 2),
        p_value_skeletal_muscle = tryCatch(wilcox.test(y ~ sex, data = pick(everything()), exact = FALSE)$p.value, error = function(e) NA),
        FDR_skeletal_muscle = p.adjust(as.numeric(p_value_skeletal_muscle), method = "bonferroni", n = nrow(gene_list)),
       ## Myocyte
       Female_skeletal_myocyte = paste(round(mean(y[sex == "Female" & cellType == "Skeletal myocyte"], na.rm = TRUE), 2),
                                                         " ± ",
                                                         round(sd(y[sex == "Female" & cellType == "Skeletal myocyte"], na.rm = TRUE), 2),
                                                         " (n:",
                                                         sum(sex == "Female" & !is.na(y) & cellType == "Skeletal myocyte"),
                                                         ")",
                                                         sep =""),
       Male_skeletal_myocyte = paste(round(mean(y[sex == "Male" & cellType == "Skeletal myocyte"], na.rm = TRUE), 2),
                                                       " ± ",
                                                       round(sd(y[sex == "Male" & cellType == "Skeletal myocyte"], na.rm = TRUE), 2),
                                                       " (n:",
                                                       sum(sex == "Male" & !is.na(y) & cellType == "Skeletal myocyte"),
                                                       ")",
                                                       sep =""),
       
       logFoldChange_skeletal_myocyte = round(mean(y[sex == "Female" & cellType == "Skeletal myocyte"], na.rm = TRUE) -
                                               mean(y[sex == "Male"  & cellType == "Skeletal myocyte"], na.rm = TRUE), 2),
        FoldChange_skeletal_myocyte = round(2^logFoldChange_skeletal_myocyte, 2),
        p_value_skeletal_myocyte = tryCatch(wilcox.test(y ~ sex, data = pick(everything()), exact = FALSE)$p.value, error = function(e) NA),
        FDR_skeletal_myocyte = p.adjust(as.numeric(p_value_skeletal_myocyte), method = "bonferroni", n = nrow(gene_list)),
        .groups = 'drop'
      ) %>%
      mutate(
        Significance_skeletal_muscle = case_when(
          FDR_skeletal_muscle < 0.001 ~ "***",
          FDR_skeletal_muscle < 0.01  ~ "**",
          FDR_skeletal_muscle < 0.05  ~ "*",
          TRUE ~ "ns"
        ),
        Significance_skeletal_myocyte = case_when(
          FDR_skeletal_myocyte < 0.001 ~ "***",
          FDR_skeletal_myocyte < 0.01  ~ "**",
          FDR_skeletal_myocyte < 0.05  ~ "*",
          TRUE ~ "ns"
        ),
        p_value_skeletal_muscle = format(p_value_skeletal_muscle, scientific = TRUE, digits = 2),
        FDR_skeletal_muscle= format(FDR_skeletal_muscle, scientific = TRUE, digits = 2),
        p_value_skeletal_myocyte = format(p_value_skeletal_myocyte, scientific = TRUE, digits = 2),
        FDR_skeletal_myocyte = format(FDR_skeletal_myocyte, scientific = TRUE, digits = 2)
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
    dat <- dat[dat$Statistics %in% c("log2(fold-change) skeletal muscle",
                                      "Fold-change skeletal muscle",
                                      "p value skeletal muscle",
                                      "FDR (Bonferroni) skeletal muscle",
                                     "Significance skeletal muscle"
                                      ),]
    colnames(dat)[1] <- "Differential Expression Analysis: Skeletal muscle"
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
  
  # Display significance table only
  output$statistics2 <- DT::renderDT({
    dat <- statisticsData()
    dat <- dat[dat$Statistics %in% c("log2(fold-change) skeletal myocyte",
                                     "Fold-change skeletal myocyte",
                                     "p value skeletal myocyte",
                                     "FDR (Bonferroni) skeletal myocyte",
                                     "Significance skeletal myocyte"
    ),]
    colnames(dat)[1] <- "Differential Expression Analysis: Skeletal myocyte"
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
  output$statistics3 <- DT::renderDT({
    dat <- statisticsData()
    dat <- dat[dat$Statistics %in% c("Female skeletal muscle",
                                     "Male skeletal muscle",
                                     "Female skeletal myocyte",
                                     "Male skeletal myocyte"
    ),]
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
      df <- pivot_wider(df, names_from = "Gene", values_from = "y")
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
