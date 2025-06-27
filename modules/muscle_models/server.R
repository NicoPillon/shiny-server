#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle cells
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
                           selected = c("SLC2A4", "PPARGC1A", "MYH7", "MYH1"),  # default genes
                           server = TRUE)
    }
  })
  
  #-----------------------------------------------
  # Reset button functionality: resets all UI inputs to their default values
  observeEvent(input$resetInputs, {
    updateSelectizeInput(session, "inputGeneSymbol", selected = character(0))
    updateCheckboxGroupInput(session, "species", selected = c("Human", "Mouse", "Rat"))
    updateCheckboxGroupInput(session, "cell_tissue", selected = c("Cell", "Tissue"))
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
    dat <- pivot_longer(dat, cols = c(6:ncol(dat)),
                        values_to = "y",
                        names_to = "Gene")
    
    # Apply same filters as for plot
    dat <- dplyr::filter(dat,
                         species %in% input$species,
                         cell_tissue %in% input$cell_tissue)
  })
  
  #-----------------------------------------------------------------
  # Boxplots
  plotDataBox <- reactive({
    dat <- selectedGeneData()
    
    # Validate that some data is available
    shiny::validate(
      need(nrow(dat) > 0, "No data available for the selected filters. Please adjust your selections.")
    )
    
    dat$cell_type <- factor(dat$cell_type,
                            levels = c("Human Primary", "Mouse C2C12", "Rat L6",
                                       "Human Tissue", "Mouse Tissue", "Rat Tissue"))
    
    ggplot(dat, aes(x = cell_type, y = y, fill = cell_type)) + 
      geom_boxplot(outlier.size = 0.1, alpha = 0.7, position = position_dodge(0.8))  + 
      geom_sina(size = 1.5, position = position_dodge(0.8), alpha = 0.1) +
      theme_bw(base_size = 17) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_rect(fill = "#5B768E", color = "black"),
            strip.text = element_text(face = "bold", color = "white")) +
      labs(x = element_blank(), 
           y = "mRNA expression (log2)",
           fill = element_blank()) +
      scale_y_continuous(breaks = round(seq(-4, 8, by = 2), 1)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      scale_fill_manual(values = c(
        "Human Primary" = "#a4281e",  # dark red
        "Human Tissue"  = "#f4a582",  # light red/orange
        "Mouse C2C12"   = "#364e63",  # dark slate blue
        "Mouse Tissue"  = "#b6d3e5",  # light blue
        "Rat L6"        = "#5b5f2b",  # dark olive
        "Rat Tissue"    = "#d9e79a"   # light sage
      )) +
      facet_wrap(.~Gene)
    
    # ggplot(dat, aes(x = Gene, y = y, fill = cell_type)) + 
    #   geom_boxplot(outlier.size = 0.1, alpha = 0.5, position = position_dodge(0.8))  + 
    #   geom_sina(size = 1.5, position = position_dodge(0.8), alpha = 0.1) +
    #   theme_minimal(base_size = 17) + 
    #   labs(x = element_blank(), 
    #        y = "mRNA expression (log2)",
    #        fill = element_blank()) +
    #   scale_y_continuous(breaks = round(seq(-4, 8, by = 2), 1)) +
    #   geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    #   scale_fill_manual(values = c("HumanCell" = "#D95F02", 
    #                                "HumanTissue" = "#FDB863", 
    #                                "MouseC2C12" = "#1B9E77", 
    #                                "MouseTissue" = "#A6D854", 
    #                                "RatL6" = "#7570B3", 
    #                                "RatTissue" = "#CAB2D6"))
    
  })
  
  output$geneBoxplot <- renderPlot({
    plotDataBox()
  })
  
  
  #-----------------------------------------------------------------
  # Heatmap
  plotDataHeatmap <- eventReactive(input$updatePlot, {

    #make matrix with found gene names
    matrixData <- selectedGeneData()
    
    # Add gene names as a column for reshaping
    matrix_long <- matrixData %>%
      rownames_to_column("gene") %>%
      pivot_longer(
        cols = -gene,
        names_to = "sample",
        values_to = "expression"
      )
    matrix_long <- matrix_long %>%
      left_join(samples_list, by = c("sample" = "sample"))
    gene_medians <- matrix_long %>%
      group_by(gene, cell_tissue) %>%
      summarise(median_expression = median(expression, na.rm = TRUE), .groups = "drop")
    gene_matrix <- gene_medians %>%
      pivot_wider(
        names_from = cell_tissue,
        values_from = median_expression
      ) %>%
      column_to_rownames("gene")  # Optional: make gene names the rownames
    return(gene_matrix)
  })
  
  
  output$geneHeatmap <- renderPlot({
    
    #df_mean <- gene_matrix
    df_mean <- plotDataHeatmap()
    
    my_sample_col <- samples_list[,c(2,5,4)]
    my_sample_col <- my_sample_col[!duplicated(my_sample_col$cell_tissue),]
    my_sample_col <- data.frame(my_sample_col, row.names = 1)
    
    my_colour = list(
      #model = setNames(unique(samples_names$model.color), unique(samples_names$model)),
      species = setNames(rev(unique(samples_list$species_color)), unique(samples_list$species))
    )
    
    plotHeatmap <- pheatmap(df_mean, scale = "row",
                            annotation_col = my_sample_col,
                            annotation_colors = my_colour,
                            display_numbers = F,
                            cluster_rows = T,
                            cluster_cols = T,
                            fontsize=11,
                            color = hcl.colors(50, "BluYl", rev = T),
                            angle_col=90)
    
    return(plotHeatmap$gtable)
  })
  
  #-----------------------------------------------------------------
  # Show references table (metadata about datasets)
  output$references <- renderDT({
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
      paste0("MuscleModels_plot_", Sys.Date(), ".png")
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
      paste0("MuscleModels_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- selectedGeneData()
      if (!is.null(df)) {
        write.table(df, file, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)
        cat("\n# Data generated on MuscleOmics.org\n", file = file, append = TRUE)
      }
    }
  )
  
}
