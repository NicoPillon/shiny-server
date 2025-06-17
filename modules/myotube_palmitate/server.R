#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle cells
#
#--------------------------------------------------------------------------------------------------------

server <- function(input, output, session) {

  # Code to send the height of the app to adjust iframe
  session$onFlushed(function() {
    session$sendCustomMessage("resizeFrame", list())
  }, once = FALSE)
    
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=gene_list$TARGET, 
                       server=TRUE, 
                       selected=c("PDK4", "PXMP4", "LSM2", "ANGPTL4", "CPT1A", "ACAA2"), 
                       options=NULL)
  
  #-----------------------------------------------------------------
  # REACTIVE: load only selected gene(s)
  selectedGeneData <- reactive({
    req(input$inputGeneSymbol)
    genename <- c("PDK4", "PXMP4", "LSM2", "ANGPTL4", "CPT1A", "ACAA2")
    genename <- toupper(input$inputGeneSymbol)
    
    # Match gene to file and row info
    file_row <- gene_list[gene_list$TARGET %in% genename, ]
    
    # Split by file
    split_rows <- split(file_row, file_row$file)
    
    # Load available data
    selected_list <- lapply(split_rows, function(rows) {
      path <- file.path("data", unique(rows$file))
      df <- arrow::read_feather(path)
      
      selected <- df[rows$row, , drop = FALSE] %>%
        data.frame()
      
      # Ensure rownames are matched properly
      rownames(selected) <- rows$TARGET[order(rows$row)]
      selected <- selected[match(rows$TARGET, rownames(selected)), , drop = FALSE]
      
      return(selected)
    })
    
    # Combine all available rows
    selected_row <- do.call(rbind, selected_list)
    rownames(selected_row) <- gsub(".*.feather\\.", "", rownames(selected_row))
    
    # Add missing genes as NA
    missing_genes <- setdiff(genename, rownames(selected_row))
    if (length(missing_genes) > 0) {
      missing_df <- matrix(NA, nrow = length(missing_genes), ncol = ncol(selected_row))
      rownames(missing_df) <- missing_genes
      colnames(missing_df) <- colnames(selected_row)
      selected_row <- rbind(selected_row, missing_df)
    }
    
    selected_row <- selected_row[genename, , drop = FALSE]  # preserve original input order
    data.frame(selected_row)
  })
  
  #-----------------------------------------------------------------
  # Boxplots
  plotDataBox <- eventReactive(input$updatePlot, {
    df <- selectedGeneData()
    
    #plot
    dat <- data.frame(metadata,
                      t(df))
    dat <- pivot_longer(dat, cols = c(12:ncol(dat)),
                        values_to = "y",
                        names_to = "Gene")
    

    
    #filter according to selected categories
    colnames(dat)
    dat <- dplyr::filter(dat,
                         concentration.micromolar >= input$concentration[1] & concentration.micromolar <= input$concentration[2],
                         time.hours >= input$duration[1] & time.hours <= input$duration[2],
                         cell.type %in% input$cell_type,
                         species %in% input$species,
                         )
    
    # Check if filtered data is empty
    validate(
      need(nrow(dat) > 0, "No data available for the selected filters. Please adjust your selections.")
    )
    
    # plot
    ggplot(dat, aes(x=Gene, y=y, fill=treatment)) + 
      geom_boxplot(outlier.size = 0.1, alpha = 0.5, position = position_dodge(0.8))  + 
      geom_sina(size = 1.5, position = position_dodge(0.8), alpha = 0.1) +
      theme_bw(16) + 
      #theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      labs(x="",
           y="mRNA expression (log2)",
           title=element_blank()) +
      scale_shape_manual(values=rep(c(15,16,17), 20)) +
      scale_fill_manual(values = c("#5B768E", "#bd1a0e")) +
      scale_y_continuous(expand = expansion(mult = c(0.05, .15)))
    
  })
  
  output$geneBoxplot <- renderPlot({
    plotDataBox()
  })
  
  #-----------------------------------------------------------------
  # Statistics table
  
  # Reactive: compute Wilcoxon test
  statisticsData <- eventReactive(input$updatePlot, {
    df <- selectedGeneData()
    
    # Merge metadata with gene expression
    dat <- data.frame(metadata, t(df))
    
    dat <- pivot_longer(dat, cols = c(12:ncol(dat)),
                        values_to = "data",
                        names_to = "Gene")
    
    # Filter according to user selection
    dat <- dplyr::filter(dat,
                         concentration.micromolar >= input$concentration[1] & concentration.micromolar <= input$concentration[2],
                         time.hours >= input$duration[1] & time.hours <= input$duration[2],
                         cell.type %in% input$cell_type,
                         species %in% input$species)
    
    # Check if data is empty
    validate(
      need(nrow(dat) > 0, "No data available for the selected filters. Please adjust your selections.")
    )
    
    # Calculate statistics
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
          TRUE ~ ""
        ),
        p_value = format(p_value, scientific = TRUE, digits = 2),
        FDR = format(FDR, scientific = TRUE, digits = 2)
      ) 
    
    # Convert Gene column into rownames (safe even if one row)
    stats_result <- stats_result %>% 
      as.data.frame() %>%
      tibble::column_to_rownames("Gene")
    
    stats_result <- data.frame(Statistics = colnames(stats_result),
                                                     t(stats_result))
    stats_result$Statistics <- gsub("_", " ", stats_result$Statistics)
    stats_result$Statistics <- gsub("logFoldChange", "log2(fold-change)", stats_result$Statistics)
    stats_result$Statistics <- gsub("FoldChange", "Fold-change", stats_result$Statistics)
    stats_result
    
    return(stats_result)
  })
  
  # Render table
  output$statistics <- DT::renderDataTable({
    dat <- statisticsData()
    DT::datatable(
      dat,
      escape = FALSE, 
      rownames = FALSE,
      options = list(
        searching = FALSE,   # Disable search bar
        paging = FALSE,      # Disable pagination
        info = FALSE,        # Hide "Showing X of Y entries"
        dom = 't'            # Only display table body
      )
      )
  })
  
  #-----------------------------------------------------------------
  # Dataset tables
  output$references <- renderDataTable({
    DT::datatable(
      references, 
      escape = FALSE, 
      rownames = FALSE,
      options = list(
        columnDefs = list(
          list(className = 'dt-center', targets = 1),  # Center column 1
          list(className = 'dt-center', targets = 2)   # Center column 2
        ),
        searching = FALSE,   # Disable search bar
        paging = FALSE,      # Disable pagination
        info = FALSE,        # Hide "Showing X of Y entries"
        dom = 't'            # Only display table body
      )
    )
  })
  
  #-----------------------------------------------------------------
  # Download button
  output$downloadGeneData <- downloadHandler(
    filename = function() {
      paste0("gene_expression_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- selectedGeneData()
      df <- data.frame(metadata,
                       t(df))

      if (!is.null(df)) {
        write.csv(df, file, row.names = TRUE)
      }
    }
  )
  
}
