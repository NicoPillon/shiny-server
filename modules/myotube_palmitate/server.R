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

  observeEvent(input$resetInputs, {
    updateSelectizeInput(session, "inputGeneSymbol", selected = character(0))
    updateSliderInput(session, "concentration", value = c(100, 500))
    updateSliderInput(session, "duration", value = c(12, 96))
    updateCheckboxGroupInput(session, "cell_type", selected = c("C2C12", "LHCN-M2", "primary"))
    updateCheckboxGroupInput(session, "species", selected = c("human", "mouse", "rat"))
  })
  
  #-----------------------------------------------------------------
  # REACTIVE: load only selected gene(s)
  selectedGeneData <- reactive({
    req(input$inputGeneSymbol)
    genename <- c("PDK4", "PXMP4", "ANGPTL4", "CPT1A", "ACAA2")
    genename <- toupper(input$inputGeneSymbol)
    
    # Match gene to file and row info
    file_row <- gene_list[gene_list$TARGET %in% genename, ]
    
    # Split by file
    split_rows <- split(file_row, file_row$file)
    
    # # Load available data
    selected_list <- lapply(split_rows, function(rows) {
      path <- file.path("data", unique(rows$file))
      
      ds <- arrow::open_dataset(path, format = "parquet")
      
      selected <- ds %>%
        filter(TARGET %in% rows$TARGET) %>%
        collect() %>%
        as.data.frame()
      
      rownames(selected) <- selected$TARGET
      selected <- selected[match(rows$TARGET, rownames(selected)), , drop = FALSE]
      selected <- selected[, !(colnames(selected) %in% "TARGET"), drop = FALSE]
      
      return(selected)
    })
    
    # Combine all available rows
    selected_row <- do.call(rbind, selected_list)
    rownames(selected_row) <- file_row$TARGET
    
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
  plotDataBox <- reactive({
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
      theme_minimal(17, base_family = "Arial") + 
      theme(axis.text.x = element_text(face = "bold", size = 14, color = "black"),
            axis.text.y = element_text("black"),
            axis.title.y = element_text(face = "bold", color = "black")) +
      labs(x="",
           y="mRNA expression (log2)",
           fill = "Treatment") +
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
  statisticsData <- reactive({
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
          TRUE ~ "ns"
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
    stats_result$Statistics <- gsub("FDR", "FDR (Bonferroni)", stats_result$Statistics)
    stats_result
    
    return(stats_result)
  })
  
  # Render table
  output$statistics1 <- DT::renderDataTable({
    dat <- statisticsData()
    dat <- dat[!dat$Statistics %in% c("mean control", "sd control", "n control", "mean palmitate", "sd palmitate", "n palmitate"),]
    colnames(dat)[1] <- "Differential Expression Analysis"
    DT::datatable(
      dat,
      escape = FALSE, 
      rownames = FALSE,
      options = list(
        searching = FALSE,   # Disable search bar
        paging = FALSE,      # Disable pagination
        info = FALSE,        # Hide "Showing X of Y entries"
        ordering = FALSE,    # Disable column sorting
        dom = 't',           # Only display table body
        columnDefs = list(
          list(
            targets = 1:(ncol(dat)-1),  # center columns 2 onward (indexing starts at 0)
            className = 'dt-center'
          )
        )
      )
    )
  })
  
  # Render table
  output$statistics2 <- DT::renderDataTable({
    dat <- statisticsData()
    dat <- dat[dat$Statistics %in% c("mean control", "sd control", "n control", "mean palmitate", "sd palmitate", "n palmitate"),]
    colnames(dat)[1] <- "Group Summary Statistics"
    DT::datatable(
      dat,
      escape = FALSE, 
      rownames = FALSE,
      options = list(
        searching = FALSE,   # Disable search bar
        paging = FALSE,      # Disable pagination
        info = FALSE,        # Hide "Showing X of Y entries"
        ordering = FALSE,    # Disable column sorting
        dom = 't',           # Only display table body
        columnDefs = list(
          list(
            targets = 1:(ncol(dat)-1),  # center columns 2 onward (indexing starts at 0)
            className = 'dt-center'
          )
        )
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
  # Download button - plot
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("MyotubePalmitate_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      plot_obj <- plotDataBox()  # This returns a ggplot object
      plot_obj <- plot_obj + labs(caption = "Plot generated on MuscleOmics.org")
      
      # Open PNG device, print plot, close device
      png(filename = file, width = 3200, height = 1800, res = 300)
      print(plot_obj)
      dev.off()
    }
  )
  
  #-----------------------------------------------------------------
  # Download button - data
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("MyotubePalmitate_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- selectedGeneData()
      df <- data.frame(metadata, t(df))
      
      if (!is.null(df)) {
        # Write the data without row names
        write.table(df, file, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)
        
        # Append the footer line
        cat("\n# Data generated on MuscleOmics.org\n", file = file, append = TRUE)
      }
    }
  )
  
  #-----------------------------------------------------------------
  # Download button - stats
  output$downloadStats <- downloadHandler(
    filename = function() {
      paste0("MyotubePalmitate_statistics_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- statisticsData()
      
      if (!is.null(df)) {
        # Write the data without row names
        write.table(df, file, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)
        
        # Append the footer line
        cat("\n# Data generated on MuscleOmics.org\n", file = file, append = TRUE)
      }
    }
  )
  
}
