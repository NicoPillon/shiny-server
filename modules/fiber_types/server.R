#-----------------------------------------------------------------------
#
# Profile of skeletal muscle fibers
#
#----------------------------------------------------------------------

server <- function(input, output, session) {

  # Code to send the height of the app to adjust iframe  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                         choices=gene_list_all, 
                         server=TRUE, 
                         selected=c("LDHA", "LDHB", "MYH7", "MYH1", "MYH2"), 
                         options=NULL)
  
  #-----------------------------------------------------------------
  # REACTIVE: load only selected gene(s) from dataset
  selectedGeneData <- reactive({
    req(input$inputGeneSymbol)
    genename <- c("LDHA", "LDHB", "MYH7", "MYH1", "MYH2")
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
    
    # Match gene to file and row info
    file_row <- gene_to_file_transcriptome[gene_to_file_transcriptome$SYMBOL %in% genename, ]
    
    # Split by file
    split_rows <- split(file_row, file_row$file)
    
    # Load available data
    selected_list <- lapply(split_rows, function(rows) {
      path <- file.path("data", unique(rows$file))
      df <- arrow::read_feather(path)
      
      selected <- df[rows$row, , drop = FALSE] %>%
        data.frame()
      
      # Ensure rownames are matched properly
      rownames(selected) <- rows$SYMBOL[order(rows$row)]
      selected <- selected[match(rows$SYMBOL, rownames(selected)), , drop = FALSE]
      
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
  # REACTIVE: load only selected gene(s) from dataset  
  selectedProteinData <- reactive({
    req(input$inputGeneSymbol)
    genename <- c("LDHA", "LDHB", "MYH7", "MYH1", "MYH2")
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
    
    # Match gene to file and row info
    file_row <- gene_to_file_proteome[gene_to_file_proteome$SYMBOL %in% genename, ]
    
    # Split by file
    split_rows <- split(file_row, file_row$file)
    
    # Load available data
    selected_list <- lapply(split_rows, function(rows) {
      path <- file.path("data", unique(rows$file))
      df <- arrow::read_feather(path)
      
      selected <- df[rows$row, , drop = FALSE] %>%
        data.frame()
      
      # Ensure rownames are matched properly
      rownames(selected) <- rows$SYMBOL[order(rows$row)]
      selected <- selected[match(rows$SYMBOL, rownames(selected)), , drop = FALSE]
      
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
  # Boxplot - transcriptome
  plotDataGene <- eventReactive(input$updatePlot, {
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))

    gene_data <- t(selected_row)
    gene_data <- as.data.frame(t(selectedGeneData()))
    colnames(gene_data) <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
    
    plotdata <- cbind(metadata_transcriptome, gene_data)
    
    # exclude mixed fibers
    plotdata <- plotdata[!grepl("Mixed", plotdata$FiberType),]
    
    # Extract and reshape data for plotting
    plotdata <- pivot_longer(plotdata, 
                             cols = genename, 
                             names_to = "gene", 
                             values_to = "data")
    
    plotdata$FiberType <- factor(plotdata$FiberType,
                                 levels = c("Type I", "Mixed Type I/II", 
                                            "Type IIA", "Mixed Type IIA/IIX", "Type IIX"))
    
    # Plot
    ggplot(plotdata, aes(x = gene, y = data, fill = FiberType)) +
      geom_boxplot(position = position_dodge(0.8), outlier.size = 0) +
      #geom_sina(size = 0.5, position = position_dodge(0.8), alpha = 0.5) +
      theme_bw(base_size = 16) + 
      labs(x = NULL, 
           y = "Relative expression, log2", 
           subtitle = "Transcriptome") +
      scale_y_continuous(expand = c(0, 4)) +
      scale_fill_manual(values = c("Type I" = "#8B0000", 
                                   "Type IIA" = "#F5DEB3", 
                                   "Type IIX"= "#D3D3D3", 
                                   "Mixed" = "#A0522D")) +
      stat_compare_means(aes(label = after_stat(p_value_formatter(..p..))),
                         parse = TRUE,
                         size = 4, 
                         vjust = -1)
    
    
  })
  
  output$GenePlot <- renderPlot({
    plotDataGene()
  })
  
  #-----------------------------------------------------------------
  # Boxplot - proteome
  plotDataProtein <- eventReactive(input$updatePlot, {
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
 
    protein_data <- as.data.frame(t(selectedProteinData()))
    colnames(protein_data) <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
    
    plotdata <- cbind(metadata_proteome,
                      protein_data)
    
    # exclude mixed fibers
    plotdata <- plotdata[!grepl("Mixed", plotdata$FiberType),]
    
    # Extract and reshape data for plotting
    plotdata <- pivot_longer(plotdata, 
                             cols = genename, 
                             names_to = "gene", 
                             values_to = "data")
    
    plotdata$FiberType <- factor(plotdata$FiberType,
                                 levels = c("Type I", "Mixed Type I/II", 
                                            "Type IIA", "Mixed Type IIA/IIX", "Type IIX"))
    
    # Plot
    ggplot(plotdata, aes(x = gene, y = data, fill = FiberType)) +
      geom_boxplot(position = position_dodge(0.8), outlier.size = 0) +
      #geom_sina(size = 0.5, position = position_dodge(0.8), alpha = 0.2) +
      theme_bw(base_size = 16) + 
      labs(x = NULL, 
           y = "Relative expression, log2", 
           subtitle = "Proteome") +
      scale_y_continuous(expand = c(0, 4)) +
      scale_fill_manual(values = c("Type I" = "#8B0000", "Type IIA" = "#F5DEB3", 
                                   "Type IIX"= "#D3D3D3", "Mixed" = "#A0522D")) +
      stat_compare_means(aes(label = after_stat(p_value_formatter(..p..))),
                         parse = TRUE,
                         size = 4, 
                         vjust = -1)
  })
  
  output$ProteinPlot <- renderPlot({
    plotDataProtein()
  })
  
  
  ##################################################################################################################
  #Dataset tables
  output$references <- renderDataTable(options=list(signif = 3),{
    DT::datatable(
      references, 
      escape = FALSE, 
      rownames = FALSE,
      options = list(
        columnDefs = list(
          list(className = 'dt-center', targets = 2),  # Align column to the center
          list(className = 'dt-center', targets = 3),  # Align column to the center
          list(className = 'dt-center', targets = 4)  # Align column to the right
          # Add more lines for additional columns if needed
        )
      )
    )
  })
  
  ##################################################################################################################
  # Download button
  output$downloadGeneData <- downloadHandler(
    filename = function() {
      paste0("FiberTypes_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      gene_data <- as.data.frame(t(selectedGeneData()))
      colnames(gene_data) <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
      gene_data <- cbind(metadata_transcriptome, gene_data)
      
      prot_data <- as.data.frame(t(selectedProteinData()))
      colnames(prot_data) <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
      prot_data <- cbind(metadata_proteome, prot_data)
      
      # Create a workbook
      wb <- createWorkbook()
      
      # Add sheets
      addWorksheet(wb, "Transcriptome")
      addWorksheet(wb, "Proteome")
      
      # Write data to sheets
      writeData(wb, "Transcriptome", gene_data)
      writeData(wb, "Proteome", prot_data)
      
      # Save workbook
      if (!is.null(wb)) {
        saveWorkbook(wb, file, overwrite = TRUE)
      }
    }
  )
  
}
