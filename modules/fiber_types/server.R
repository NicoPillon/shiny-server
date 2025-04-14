#-----------------------------------------------------------------------
#
# Profile of skeletal muscle fibers
#
#----------------------------------------------------------------------

server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=gene_list_all, 
                       server=TRUE, 
                       selected=c("LDHA", "LDHB", "MYH7", "MYH1", "MYH2"), 
                       options=NULL)
  
  #-----------------------------------------------------------------
  # REACTIVE: load only selected gene(s) from dataset
  selectedGeneData <- reactive({
    req(input$inputGeneSymbol)
    
    # Parse gene names
    genename <- c("LDHA", "LDHB", "MYH7", "MYH1", "MYH2")
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
    
    # Match gene to file and row info
    file_row <- gene_to_file_transcriptome[gene_to_file_transcriptome$SYMBOL %in% genename, ]
    req(nrow(file_row) > 0)
    
    # Split by file (each group is a set of rows from one file)
    split_rows <- split(file_row, file_row$file)
    
    # Read and collect rows from each file
    selected_list <- lapply(split_rows, function(rows) {
      path <- file.path("data", unique(rows$file))
      df <- arrow::read_feather(path)
      selected <- df[rows$row, ]
      rownames(selected) <- rows$SYMBOL
      selected
    })
    
    # Combine into single data frame
    selected_row <- do.call(rbind, selected_list) %>%
      data.frame()
    rownames(selected_row) <- genename
    
    return(selected_row)
  })
  
  selectedProteinData <- reactive({
    req(input$inputGeneSymbol)
    
    # Parse gene names
    genename <- c("LDHA", "LDHB", "MYH7", "MYH1", "MYH2")
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
    
    # Match gene to file and row info
    file_row <- gene_to_file_proteome[gene_to_file_proteome$SYMBOL %in% genename, ]
    req(nrow(file_row) > 0)
    
    # Split by file (each group is a set of rows from one file)
    split_rows <- split(file_row, file_row$file)
    
    # Read and collect rows from each file
    selected_list <- lapply(split_rows, function(rows) {
      path <- file.path("data", unique(rows$file))
      df <- arrow::read_feather(path)
      selected <- df[rows$row, ]
      rownames(selected) <- rows$SYMBOL
      selected
    })
    
    # Combine into single data frame
    selected_row <- do.call(rbind, selected_list) %>%
      data.frame()
    rownames(selected_row) <- genename
    
    return(selected_row)
  })
  
  
  #-----------------------------------------------------------------
  # Boxplot - transcriptome
  plotDataGene <- eventReactive(input$updatePlot, {
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))

    plotdata <- data.frame(metadata_transcriptome, 
                           t(selectedGeneData()))
    
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
      geom_sina(size = 0.5, position = position_dodge(0.8), alpha = 0.5) +
      theme_bw(base_size = 16) + 
      labs(x = NULL, 
           y = "Relative expression, log2", 
           subtitle = "Transcriptome") +
      scale_y_continuous(expand = c(0, 4)) +
      scale_fill_manual(values = c("Type I" = "#8B0000", "Type IIA" = "#F5DEB3", "Type IIX"= "#D3D3D3", "Mixed" = "#A0522D"))
  })
  
  output$GenePlot <- renderPlot({
    plotDataGene()
  })
  
  #-----------------------------------------------------------------
  # Boxplot - proteome
  plotDataProtein <- eventReactive(input$updatePlot, {
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
 
    plotdata <- data.frame(metadata_proteome, 
                           t(selectedProteinData()))
    
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
      geom_sina(size = 0.5, position = position_dodge(0.8), alpha = 0.5) +
      theme_bw(base_size = 16) + 
      labs(x = NULL, 
           y = "Relative expression, log2", 
           subtitle = "Proteome") +
      scale_y_continuous(expand = c(0, 4)) +
      scale_fill_manual(values = c("Type I" = "#8B0000", "Type IIA" = "#F5DEB3", "Type IIX"= "#D3D3D3", "Mixed" = "#A0522D"))
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
  
}
