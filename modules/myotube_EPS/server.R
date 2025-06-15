#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic myotube response to electrical pulse stimulation
#
#--------------------------------------------------------------------------------------------------------

server <- function(input, output, session) {

  # Code to send the height of the app to adjust iframe
  session$onFlushed(function() {
    session$sendCustomMessage("resizeFrame", list())
  }, once = FALSE)
    
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=gene_to_file$TARGET, 
                       server=TRUE, 
                       selected=c("HMOX1", "NR4A3", "IL6", "KLF15", "TXNIP"), 
                       options=NULL)
  
  #-----------------------------------------------------------------
  # REACTIVE: load only selected target from dataset
  selectedGeneData <- reactive({
    req(input$inputGeneSymbol)
    genename <- c("HMOX1", "NR4A3", "IL6", "KLF15", "TXNIP")
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
    
    # Find which file contains the gene
    file_row <- gene_to_file[gene_to_file$TARGET %in% genename, ]
    req(nrow(file_row) > 0)  # stop if gene not found
    
    # get data
    df <- datamatrix[genename, ]
    df <- data.frame(geo_accession = colnames(df),
                     t(df))
    return(df)
  })
  
  #-----------------------------------------------------------------
  # Boxplots
  plotDataBox <- eventReactive(input$updatePlot, {
    df <- selectedGeneData()
    
    #plot
    plotdata <- right_join(metadata, 
                           df)
    plotdata$Condition <- factor(plotdata$Condition, levels = c("Control", "EPS"))
    
    # Identify pairIDs where POST Time is within the selected time range
    valid_pairs <- plotdata %>%
      filter(Condition == "EPS",
             time_hours %in% input$duration) %>%
      pull(pairID)
    
    # Filter the full dataset by valid pairIDs and other filters
    plotdata <- plotdata %>%
      filter(cell.type %in% input$cell_type,
             Pulse_parameters %in% input$pulse_param,
             pairID %in% valid_pairs)
    
    # Check if filtered data is empty
    validate(
      need(nrow(plotdata) > 0, "No data available for the selected filters. Please adjust your selections.")
    )
    
    # pivot
    plotdata <- pivot_longer(plotdata, cols = c(15:ncol(plotdata)),
                             values_to = "data",
                             names_to = "Gene")
    
    # plot
    ggplot(plotdata, aes(x = Gene, y = data, fill=Condition)) +
      geom_boxplot(outlier.size = 0.1, alpha = 0.5, position = position_dodge(0.8))  + 
      geom_sina(size = 1.5, position = position_dodge(0.8), alpha = 0.1) +
      theme_bw(16) + 
      labs(x="",
           y="mRNA expression (log2)",
           title=element_blank()) +
      scale_shape_manual(values=rep(c(15,16,17), 20)) +
      scale_fill_manual(values = c("#5B768E", "#bd1a0e")) +
      scale_y_continuous(expand = expansion(mult = c(0.05, .15))) +
      stat_compare_means(aes(label = after_stat(p_value_formatter(..p..))),
                         parse = TRUE,
                         size = 4, 
                         vjust = -1)
  })
  
  output$geneBoxplot <- renderPlot({
    plotDataBox()
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
      df <- right_join(metadata, 
                       df)

      if (!is.null(df)) {
        write.csv(df, file, row.names = TRUE)
      }
    }
  )
  
}
