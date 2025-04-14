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
                       choices=gene_list, 
                       server=TRUE, 
                       selected=c("GENE1", "GENE10", "GENE22"), 
                       options=NULL)
  
  #-----------------------------------------------------------------
  # REACTIVE: load only selected gene(s) from dataset
  selectedGeneData <- reactive({
    req(input$inputGeneSymbol)
    genename <- toupper(input$inputGeneSymbol)
    
    # Match gene names
    matched_ids <- which(gene_list %in% genename)
    
    # Load only the required genes (rows)
    df <- arrow::read_feather("data/datamatrix.feather", as_data_frame = FALSE)[matched_ids, ] %>%
      as.data.frame()
    rownames(df) <- gene_list[matched_ids]
    
    return(df)
  })
  
  #-----------------------------------------------------------------
  # Boxplot
  plotDataBox <- eventReactive(input$updatePlot, {
    df <- selectedGeneData()
    
    dat <- data.frame(sample_list,
                      t(df))
    dat <- pivot_longer(dat, cols = rownames(df),
                        values_to = "data",
                        names_to = "gene")
    
    ggplot(dat, aes(x=condition, y=data, fill=species)) + 
      geom_boxplot(outlier.size = 0.1, fill = "gray80", alpha = 0.5)  + 
      geom_sina(aes(color = species), size = 1.5, position = position_dodge(0), alpha = 0.1) +
      facet_wrap(~gene) +
      theme_bw(16) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(x="",
           y="mRNA expression (log2)",
           title=element_blank())
    
  })
  
  output$geneBoxplot <- renderPlot({
    plotDataBox()
  })
  
  #-----------------------------------------------------------------
  # Dataset tables
  output$references <- renderDataTable(options=list(signif = 3),{
    DT::datatable(
      references, 
      escape = FALSE, 
      rownames = FALSE,
      options = list(
        columnDefs = list(
          list(className = 'dt-center', targets = 1),  # Align column to the center
          list(className = 'dt-center', targets = 2)  # Align column to the center
          # Add more lines for additional columns if needed
        )
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
      df <- data.frame(sample_list,
                       t(df))
      if (!is.null(df)) {
        write.csv(df, file, row.names = TRUE)
      }
    }
  )
  
}
