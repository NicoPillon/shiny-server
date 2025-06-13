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
                       choices=genelist$human_SYMBOL, 
                       server=TRUE, 
                       selected=c("PDK4", "PXMP4", "LSM2", "ANGPTL4", "CPT1A", "ACAA2"), 
                       options=NULL)
  
  # REACTIVE: load only selected gene(s)
  selectedGeneData <- reactive({
    req(input$inputGeneSymbol)
    genename <- c("PDK4", "PXMP4", "LSM2", "ANGPTL4", "CPT1A", "ACAA2")
    genename <- toupper(input$inputGeneSymbol)
    
    # Match gene names
    matched_ids <- which(genelist$human_SYMBOL %in% genename)
    
    if (length(matched_ids) == 0) {
      showNotification("No genes matched the database.", type = "error")
      return(NULL)
    }
    
    # Load only the required genes (rows)
    df <- arrow::read_feather("data/datamatrix.feather", as_data_frame = FALSE)[matched_ids, ] %>%
      as.data.frame()
    rownames(df) <- genelist$human_SYMBOL[matched_ids]
    
    return(df)
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
      df <- data.frame(metadata,
                       t(df))

      if (!is.null(df)) {
        write.csv(df, file, row.names = TRUE)
      }
    }
  )
  
}
