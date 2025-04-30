#-----------------------------------------------------------------------
#
# Human Aging in skeletal muscle
#
#----------------------------------------------------------------------

server <- function(input, output, session) {

  # Code to send the height of the app to adjust iframe
  session$onFlushed(function() {
    session$sendCustomMessage("resizeFrame", list())
  }, once = FALSE)

  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=gene_list, 
                       server=TRUE, 
                       selected=c("MYH8"), 
                       options=NULL)
  
  #-----------------------------------------------------------------
  # REACTIVE: load only selected gene(s) from dataset
  selectedGeneData <- reactive({
    req(input$inputGeneSymbol)
    genename <- "MYH8"
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
    
    # Find which file contains the gene
    file_row <- gene_to_file[gene_to_file$SYMBOL == genename, ]
    req(nrow(file_row) > 0)  # stop if gene not found
    
    file_path <- file.path("data", file_row$file)
    
    # Read the file and extract the gene row
    df <- arrow::read_feather(file_path, as_data_frame = FALSE)
    
    gene_index <- file_row$row
    req(length(gene_index) == 1)
    
    selected_row <- df[gene_index, ] %>% as.data.frame()
    rownames(selected_row) <- genename
    
    return(selected_row)
  })
  
    #-----------------------------------------------------------------
    # Boxplot
    output$geneBoxPlotAging <- renderPlot({
      
      plotdata <- data.frame(metadata, 
                             genedata = as.numeric(selectedGeneData()))
      
      #filter according to selected categories
      plotdata <- dplyr::filter(plotdata,
                                diagnosis %in% input$diagnosis_sarcopenia)
      
      
      # label with n size
      plotdata$sex <- gsub("^Male", paste0("Male, n = ", nrow(plotdata[plotdata$sex == "Male",])), plotdata$sex)
      plotdata$sex <- gsub("^Female", paste0("Female, n = ", nrow(plotdata[plotdata$sex == "Female",])), plotdata$sex)
      
      ggplot(plotdata, aes(x=age_group, y=genedata)) +
          geom_boxplot(outlier.size = 0.1, fill = "gray80", alpha = 0.5)  + 
          geom_sina(aes(color = diagnosis, shape = diagnosis), 
                    size = 1.5, position = position_dodge(0), alpha = 0.25) +
          theme_bw(16) + 
          theme(legend.position = "right", 
                legend.title = element_blank()) +
          facet_wrap(.~sex, ncol = 2) +
          labs(x="Age group",
               y="mRNA expression, log2") +
          scale_shape_manual(values=rep(c(15,16,17), 20)) +
          scale_color_manual(values = c("#5B768E", "#bd1a0e")) +
          scale_y_continuous(expand = expansion(mult = c(0, .15))) +
          stat_compare_means(aes(label = after_stat(p_value_formatter(..p..))),
                             ref.group = "40-60",
                             parse = TRUE,
                             size = 4, 
                             vjust = -1)
    })
    
    #-----------------------------------------------------------------
    # Correlation
    output$geneCorrelationPlotAging <- renderPlot({
      
      plotdata <- data.frame(metadata, 
                             genedata = as.numeric(selectedGeneData()))
      
      #filter according to selected categories
      plotdata <- dplyr::filter(plotdata,
                                diagnosis %in% input$diagnosis_sarcopenia)
      
      
      # label with n size
      plotdata$sex <- gsub("^Male", paste0("Male, n = ", nrow(plotdata[plotdata$sex == "Male",])), plotdata$sex)
      plotdata$sex <- gsub("^Female", paste0("Female, n = ", nrow(plotdata[plotdata$sex == "Female",])), plotdata$sex)
      
      ggplot(plotdata, aes(x=age, y=genedata)) +
          geom_jitter(aes(color = diagnosis, shape = diagnosis), 
                      width = 5, size = 3, alpha = 0.5)  + 
          geom_smooth(method = "lm", color = "black", se = FALSE) +
          theme_bw(16) +  
        theme(legend.position = "right", 
              legend.title = element_blank()) +
          facet_wrap(.~sex, ncol = 1) +
          labs(x="Age, years",
               y="mRNA expression, log2") +
          scale_shape_manual(values=rep(c(15,16,17), 20)) +
          scale_color_manual(values = c("#5B768E", "#bd1a0e")) +
          scale_y_continuous(expand = expansion(mult = c(0, .15))) +
          stat_cor(size = 4, 
                   vjust = -1, 
                   label.x = 20)
      
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
      df <- data.frame(metadata,
                       t(selectedGeneData()))
      if (!is.null(df)) {
        write.csv(df, file, row.names = TRUE)
      }
    }
  )
  
}