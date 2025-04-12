#-----------------------------------------------------------------------
#
# Human Obesity
#
#----------------------------------------------------------------------

server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=gene_list, 
                       server=TRUE, 
                       selected=c("LEP"), 
                       options=NULL)
  
  #-----------------------------------------------------------------
  # REACTIVE: load only selected gene(s) from dataset
  selectedGeneData <- reactive({
    req(input$inputGeneSymbol)
    genename <- "PDK4"
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
  output$genePlotObesity <- renderPlot({
    
    plotdata <- data.frame(metadata, 
                           genedata = as.numeric(selectedGeneData()))
    
    plotdata$bmi_category <- factor(plotdata$bmi_category, 
                                    levels=c("Lean", "Overweight", "Obesity"))

    #filter according to selected categories
    plotdata <- dplyr::filter(plotdata,
                              diagnosis %in% input$diagnosis_diabetes ,
                              age >= input$age[1] & age <= input$age[2])
    
    # label with n size
    plotdata$sex <- gsub("^male", paste0("Male, n = ", nrow(plotdata[plotdata$sex == "male",])), plotdata$sex)
    plotdata$sex <- gsub("^female", paste0("Female, n = ", nrow(plotdata[plotdata$sex == "female",])), plotdata$sex)
    
    # plot
    cowplot::plot_grid(
      ggplot(plotdata, aes(x=bmi, y=genedata)) +  
        geom_point(aes(color = diagnosis, shape = diagnosis), size = 3, alpha = 0.25)  + 
        geom_smooth(method = "lm", color = "black", se = FALSE) +
        theme_bw(16) +  
        theme(legend.position = "none") +
        facet_wrap(.~sex, ncol = 1) +
        labs(x="BMI, kg/m2",
             y="mRNA expression, log2") +
        scale_shape_manual(values=rep(c(15,16,17), 20)) +
        scale_color_manual(values = c("darkgreen", "orange", "darkred")) +
        scale_y_continuous(expand = expansion(mult = c(0, .15))) +
        stat_cor(size = 4, 
                 vjust = -1, 
                 label.x = 20),
      
      ggplot(plotdata, aes(x=bmi_category, y=genedata)) +  
        geom_boxplot(outlier.size = 0.1, fill = "gray80", alpha = 0.5)  + 
        geom_sina(aes(color = diagnosis, shape = diagnosis), 
                  size = 1.5, position = position_dodge(0), alpha = 0.25) +
        theme_bw(16) + 
        theme(legend.position = "right",
              legend.title = element_blank()) +
        facet_wrap(.~sex, ncol = 1) +
        labs(x="BMI category",
             y="mRNA expression, log2") +
        scale_shape_manual(values=rep(c(15,16,17), 20)) +
        scale_color_manual(values = c("darkgreen", "orange", "darkred")) +
        scale_y_continuous(expand = expansion(mult = c(0, .15))) +
        stat_compare_means(aes(label = after_stat(p_value_formatter(..p..))),
                           ref.group = "Lean",
                           parse = TRUE,
                           size = 4, 
                           vjust = -1),
      ncol = 2
    )
    
  })
  
  #-----------------------------------------------------------------
  #Dataset tables
  output$datasets <- renderDataTable(options=list(signif = 3),{
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
