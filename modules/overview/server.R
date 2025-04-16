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

  # Code to collect gene name from loading page    
  observe({
    query <- parseQueryString(session$clientData$url_search)
    
    if (!is.null(query$gene)) {
      selected_gene <- toupper(query$gene)  # normalize to uppercase
      updateSelectizeInput(session, "inputGeneSymbol",
                           choices = gene_list,
                           selected = selected_gene,
                           server = TRUE)
    } else {
      updateSelectizeInput(session, "inputGeneSymbol",
                           choices = gene_list,
                           selected = "NR4A3",  # default gene
                           server = TRUE)
    }
  })
  

  #----------------------------------------------------------------------------------------
  # Gene description
  get_gene_description <- function(gene_symbol) {
    req(input$inputGeneSymbol)  # Makes sure input is available
    symbol <- input$inputGeneSymbol

    # Step 1: Query gene symbol
    res <- GET(paste0("https://mygene.info/v3/query?q=", symbol, "&species=human"))
    data <- fromJSON(content(res, "text", encoding = "UTF-8"))
    
    # Step 2: Extract the gene ID (select the first EntrezID)
    gene_id <- data$hits[[1]][1]
    
    # Step 3: Use the ID to get full gene info
    res2 <- GET(paste0("https://mygene.info/v3/gene/", gene_id))
    info <- fromJSON(content(res2, "text", encoding = "UTF-8"))
    
    # Return summary, or a fallback message
    if (!is.null(info$summary)) {
      return(info$summary)
    } else {
      return("No summary available.")
    }
  }
  
  output$gene_description <- renderText({
    req(input$inputGeneSymbol)
    desc <- get_gene_description(input$inputGeneSymbol)
    paste("Description:", desc)
  })
  
  #----------------------------------------------------------------------------------------
  # Human muscle aging
  dat_human_muscle_aging <- reactive({
    req(input$inputGeneSymbol)
    genename <- "NR4A3"
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
    
    # Find which file contains the gene
    gene_to_file <- readRDS("../human_muscle_aging/data/gene_list.Rds")
    file_row <- gene_to_file[gene_to_file$SYMBOL == genename, ]
    req(nrow(file_row) > 0)  # stop if gene not found
    
    file_path <- file.path("../human_muscle_aging/data", file_row$file)
    
    # Read the file and extract the gene row
    df <- arrow::read_feather(file_path, as_data_frame = FALSE)
    
    gene_index <- file_row$row
    req(length(gene_index) == 1)
    
    selected_row <- df[gene_index, ] %>% as.data.frame()
    rownames(selected_row) <- genename
    
    metadata <- readRDS("../human_muscle_aging/data/metadata.Rds")
    plotdata <- data.frame(metadata, 
                           genedata = as.numeric(selected_row))

    ggplot(plotdata, aes(x=age_group, y=genedata)) +
      geom_boxplot(outlier.size = 0.1, fill = "gray80", alpha = 0.5)  + 
      geom_sina(size = 1.5, position = position_dodge(0), alpha = 0.25) +
      theme_bw(16) + 
      theme(legend.position = "right", 
            legend.title = element_blank()) +
      #facet_wrap(.~sex, ncol = 2) +
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

  output$plot_human_muscle_aging <- renderPlot({
    dat_human_muscle_aging()
  })

  #-----------------------------------------------------------------------------------------------------------------
  # Human muscle obesity
  dat_human_muscle_obesity <- reactive({
    req(input$inputGeneSymbol)
    genename <- "NR4A3"
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
    
    # Find which file contains the gene
    gene_to_file <- readRDS("../human_muscle_obesity/data/gene_list.Rds")
    file_row <- gene_to_file[gene_to_file$SYMBOL == genename, ]
    req(nrow(file_row) > 0)  # stop if gene not found
    
    file_path <- file.path("../human_muscle_obesity/data", file_row$file)
    
    # Read the file and extract the gene row
    df <- arrow::read_feather(file_path, as_data_frame = FALSE)
    
    gene_index <- file_row$row
    req(length(gene_index) == 1)
    
    selected_row <- df[gene_index, ] %>% as.data.frame()
    rownames(selected_row) <- genename
    
    metadata <- readRDS("../human_muscle_obesity/data/metadata.Rds")
    plotdata <- data.frame(metadata, 
                           genedata = as.numeric(selected_row))
    
    plotdata$bmi_category <- factor(plotdata$bmi_category, 
                                    levels=c("Lean", "Overweight", "Obesity"))
    
    ggplot(plotdata, aes(x=bmi_category, y=genedata)) +  
      geom_boxplot(outlier.size = 0.1, fill = "gray80", alpha = 0.5)  + 
      geom_sina(size = 1.5, position = position_dodge(0), alpha = 0.25) +
      theme_bw(16) + 
      theme(legend.position = "right",
            legend.title = element_blank()) +
      #facet_wrap(.~sex, ncol = 2) +
      labs(x="BMI category",
           y="mRNA expression, log2") +
      scale_shape_manual(values=rep(c(15,16,17), 20)) +
      scale_color_manual(values = c("#5B768E", "orange", "darkred")) +
      scale_y_continuous(expand = expansion(mult = c(0, .15))) +
      stat_compare_means(aes(label = after_stat(p_value_formatter(..p..))),
                         ref.group = "Lean",
                         parse = TRUE,
                         size = 4, 
                         vjust = -1)
  })
  
  output$plot_human_muscle_obesity <- renderPlot({
    dat_human_muscle_obesity()
  })
  
  
}
