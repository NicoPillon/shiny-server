#--------------------------------------------------------------------------------------------------------
#
# Overview
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
  # Load sections one after the next
  
  # Triggers
  start_aging <- reactiveVal(FALSE)
  start_obesity <- reactiveVal(FALSE)
  start_fibers <- reactiveVal(FALSE)
  
    # Once aging plot is shown, start obesity computation
  observeEvent(input$plot_aging_done, {
    start_obesity(TRUE)
  })
  
  # Once obesity plot is shown, start fibers
  observeEvent(input$plot_obesity_done, {
    start_fibers(TRUE)
  })
  
  #----------------------------------------------------------------------------------------
  # Gene description
  get_gene_description <- function(gene_symbol) {
    isolate({
      start_aging(FALSE)
      start_obesity(FALSE)
      start_fibers(FALSE)
    })
    
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
    isolate({
      start_obesity(FALSE)
      start_fibers(FALSE)
    })
    
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

    # plot
    ggplot(plotdata, aes(x=age_group, y=genedata)) +
      geom_boxplot(outlier.size = 0.1, fill = "gray80", alpha = 0.5)  + 
      geom_sina(size = 0.5, position = position_dodge(0), alpha = 0.25) +
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
    p <- dat_human_muscle_aging()
    htmlwidgets::onRender(p, "
    function(el, x) {
      Shiny.setInputValue('plot_aging_done', Math.random());
    }
  ")
  })

  #----------------------------------------------------------------------------------------
  # Human muscle obesity
  dat_human_muscle_obesity <- reactive({
    req(input$inputGeneSymbol)
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

    # plot
    ggplot(plotdata, aes(x=bmi_category, y=genedata)) +  
      geom_boxplot(outlier.size = 0.1, fill = "gray80", alpha = 0.5)  + 
      geom_sina(size = 0.5, position = position_dodge(0), alpha = 0.25) +
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
    p <- dat_human_muscle_obesity()
    htmlwidgets::onRender(p, "
    function(el, x) {
      Shiny.setInputValue('plot_obesity_done', Math.random());
    }
  ")
  })
  
  #----------------------------------------------------------------------------------------
  # Fiber types
  dat_fiber_types <- reactive({
    req(input$inputGeneSymbol)
    genename <- "NR4A3"
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
    
    # Find which file contains the gene
    gene_to_file <- readRDS("../fiber_types/data/gene_list_transcriptome.Rds")
    file_row <- gene_to_file[gene_to_file$SYMBOL == genename, ]
    req(nrow(file_row) > 0)  # stop if gene not found
    
    file_path <- file.path("../fiber_types/data", file_row$file)
    
    # Read the file and extract the gene row
    df <- arrow::read_feather(file_path, as_data_frame = FALSE)
    
    gene_index <- file_row$row
    req(length(gene_index) == 1)
    
    selected_row <- df[gene_index, ] %>% as.data.frame()
    rownames(selected_row) <- genename
    
    metadata <- readRDS("../fiber_types/data/metadata_transcriptome.Rds")
    plotdata <- data.frame(metadata, 
                           genedata = as.numeric(selected_row))
    
    # exclude mixed fibers
    plotdata <- plotdata[!grepl("Mixed", plotdata$FiberType),]
    
    ggplot(plotdata, aes(x = FiberType, y = genedata, fill = FiberType)) +
      geom_boxplot(position = position_dodge(0.8), outlier.size = 0) +
      geom_sina(size = 0.5, position = position_dodge(0.8), alpha = 0.5) +
      theme_bw(base_size = 16) + 
      theme(legend.position = "none") +
      labs(x = "Fiber Type", 
           y = "Relative expression, log2") +
      scale_y_continuous(expand = c(0, 4)) +
      scale_fill_manual(values = c("Type I" = "#8B0000", 
                                   "Type IIA" = "#F5DEB3", 
                                   "Type IIX"= "#D3D3D3"))
  })
  
  output$plot_fiber_types <- renderPlot({
    dat_fiber_types()
  })
  
}
