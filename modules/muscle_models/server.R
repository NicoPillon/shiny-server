#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle cells
#
#--------------------------------------------------------------------------------------------------------

server <- function(input, output, session) {

  #-----------------------------------------------
  # Resize the iframe containing the app
  # Sends a custom message to the parent HTML to adjust height dynamically
  session$onFlushed(function() {
    session$sendCustomMessage("resizeFrame", list())
  }, once = FALSE)
    
  #-----------------------------------------------
  # Code to collect target name from loading page    
  observe({
    query <- parseQueryString(session$clientData$url_search)
    
    if (!is.null(query$target)) {
      selected_target <- toupper(query$target)  # normalize to uppercase
      updateSelectizeInput(session, "inputGeneSymbol",
                           choices = gene_list$TARGET,
                           selected = selected_target,
                           server = TRUE)
    } else {
      updateSelectizeInput(session, "inputGeneSymbol",
                           choices = gene_list$TARGET,
                           selected = c("FN1", "SPP1", "MYL4", "MYH7", "ATP2A1", "MYL3"),  # default genes
                           server = TRUE)
    }
  })
  
  # REACTIVE: load only selected gene(s)
  selectedGeneData <- reactive({
    req(input$inputGeneSymbol)
    genename <- c("FN1", "SPP1", "MYL4", "MYH7", "ATP2A1", "MYL3")
    genename <- toupper(input$inputGeneSymbol)
    
    # Match gene names
    matched_ids <- which(genelist %in% genename)
    
    if (length(matched_ids) == 0) {
      showNotification("No genes matched the database.", type = "error")
      return(NULL)
    }
    
    # Load only the required genes (rows)
    df <- arrow::read_feather("data/datamatrix.feather", as_data_frame = FALSE)[matched_ids, ] %>%
      as.data.frame()
    rownames(df) <- genelist[matched_ids]
    
    return(df)
  })
  
  #-----------------------------------------------------------------
  # Boxplots
  plotDataBox <- eventReactive(input$updatePlot, {
    df <- selectedGeneData()
    
    #plot
    dat <- data.frame(samples_list,
                      t(df))
    dat <- pivot_longer(dat, cols = c(7:ncol(dat)),
                        values_to = "y",
                        names_to = "Gene")
    
    
    ggplot(dat, aes(x=cell_tissue, y=y, fill=species)) + 
      geom_boxplot(outlier.size = 0.1, fill = "gray80", alpha = 0.5)  + 
      geom_sina(aes(color = species), size = 1.5, position = position_dodge(0), alpha = 0.1) +
      facet_wrap(~Gene) +
      theme_bw(16) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(x="",
           y="mRNA expression (log2)",
           title=element_blank()) +
      scale_y_continuous(breaks = round(seq(-4, 8, by=2),1)) +
      geom_hline(aes(yintercept=0), linetype="dashed", show.legend=F, color="gray60") +
      scale_color_manual(values=samples_list$species_colors)
    
  })
  
  output$geneBoxplot <- renderPlot({
    plotDataBox()
  })
  
  
  #-----------------------------------------------------------------
  # Heatmap
  plotDataHeatmap <- eventReactive(input$updatePlot, {

    #make matrix with found gene names
    matrixData <- selectedGeneData()
    
    # Add gene names as a column for reshaping
    matrix_long <- matrixData %>%
      rownames_to_column("gene") %>%
      pivot_longer(
        cols = -gene,
        names_to = "sample",
        values_to = "expression"
      )
    matrix_long <- matrix_long %>%
      left_join(samples_list, by = c("sample" = "sample"))
    gene_medians <- matrix_long %>%
      group_by(gene, cell_tissue) %>%
      summarise(median_expression = median(expression, na.rm = TRUE), .groups = "drop")
    gene_matrix <- gene_medians %>%
      pivot_wider(
        names_from = cell_tissue,
        values_from = median_expression
      ) %>%
      column_to_rownames("gene")  # Optional: make gene names the rownames
    return(gene_matrix)
  })
  
  
  output$geneHeatmap <- renderPlot({
    
    #df_mean <- gene_matrix
    df_mean <- plotDataHeatmap()
    
    my_sample_col <- samples_list[,c(2,5,4)]
    my_sample_col <- my_sample_col[!duplicated(my_sample_col$cell_tissue),]
    my_sample_col <- data.frame(my_sample_col, row.names = 1)
    
    my_colour = list(
      #model = setNames(unique(samples_names$model.color), unique(samples_names$model)),
      species = setNames(rev(unique(samples_list$species_color)), unique(samples_list$species))
    )
    
    plotHeatmap <- pheatmap(df_mean, scale = "row",
                            annotation_col = my_sample_col,
                            annotation_colors = my_colour,
                            display_numbers = F,
                            cluster_rows = T,
                            cluster_cols = T,
                            fontsize=11,
                            color = hcl.colors(50, "BluYl", rev = T),
                            angle_col=90)
    
    return(plotHeatmap$gtable)
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
      df <- data.frame(samples_list,
                       t(df))
      df$sample <- NULL
      df$species_colors <- NULL
      if (!is.null(df)) {
        write.csv(df, file, row.names = TRUE)
      }
    }
  )
  
}
