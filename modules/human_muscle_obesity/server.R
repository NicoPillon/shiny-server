#-----------------------------------------------------------------------
#
# Human Obesity
#
#----------------------------------------------------------------------

server <- function(input, output, session) {

  # Code to send the height of the app to adjust iframe
  session$onFlushed(function() {
    session$sendCustomMessage("resizeFrame", list())
  }, once = FALSE)
    
  updateSelectizeInput(session, 'inputTarget', 
                       choices=gene_list, 
                       server=TRUE, 
                       selected=c("LEP"), 
                       options=NULL)
  
  updateSelectizeInput(session, 'inputPathway', 
                       choices=gene_list, 
                       server=TRUE, 
                       selected=c("GRB14", "LPL", "MRC2", "FRZB", "NMNAT1", "ACAP1", "TIMM50"),
                       options=NULL)
  
  #-----------------------------------------------------------------
  # REACTIVE: load only selected target from dataset
  selectedTargetData <- reactive({
    req(input$inputTarget)
    genename <- "PDK4"
    genename <- toupper(unlist(strsplit(input$inputTarget, "[,;\\s]+")))
    
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
  # REACTIVE: load only selected target from dataset
  selectedPathwayData <- reactive({
    req(input$inputPathway)
    genename <- c("LEP", "ADIPOQ", "MRC1", "ITGAM")
    genename <- toupper(unlist(strsplit(input$inputPathway, "[,;\\s]+")))
    
    # Match gene to file and row info
    file_row <- gene_to_file[gene_to_file$SYMBOL %in% genename, ]
    
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
  # Boxplot
  output$boxplot <- renderPlot({

    plotdata <- data.frame(metadata, 
                           genedata = as.numeric(selectedTargetData()))
    
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
    ggplot(plotdata, aes(x=bmi_category, y=genedata)) +  
        geom_boxplot(outlier.size = 0.1, fill = "gray80", alpha = 0.5)  + 
        geom_sina(aes(color = diagnosis, shape = diagnosis), 
                  size = 1.5, position = position_dodge(0), alpha = 0.25) +
        theme_bw(16) + 
        theme(legend.position = "right",
              legend.title = element_blank()) +
        facet_wrap(.~sex, ncol = 2) +
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
 
  #----------------------------------------------------------------- 
  # Correlation
  output$correlationPlot <- renderPlot({
    
    plotdata <- data.frame(metadata, 
                           genedata = as.numeric(selectedTargetData()))
    
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
    ggplot(plotdata, aes(x=bmi, y=genedata)) +  
        geom_point(aes(color = diagnosis, shape = diagnosis), size = 3, alpha = 0.5)  + 
        geom_smooth(method = "lm", color = "black", se = FALSE) +
        theme_bw(16) +  
      theme(legend.position = "right", 
            legend.title = element_blank()) +
        facet_wrap(.~sex, ncol = 1) +
        labs(x="BMI, kg/m2",
             y="mRNA expression, log2") +
        scale_shape_manual(values=rep(c(15,16,17), 20)) +
        scale_color_manual(values = c("#5B768E", "orange", "darkred")) +
        scale_y_continuous(expand = expansion(mult = c(0, .15))) +
        stat_cor(size = 4, 
                 vjust = -1, 
                 label.x = 20)
    
  })
  
  #----------------------------------------------------------------- 
  # Pathway
  output$transcriptomicPlot <- renderPlot({
    
    target_names <- c("GRB14", "LPL", "MRC2", "FRZB", "NMNAT1", "ACAP1", "TIMM50")
    target_names <- input$inputPathway
    transcriptomics <- readRDS("data/stats.transcriptomics.Rds")
    transcriptomics <- transcriptomics %>%
      subset(SYMBOL %in% target_names)
    
    # Create the data frame
    transcriptomics <- data.frame(
      SYMBOL = transcriptomics$SYMBOL,
      logFC = transcriptomics$logFC,
      adj.P.Val = transcriptomics$adj.P.Val
    )
    
    # Create a factor variable for coloring
    transcriptomics <- transcriptomics %>%
      mutate(
        pval_category = case_when(
          adj.P.Val < 0.01 ~ "FDR < 0.01",
          adj.P.Val < 0.05 ~ "FDR < 0.05",
          adj.P.Val < 0.1  ~ "FDR < 0.1",
          TRUE             ~ "ns"
        ),
        pval_category = factor(pval_category, levels = c("ns", "FDR < 0.1", "FDR < 0.05", "FDR < 0.01"))
      )
    
    # Plot
    ggplot(transcriptomics, aes(x = logFC, y = reorder(SYMBOL, logFC), fill = pval_category)) +
      geom_col(color = "black") +
      geom_vline(xintercept = 0, color = "black") +
      scale_fill_manual(
        values = c(
          "FDR < 0.01" = "#7f0000",  # deep dark red
          "FDR < 0.05" = "#d7301f",  # strong red
          "FDR < 0.1"  = "#fc9272",  # light red/pink
          "ns"         = "lightgray"
        ),
        name = element_blank()
      ) +
      labs(x = "log2(Fold-Change)",
           y = NULL,
           title = "mRNA") +
      theme_minimal(base_size = 16) +
      theme(legend.position = "bottom")
  })
  
  
  #-----------------------------------------------------------------
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
                       t(selectedTargetData()))
      if (!is.null(df)) {
        write.csv(df, file, row.names = TRUE)
      }
    }
  )
  
}
