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
  # Detect unavailable genes
  observe({
    req(input$inputGeneSymbol)
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))

    # Check across all datasets
    found_in_any <- any(genename %in% gene_list)
    
    if (!found_in_any) {
      output$gene_warning <- renderUI({
        tags$div(style = "color: red; font-weight: bold; padding-top: 10px;",
                 paste("⚠️ Sorry, the gene", genename, "was not found in any dataset. Please check your spelling."))
      })
      
      # optionally: prevent downstream analysis
      isolate({
        start_aging(FALSE)
        start_obesity(FALSE)
        start_fibers(FALSE)
      })
      
    } else {
      output$gene_warning <- renderUI({ NULL })  # Clear the warning
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
  # Gene Structure from alphafold
  # output$structureViewer <- renderR3dmol({
  #   # Step 1: Get UniProt ID for input gene
  #   symbol <- input$inputGeneSymbol
  #   ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #   uniprot_id <- gene_map$uniprot_gn_id[gene_map$hgnc_symbol == input$inputGeneSymbol][1]
  # 
  #   # Step 2: Download AlphaFold PDB
  #   pdb_path <- paste0(tempdir(), "/", uniprot_id, ".pdb")
  #   url <- paste0("https://alphafold.ebi.ac.uk/files/AF-", uniprot_id, "-F1-model_v4.pdb")
  #   tryCatch({
  #     download.file(url, destfile = pdb_path, quiet = TRUE)
  #     pdb_text <- paste(readLines(pdb_path), collapse = "\n")
  #     
  #     # Step 3: Display it with r3dmol
  #     r3dmol() %>%
  #       m_add_model(data = pdb_text, format = "pdb") %>%
  #       m_set_style(style = m_style_cartoon(color = "spectrum")) %>%
  #       m_zoom_to()
  #   }, error = function(e) {
  #     showNotification("Could not find structure.", 
  #                      type = "error")
  #     return(NULL)
  #   })
  # })
  
  #----------------------------------------------------------------------------------------
  # Module fiber_types
  dat_fiber_types_RNA <- reactive({
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
    
    plot_RNA <- ggplot(plotdata, aes(x = FiberType, y = genedata, fill = FiberType)) +
      geom_boxplot(position = position_dodge(0.8), 
                   outlier.size = 0) +
      geom_sina(position = position_dodge(0.8),
                size = 0.5,
                alpha = 0.25) +
      theme_bw(base_size = 16) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(x = element_blank(),
           y = "mRNA expression, log2") +
      scale_y_continuous(expand = expansion(mult = c(0, .15))) +
      scale_fill_manual(values = c("Type I" = "#8B0000", 
                                   "Type IIA" = "#F5DEB3", 
                                   "Type IIX"= "#D3D3D3")) +
      stat_compare_means(aes(label = after_stat(p_value_formatter(..p..))),
                         parse = TRUE,
                         label.x.npc = "center",
                         size = 4, 
                         vjust = -1,
                         hjust = 0.5)
  })
  
  dat_fiber_types_Prot <- reactive({
    req(input$inputGeneSymbol)
    genename <- "CKM"
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
    
    # Find which file contains the gene
    gene_to_file <- readRDS("../fiber_types/data/gene_list_proteome.Rds")
    file_row <- gene_to_file[gene_to_file$SYMBOL == genename, ]
    req(nrow(file_row) > 0)  # stop if gene not found
    
    file_path <- file.path("../fiber_types/data", file_row$file)
    
    # Read the file and extract the gene row
    df <- arrow::read_feather(file_path, as_data_frame = FALSE)
    
    gene_index <- file_row$row
    req(length(gene_index) == 1)
    
    selected_row <- df[gene_index, ] %>% as.data.frame()
    rownames(selected_row) <- genename
    
    metadata <- readRDS("../fiber_types/data/metadata_proteome.Rds")
    plotdata <- data.frame(metadata, 
                           genedata = as.numeric(selected_row))
    
    # exclude mixed fibers
    plotdata <- plotdata[!grepl("Mixed", plotdata$FiberType),]
    
    plot_Prot <- ggplot(plotdata, aes(x = FiberType, y = genedata, fill = FiberType)) +
      geom_boxplot(position = position_dodge(0.8), 
                   outlier.size = 0) +
      geom_sina(position = position_dodge(0.8),
                size = 0.5,
                alpha = 0.25) +
      theme_bw(base_size = 16) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(x = element_blank(),
           y = "Protein expression, log2") +
      scale_y_continuous(expand = expansion(mult = c(0, .15))) +
      scale_fill_manual(values = c("Type I" = "#8B0000", 
                                   "Type IIA" = "#F5DEB3", 
                                   "Type IIX"= "#D3D3D3")) +
      stat_compare_means(aes(label = after_stat(p_value_formatter(..p..))),
                         parse = TRUE,
                         label.x.npc = "center",
                         size = 4, 
                         vjust = -1,
                         hjust = 0.5)
  })
  
  output$plot_fiber_types <- renderPlot({
    req(input$inputGeneSymbol)
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))[1]
    
    # Safe function to wrap a plot call and fallback to message if error
    safe_plot <- function(expr, label) {
      tryCatch(
        expr,
        error = function(e) {
          ggplot() + 
            annotate("text", x = 0.5, y = 0.5, label = paste(label, "\nnot available"), size = 4, hjust = 0.5) +
            theme_void()
        }
      )
    }
    
    # Run both plot functions safely
    plot_rna <- safe_plot(dat_fiber_types_RNA(), "RNA data")
    plot_prot <- safe_plot(dat_fiber_types_Prot(), "Protein data")
    
    # Combine plots
    cowplot::plot_grid(plot_rna, plot_prot, ncol = 2)
  })
  
  
  
  #----------------------------------------------------------------------------------------
  # Module human_muscle_aging
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
    ggplot(plotdata, aes(x=sex, y=genedata, fill = age_group)) +
      geom_boxplot(position = position_dodge(0.8), 
                   outlier.size = 0) +
      geom_sina(position = position_dodge(0.8),
                size = 0.5,
                alpha = 0.25) +
      theme_bw(16) + 
      theme(legend.position = "right") +
      #facet_wrap(.~sex, ncol = 2) +
      labs(x=NULL,
           y="mRNA expression, log2",
           fill = "Age") +
      scale_shape_manual(values=rep(c(15,16,17), 20)) +
      scale_fill_manual(values = c("orangered1", "orangered3", "orangered4")) +
      scale_y_continuous(expand = expansion(mult = c(0, .15))) +
      stat_compare_means(aes(label = after_stat(p_value_formatter(..p..))),
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
  # Module human_muscle_obesity
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
    ggplot(plotdata, aes(x=sex, y=genedata, fill=bmi_category)) +  
      geom_boxplot(position = position_dodge(0.8), 
                   outlier.size = 0) +
      geom_sina(position = position_dodge(0.8),
                size = 0.5,
                alpha = 0.25) +
      theme_bw(16) + 
      theme(legend.position = "right",
            legend.title = element_blank()) +
      #facet_wrap(.~sex, ncol = 2) +
      labs(x=NULL,
           y="mRNA expression, log2") +
      scale_shape_manual(values=rep(c(15,16,17), 20)) +
      scale_fill_manual(values = c("gold1", "gold3", "gold4")) +
      scale_y_continuous(expand = expansion(mult = c(0, .15))) +
      stat_compare_means(aes(label = after_stat(p_value_formatter(..p..))),
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
  # Module MetaMEx
  dat_metamex_human <- reactive({
    req(input$inputGeneSymbol)
    gene_symbol <- "NR4A3"
    gene_symbol <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))
    
    # Collect data
    limma_list <- MetaMEx_human
    plot_data <- lapply(names(limma_list), function(condition) {
      res <- limma_list[[condition]]
      if (gene_symbol %in% rownames(res)) {
        gene_row <- res[gene_symbol, ]
        data.frame(
          Condition = condition,
          logFC = gene_row$logFC,
          CI.Lower = gene_row$CI.L,
          CI.Upper = gene_row$CI.R,
          P.Value = gene_row$P.Value
        )
      } else {
        data.frame(
          Condition = condition,
          logFC = NA,
          CI.Lower = NA,
          CI.Upper = NA,
          P.Value = NA
        )
      }
    }) %>% bind_rows()
    
    # Plot
    ggplot(plot_data, aes(x = Condition, y = logFC, fill = Condition)) +
      geom_bar(stat = "identity", na.rm = TRUE, linewidth = 0.2, color = "black") +
      geom_errorbar(aes(ymin = CI.Lower, ymax = CI.Upper), width = 0.2, na.rm = TRUE) +
      geom_hline(yintercept = 0) +
      theme_bw(16) + theme(legend.position = "none") +
      labs(
        y = "log2 Fold Change (±95% CI)",
        x = element_blank()
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_viridis_d()
  })
  
  output$plot_metamex <- renderPlot({
    p <- dat_metamex_human()
    print(p)
  })
  
  
  #----------------------------------------------------------------------------------------
  # Module muscle_models
  dat_muscle_models <- reactive({
    req(input$inputGeneSymbol)
    genename <- c("NR4A3")
    genename <- toupper(input$inputGeneSymbol)
    
    # Match gene names
    genelist <- readRDS("../muscle_models/data/gene_names.Rds")
    matched_ids <- which(genelist %in% genename)
    
    if (length(matched_ids) == 0) {
      showNotification("No genes matched the database.", type = "error")
      return(NULL)
    }
    
    # Load only the required genes (rows)
    df <- arrow::read_feather("../muscle_models/data/datamatrix.feather", as_data_frame = FALSE)[matched_ids, ] %>%
      as.data.frame()
    rownames(df) <- genelist[matched_ids]
    
    #plot
    samples_list <- readRDS("../muscle_models/data/samples_list.Rds")
    dat <- data.frame(samples_list,
                      t(df))
    dat <- pivot_longer(dat, cols = c(7:ncol(dat)),
                        values_to = "y",
                        names_to = "Gene")
    
    
    ggplot(dat, aes(x=cell_tissue, y=y, fill=species)) + 
      geom_boxplot(position = position_dodge(0.8), 
                   outlier.size = 0) +
      geom_sina(position = position_dodge(0.8),
                size = 0.5,
                alpha = 0.25) +
      theme_bw(16) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(x="",
           y="mRNA expression, log2") +
      scale_y_continuous(breaks = round(seq(-4, 8, by=2),1)) +
      geom_hline(aes(yintercept=0), linetype="dashed", show.legend=F, color="gray60") +
      scale_color_manual(values=samples_list$species_colors)
  })
  
  output$plot_muscle_models <- renderPlot({
    p <- dat_muscle_models()
    htmlwidgets::onRender(p, "
    function(el, x) {
      Shiny.setInputValue('plot_obesity_done', Math.random());
    }
  ")
  })
  
}
