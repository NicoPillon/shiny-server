#--------------------------------------------------------------------------------------------------------
#
# Overview
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
  # Code to collect gene name from loading page    
  observe({
    query <- parseQueryString(session$clientData$url_search)
    
    if (!is.null(query$gene)) {
      selected_gene <- toupper(query$gene)  # normalize to uppercase
      updateSelectizeInput(session, "inputGeneSymbol",
                           choices = gene_list$SYMBOL,
                           selected = selected_gene,
                           server = TRUE)
    } else {
      updateSelectizeInput(session, "inputGeneSymbol",
                           choices = gene_list$SYMBOL,
                           selected = "PDK4",  # default gene
                           server = TRUE)
    }
  })
  
  
  #----------------------------------------------------------------------------------------
  # Detect unavailable genes
  observe({
    req(input$inputGeneSymbol)
    genename <- toupper(unlist(strsplit(input$inputGeneSymbol, "[,;\\s]+")))

    # Check across all datasets
    found_in_any <- any(genename %in% gene_list$SYMBOL)
    
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
  
  #-----------------------------------------------------------------
  # REACTIVE: Load only selected gene(s) from on-disk parquet files
  selectedGeneData <- reactive({
    req(input$inputGeneSymbol)  # Ensure input is not empty
    genename <- "PDK4"
    genename <- toupper(input$inputGeneSymbol)  # Convert input to uppercase for consistency
    
    # Find matching file(s) and row(s) for each gene
    file_row <- gene_list[gene_list$SYMBOL %in% genename, ]
    split_rows <- split(file_row, file_row$file)  # Group by file for efficient loading
    split_rows
    
    # Load selected genes from each file using lazy loading via Arrow
    selected_list <- lapply(split_rows, function(rows) {
      path <- file.path("data", unique(rows$file))  # Path to parquet dataset
      ds <- arrow::open_dataset(path, format = "parquet")  # Lazy load dataset
      selected <- ds %>%
        filter(TARGET %in% rows$TARGET) %>%
        collect() %>%
        as.data.frame()
      
      # Reorder rows to match input gene order
      rownames(selected) <- selected$TARGET
      selected <- selected[match(rows$TARGET, rownames(selected)), , drop = FALSE]
      selected <- selected[, !(colnames(selected) %in% "TARGET"), drop = FALSE]
      return(selected)
    })
    
    # Combine rows across files
    selected_row <- do.call(rbind, selected_list)
    rownames(selected_row) <- file_row$TARGET
    return(selected_row)
  })
  
  
  #----------------------------------------------------------------------------------------
  # Overview bar plot
  output$overviewBarHighcharter <- renderHighchart({
    req(selectedGeneData())
    
    #dat <- selected_row
    dat <- selectedGeneData()
    
    # Add Significance category
    dat <- dat %>%
      mutate(Significance = case_when(
        adj.P.Val < 0.001 ~ "FDR < 0.001",
        adj.P.Val < 0.01  ~ "FDR < 0.01",
        adj.P.Val < 0.05  ~ "FDR < 0.05",
        TRUE              ~ "ns"
      )) %>%
      mutate(Significance = factor(Significance, levels = c("FDR < 0.001", "FDR < 0.01", "FDR < 0.05", "ns")))
    
    # Sort and preserve order
    dat <- dat %>%
      arrange(desc(logFC)) %>%
      mutate(experiment = factor(experiment, levels = experiment))
    
    # Define color mapping for significance
    signif_colors <- c(
      "FDR < 0.001" = "darkgreen",
      "FDR < 0.01"  = "orange",
      "FDR < 0.05"  = "yellow",
      "ns"          = "gray"
    )
    
    # Convert to list of points for highcharter
    data_points <- purrr::pmap(dat, function(experiment, logFC, adj.P.Val, url, Significance, ...) {
      list(
        y = round(logFC,2),
        FDR = signif(adj.P.Val, 2),
        name = experiment,
        color = signif_colors[[Significance]],
        #url = url,
        url = paste0(url, "?gene=", URLencode(input$inputGeneSymbol))
      )
    })
    
    hc_theme_custom <- hc_theme(
      chart = list(
        style = list(fontFamily = "Arial")
      ),
      title = list(
        style = list(fontFamily = "Arial")
      )
    )
    
    highchart() %>%
      hc_add_theme(hc_theme_custom) %>%
      hc_chart(
        type = "bar",
        plotBorderWidth = 1,
        plotBorderColor = NULL
      ) %>%
      hc_xAxis(
        type = "category",
        gridLineWidth = 1,
        gridLineColor = "#e0e0e0",
        labels = list(style = list(color = "black", 
                                   fontSize = "14px",
                                   whiteSpace = "nowrap",
                                   textOverflow = "none",
                                   overflow = "allow"))
      ) %>%
      hc_yAxis(
        title = list(text = "Fold-change (log2)", style = list(color = "black", fontSize = "14px")),
        labels = list(style = list(color = "gray", fontSize = "12px")),
        gridLineWidth = 1,
        gridLineColor = "#e0e0e0",
        plotLines = list(
          list(value = 0, width = 2, color = "black", zIndex = 5)
        )
      ) %>%
      hc_add_series(
        name = "Experiments",
        data = data_points,
        showInLegend = FALSE
      ) %>%
      hc_plotOptions(series = list(
        cursor = "pointer",
        point = list(events = list(
          click = JS("function() { window.open(this.url, '_blank'); }")
        ))
      )) %>%
      hc_tooltip(
        useHTML = TRUE,
        pointFormat = "Click to explore this dataset!</b>"
      )
    
  })
  
  #----------------------------------------------------------------------------------------
  # Gene description
  output$gene_description <- renderUI({
    req(input$inputGeneSymbol)
    symbol <- toupper(input$inputGeneSymbol)
    
    # Initialize defaults
    ncbi_text <- "🔬 <b>NCBI</b>: No match found."
    uniprot_text <- "🔬 <b>UniProt</b>: No match found."
    genecards_link <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", symbol)
    genecards_text <- paste0("🔬 <b>GeneCards</b>: <a href='", genecards_link, "' target='_blank'>View entry for ", symbol, "</a>")
    
    # --- NCBI Summary ---
    try({
      search <- entrez_search(db = "gene", term = paste0(symbol, "[Gene Name] AND Homo sapiens[Organism]"))
      if (length(search$ids) > 0) {
        summary <- entrez_summary(db = "gene", id = search$ids[1])
        ncbi_text <- paste0("🔬 <b>NCBI</b>:<ul><li>", summary$summary, "</li></ul>")
      }
    }, silent = TRUE)
    
    # --- UniProt Comments (all values) ---
    try({
      search_url <- paste0(
        "https://rest.uniprot.org/uniprotkb/search?",
        "query=(", symbol, ")+AND+organism_id:9606+AND+reviewed:true",
        "&fields=accession&format=json"
      )
      res <- httr::GET(search_url)
      if (res$status_code == 200) {
        hits <- content(res, as = "parsed", simplifyVector = FALSE)$results
        if (length(hits) > 0) {
          acc <- hits[[1]]$primaryAccession
          entry_url <- paste0("https://rest.uniprot.org/uniprotkb/", acc, ".json")
          entry_res <- httr::GET(entry_url)
          if (entry_res$status_code == 200) {
            comments <- content(entry_res, as = "parsed", simplifyVector = FALSE)$comments
            if (!is.null(comments)) {
              values <- unlist(lapply(comments, function(cmt) {
                if (!is.null(cmt$texts)) cmt$texts[[1]]$value else NULL
              }))
              if (length(values) > 0) {
                uniprot_text <- paste0(
                  "🔬 <b>UniProt</b>:<ul><li>",
                  paste(values, collapse = "</li><li>"),
                  "</li></ul>"
                )
              }
            }
          }
        }
      }
    }, silent = TRUE)
    
    # --- Render all HTML ---
    HTML(paste(
      "<div style='margin-top:15px;'>",
      ncbi_text, "<br>",
      uniprot_text, "<br>",
      genecards_text,
      "</div>"
    ))
  })
  
  
  
}
