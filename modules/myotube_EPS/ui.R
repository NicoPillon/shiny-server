#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle cells - UI
#
#--------------------------------------------------------------------------------------------------------

ui <- fluidPage(
  title = "MyotubeEPS",         # Title of the browser tab
  style = "padding:0%",              # Remove default padding around the page
  
  #---------------------------------------------------------
  # Load custom CSS for consistent style across the app
  tags$head(includeCSS("../../www/style.css")),
  
  #---------------------------------------------------------
  # Custom JS for iframe resizing
  tags$head(
    tags$script(HTML("
        Shiny.addCustomMessageHandler('resizeFrame', function(message) {
          const height = document.documentElement.scrollHeight;
          parent.postMessage({ frameHeight: height }, '*');
        });

        window.addEventListener('resize', function() {
          const height = document.documentElement.scrollHeight;
          parent.postMessage({ frameHeight: height }, '*');
        });
      "))
  ),
  
  #---------------------------------------------------------
  # Navigation bar layout with multiple tabs
  navbarPage(
    # Custom title: logo image + text
    title = HTML('
      <div style="display: flex; align-items: center;margin: -10px;">
        <img src="../../www/img/snippet/myotube_EPS.png" style="height: 40px; margin-right: 10px;">
        <span style="font-size: 20px; font-weight: bold; margin-right: 10px;">Electrical Pulse Stimulation</span>
      </div>
    '),
    
    #=====================================================================
    #============================ TAB 1 =================================
    #=====================================================================
    tabPanel("Analysis",  # Main tab for exploring data and results
             
             # Introductory text with link to the methods section
             p(HTML('Skeletal muscle cells response to electrical pulse stimulation using integrated transcriptomic data. \n To understand how the data was generated, read the 
             <a href="#" onclick="$(\'.navbar-nav a:contains(\\\'Methods\\\')\').click()">methods</a>.')),
             
             tags$br(),
             
             # Layout: Sidebar on the left + Main content on the right
             sidebarLayout(
               
               #-----------------------------
               # Sidebar panel with input filters
               sidebarPanel(width = 3,
                            h4(tags$b("Select your criteria")),
                            
                            # Gene selector (multi-select)
                            selectizeInput("inputGeneSymbol", 
                                           "Genes of Interest", 
                                           choices = NULL, 
                                           multiple = TRUE, 
                                           width = 600),
                            
                            # Sliders for EPS parameters
                            checkboxGroupInput("duration", 
                                               label = "Duration (hours)", 
                                               selected = c("3", "3+Rest3h", "6", "24", "48", "72", "72(1h+3rest)"),
                                               choices = c("3", "3+Rest3h", "6", "24", "48", "72", "72(1h+3rest)")),
                            checkboxGroupInput("pulse_param", 
                                               label = "Pulse Parameters", 
                                               selected = c("14V, 5Hz, 2ms", "3V, 1Hz, 6ms", "13V, 66Hz, 2ms",
                                                            "10V, 1Hz, 2ms", "30V, 1Hz, 2ms", "10V, 1Hz, 10ms", 
                                                            "11.5V, 1Hz, 2ms"),
                                               choices = c("14V, 5Hz, 2ms", "3V, 1Hz, 6ms", "13V, 66Hz, 2ms",
                                                           "10V, 1Hz, 2ms", "30V, 1Hz, 2ms", "10V, 1Hz, 10ms", 
                                                           "11.5V, 1Hz, 2ms")),
                            
                            # Checkbox filters for cell type and species
                            checkboxGroupInput("cell_type", 
                                               label = "Cell Type", 
                                               selected = c("mouse C2C12", "human primary"),
                                               choices = c("mouse C2C12", "human primary")),
                            
                            # Button to reset all filters
                            actionButton("resetInputs", "Reset Filters", icon = icon("undo")),
                            
                            tags$hr(),
                            
                            # Download section
                            h4(tags$b("Download")),
                            downloadButton("downloadPlot", "Plot (.png)"), tags$br(), 
                            downloadButton("downloadData", "Data (.csv)"), tags$br(), 
                            downloadButton("downloadStats", "Statistics (.csv)")
               ),
               
               #-----------------------------
               # Main panel with output: plot + tables
               mainPanel(width = 9, style = "padding:0% 4% 0% 1%;",
                         
                         # Plot of gene expression (boxplot)
                         plotOutput("geneBoxplot", height = "550px") %>% withSpinner(color = "#5B768E", type = 8),
                         
                         #tags$hr(),
                         tags$em(textOutput("filterSummaryText")),
                         
                         tags$hr(),
                         
                         # Table 1: Differential expression results
                         DT::DTOutput("statistics1"),
                         tags$br(),
                         
                         # Table 2: Summary statistics (mean, SD, n)
                         DT::DTOutput("statistics2"),
                         tags$br()
               )
             )
    ),
    
    #=====================================================================
    #============================ TAB 2 =================================
    #=====================================================================
    tabPanel("Methods",  # Tab describing methodology
             
             style = "padding:0% 5% 1% 5%;",
             
             # Introduction
             h3("Why this app?"),
             p("This app provides a centralized and user-friendly platform to explore how electrical pulse stimulation (EPS) influences the transcriptome of skeletal muscle cells. 
EPS mimics the physiological trigger of muscle contraction — a key event that reshapes gene expression and cellular function. By compiling and harmonizing data from multiple independent studies, 
this tool facilitates hypothesis generation, cross-species comparisons, and mechanistic insights relevant to muscle physiology, exercise biology, and metabolic health."),
             
             # Data processing and integration methodology
             h3("Methods"),
             p("Transcriptomic datasets from human, mouse, and rat skeletal myotubes exposed to electrical pulse stimulation (EPS) were collected from public repositories, including both microarray (Affymetrix) and RNA-sequencing platforms. 
Orthologous gene annotation was performed using the biomaRt package to harmonize gene identifiers across species. Raw expression values were median-centered within each dataset to mitigate platform-specific effects. 
All datasets were merged by gene symbol using ortholog mappings. Genes with low expression or missing values in most samples were filtered out, and outlier samples were removed based on average pairwise Pearson correlation across genes. 
Batch effects due to different experimental sources (GEO accession) were corrected using the `removeBatchEffect` function from the limma package. Differential expression between EPS and control conditions was evaluated using a paired linear model followed by moderated t-tests and Bonferroni correction. 
Study metadata, such as species, EPS parameters, and exposure duration, were extracted from GEO records and incorporated into the analysis."),
             
             # Statistical tests
             h3("Statistics"),
             p("Statistical analysis is based on the Wilcoxon signed-rank test comparing palmitate-treated samples to controls. 
        The results table presents both unadjusted and Bonferroni-adjusted p-values to correct for multiple testing across 
        all transcripts in the database, ensuring a highly conservative approach. Significance is indicated as follows: 
        * for FDR < 0.05, ** for FDR < 0.01, and *** for FDR < 0.001. 'ns' denotes non-significant results."),
             
             # Citation
             h3("Citation"),
             p("Pillon NJ et al. Unpublished."),
             
             # References table for included datasets
             h3("Datasets"),
             tags$p(
               tags$em("Have we missed a relevant study? Please ",
                       a("let us know!", href = "mailto:nicolas.pillon@ki.se", target = "_blank")
               )
             ),
             dataTableOutput("references")
    )
  )
)

