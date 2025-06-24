#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle cells - UI
#
#--------------------------------------------------------------------------------------------------------

ui <- fluidPage(
  title = "MyotubeT2D",         # Title of the browser tab
  style = "padding:0%",              # Remove default padding around the page
  
  #---------------------------------------------------------
  # Load custom CSS for consistent style across the app
  #tags$head(includeCSS("../../www/style.css")),
  
  #---------------------------------------------------------
  # Navigation bar layout with multiple tabs
  navbarPage(
    # Custom title: logo image + text
    title = HTML('
      <div style="display: flex; align-items: center;margin: -10px;">
        <img src="../../www/img/snippet/myotube_T2D.png" style="height: 40px; margin-right: 10px;">
        <span style="font-size: 20px; font-weight: bold; margin-right: 10px;">Type 2 Diabetes</span>
      </div>
    '),
    
    #=====================================================================
    #============================ TAB 1 =================================
    #=====================================================================
    tabPanel("Analysis",  # Main tab for exploring data and results
             
             # Introductory text with link to the methods section
             p(HTML('Transcriptomic profile of skeletal muscle cells grown from individuals with type 2 diabetes or healthy controls. \n To understand how the data was generated, read the 
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
                            
                            # Sliders for palmitate concentration and exposure duration
                            checkboxGroupInput("sex", 
                                               label = "Sex",
                                               choices = c("Male", "Female"),
                                               selected = c("Male", "Female")),

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
                         DT::dataTableOutput("statistics1"),
                         tags$br(),
                         
                         # Table 2: Summary statistics (mean, SD, n)
                         DT::dataTableOutput("statistics2"),
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
             p("This app provides an interactive resource to investigate how gene expression differs in primary human skeletal muscle cells from individuals with and without type 2 diabetes (T2D). 
         By pooling and harmonizing multiple publicly available transcriptomic studies—including both myoblasts and myotubes—this tool enables exploration of diabetes-related molecular signatures 
         in skeletal muscle cells, with applications for research, diagnosis, and therapeutic development."),
             
             # Data processing and integration methodology
             h3("Methods"),
             p("Publicly available transcriptomic datasets were downloaded from the Gene Expression Omnibus (GEO), including raw data from microarray and RNA-seq platforms. All datasets 
        represent primary human skeletal muscle cells (myoblasts and/or myotubes) derived from individuals with type 2 diabetes and healthy controls."),
             p("Microarray data were normalized using Robust Multi-array Average (RMA) and filtered based on expression intensity. RNA-seq datasets were processed using variance-stabilizing 
        transformation (VST) following filtering of low-count genes (≥10 counts in ≥50% of samples). For Illumina bead arrays, detection p-values were used to retain reliably expressed genes."),
             p("Gene annotation was standardized to gene symbols, and data were aggregated at the gene level. Each dataset was median-centered and quantile-normalized. 
        Only genes detected in at least 50% of the studies were retained. Batch correction was performed using study ID (GEO accession) as a covariate, and outlier removal was based 
        on low inter-sample correlation (<1st percentile). Quality control included boxplots, histograms, correlation matrices, and PCA."),
             
             
             # Statistical tests
             h3("Statistics"),
             p("Statistical analysis is based on the Wilcoxon signed-rank test comparing palmitate-treated samples to controls. 
        The results table presents both unadjusted and Bonferroni-adjusted p-values to correct for multiple testing across 
        all transcripts in the database, ensuring a highly conservative approach. Significance is indicated as follows: 
        * for FDR < 0.05, ** for FDR < 0.01, and *** for FDR < 0.001. 'ns' denotes non-significant results."),
             
             # Citation
             h3("Citation"),
             p("Pillon NJ. Unpublished."),
             
             # References table for included datasets
             h3("Datasets"),
             tags$p(
               tags$em("Have we missed a relevant study? Please ",
                       a("let us know!", href = "mailto:nicolas.pillon@ki.se", target = "_blank")
               )
             ),
             dataTableOutput("references")
    ),
    
    #=====================================================================
    #======================= IFRAME HEIGHT SYNC ==========================
    #=====================================================================
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
    )
  )
)
