#-----------------------------------------------------------------------
#
# Profile of skeletal muscle fibers
#
#----------------------------------------------------------------------

ui <- fluidPage(
  title="FiberTypes",    # Title of the browser tab
  style="padding:0%",    # Remove default padding around the page
  
  #---------------------------------------------------------
  # Load custom CSS for consistent style across the app
  tags$head(includeCSS("../../www/style.css")),
  
  #---------------------------------------------------------
  # Navigation bar layout with multiple tabs
  navbarPage(
    # Custom title: logo image + text
    title = HTML('
      <div style="display: flex; align-items: center;margin: -10px;">
        <img src="../../www/img/snippet/fiber_types.png" style="height: 40px; margin-right: 10px;">
        <span style="font-size: 20px; font-weight: bold; margin-right: 10px;">Fiber Types</span>
      </div>
    '),
    
    #=====================================================================
    #============================ TAB 1 =================================
    #=====================================================================
    tabPanel("Analysis",
             
             # Introductory text with link to the methods section
             p(HTML('Proteome and transcriptome of human skeletal muscle fibers. \n To understand how the data was generated, read the 
             <a href="#" onclick="$(\'.navbar-nav a:contains(\\\'Methods\\\')\').click()">methods</a>.')),
             
             tags$br(),
             
             # Layout: Sidebar on the left + Main content on the right
             sidebarLayout(
               
               #-----------------------------
               # Sidebar panel with input filters
               sidebarPanel(width = 3,
                            h4(tags$b("Select your criteria")),
                            
                            # Gene selector (multi-select)
                            selectizeInput("inputTarget", "Targets of Interest", choices=NULL, multiple=T, width=600),

                            # select dataset
                            selectInput(
                              inputId = "inputOmics",
                              label = "OMICS",
                              choices = c("Proteome", "Transcriptome"),
                              selected = c("Transcriptome"),
                              multiple = FALSE,
                              selectize = TRUE
                            ),
                            
                            # Checkbox filters for fibers of interest
                            checkboxGroupInput("fibers", 
                                               label = "Fiber Types", 
                                               selected = c("Type I", "Type IIA", "Type IIX"),
                                               choices = c("Type I", "Type IIA", "Type IIX")),
                            
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
                         
                         plotOutput("TargetPlot", height="400px") %>% withSpinner(color="#5b768e", type = 8),

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
             p("This analysis integrates proteomic and transcriptomic datasets of human skeletal muscle fibers to identify and compare fiber-type specific signatures."),
             
             # Data processing and integration methodology
             h3("Methods"),
             p("For proteomics, raw intensity values were retrieved from supplementary Excel files, and metadata were cleaned and harmonized. Protein expression matrices were constructed by aligning sample identifiers and names. Rows with multiple name mappings were split and only those with fewer missing values were retained. Fiber types were annotated based on the relative abundance of MYH isoforms, using published criteria: Type I (MYH7 ≥ 80%), Type IIA (MYH2 ≥ 80%), and Type IIX (MYH1 ≥ 20%)."),
             
             p("For transcriptomics, raw UMI counts or CEL files were preprocessed with standard Bioconductor pipelines. Lowly expressed genes were filtered, normalized using TMM or RMA, and batch effects were removed with `limma::removeBatchEffect`. Gene identifiers were mapped to gene symbols. Metadata were extracted and harmonized to enable sample-wise comparison. Single muscle fibers were annotated based on the expression levels of MYH7, MYH1, and MYH2. Normalized expression values were back-transformed to linear scale and used to compute isoform contributions. Fibers were classified into 'Type I', 'Type IIA', or 'Type IIX' based on tertiles of isoform expression, with intermediate profiles labeled as 'Mixed' and excluded."),
             
             # Statistical tests
             h3("Statistics"),
             p("When all fiber types are selected, Type IIA and Type IIX fibers are grouped together as 'Type II' to enable a comparison with Type I fibers. If only Type IIA or only Type IIX is selected along with Type I, the comparison is performed between those specific fiber types. When both Type IIA and Type IIX are selected, a dedicated comparison is performed between these two fast-twitch fiber types. For each target, summary statistics (mean, standard deviation, and sample size) are computed per fiber type. Differential expression is assessed using the Wilcoxon rank-sum test, and p-values are adjusted using the Bonferroni method to correct for multiple testing across all genes in the dataset. Adjusted p-values are shown as FDR (Bonferroni). Statistical significance is indicated as follows: * for FDR < 0.05, ** for FDR < 0.01, and *** for FDR < 0.001. 'ns' indicates a non-significant result."),
             
             # Citation
             h3("Citation"),
             "Pillon NJ. Unpublished analysis",
             
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