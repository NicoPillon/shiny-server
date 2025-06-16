#-----------------------------------------------------------------------
#
# Profile of skeletal muscle fibers
#
#----------------------------------------------------------------------

ui <- fluidPage(title="FiberTypes",
                style="padding:0%",
                
                # CSS for style
                tags$head(includeCSS("../../www/style.css")),
                
                # main page
                navbarPage("FiberTypes",
                           
                           tabPanel(
                             title = "Plots",
                             
                             # main page
                             sidebarLayout(
                               sidebarPanel(width = 3,
                                            selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                                            actionButton("updatePlot", "Refresh plot", icon("refresh")),
                                            tags$hr(),
                                            tags$b("Statistics:"),
                                            em("Kruskal–Wallis tests to assess whether at least one fiber type differs from the others. Reported p-values are not adjusted for multiple comparisons."),
                                            tags$hr(),
                                            downloadButton("downloadGeneData", "Download Data")
                               ),
                               mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                                         plotOutput("GenePlot", height="400px") %>% withSpinner(color="#5b768e", type = 8),
                                         plotOutput("ProteinPlot", height="400px") %>% withSpinner(color="#5b768e", type = 8)
                               )
                             )
                           ),
                           
                           # description of methods
                           tabPanel("Description",
                                    
                                    h3("Citation"),
                                    "Unpublished analysis",
                                    
                                    tags$hr(),
                                    
                                    h3("Methods"),
                                    
                                    p("This analysis integrates proteomic and transcriptomic datasets of human skeletal muscle fibers to identify and compare fiber-type specific signatures."),
                                    
                                    p("For proteomics, raw intensity values were retrieved from supplementary Excel files, and metadata were cleaned and harmonized. Protein expression matrices were constructed by aligning sample identifiers and gene symbols. Rows with multiple gene mappings were split and only those with fewer missing values were retained. Fiber types were annotated based on the relative abundance of MYH isoforms, using published criteria: Type I (MYH7 ≥ 80%), Type IIA (MYH2 ≥ 80%), and Type IIX (MYH1 ≥ 20%)."),
                                    
                                    p("For transcriptomics, raw UMI counts or CEL files were preprocessed with standard Bioconductor pipelines. Lowly expressed genes were filtered, normalized using TMM or RMA, and batch effects were removed with `limma::removeBatchEffect`. Gene identifiers were mapped to gene symbols. Metadata were extracted and harmonized to enable sample-wise comparison. Single muscle fibers were annotated based on the expression levels of MYH7, MYH1, and MYH2. Normalized expression values were back-transformed to linear scale and used to compute isoform contributions. Fibers were classified into 'Type I', 'Type IIA', or 'Type IIX' based on tertiles of isoform expression, with intermediate profiles labeled as 'Mixed' and excluded."),
                                    
                                    tags$hr(),
                                    
                                    h3("Datasets Included in the Analysis"),
                                    tags$p(
                                      tags$b("Are we missing a relevant study? Please "),
                                      a("let us know!", href = "mailto:nicolas.pillon@ki.se", target = "_blank")
                                    ),
                                    dataTableOutput("references")
                           ),
                           
                           # Code to send height to resizing iframe
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