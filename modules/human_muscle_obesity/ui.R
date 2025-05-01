#-----------------------------------------------------------------------
#
# Human Obesity
#
#----------------------------------------------------------------------

ui <- fluidPage(title="Obesity",
                style="padding:0%",
                
                # CSS for style
                tags$head(includeCSS("../../www/style.css")),
                
                # main page
                navbarPage("Obesity",
                           
                           tabPanel(
                             title = "Single target",
                             sidebarLayout(
                               sidebarPanel(width = 3,
                                            selectizeInput("inputTarget", "Target Name:", choices=NULL, multiple=F, width=1000),
                                            tags$em("Type the name of a gene, protein or metabolite and select your population of interest."),
                                            tags$hr(),
                                            sliderInput("age", tags$b("Age (years)"),
                                                        min = 18, max = 90, value = c(18,90), step = 1, sep = ""),
                                            checkboxGroupInput("diagnosis_diabetes", 
                                                               label = "Diabetes diagnosis", 
                                                               selected = c("Healthy", "Prediabetes", "T2D"),
                                                               choices = c("Healthy", "Prediabetes", "T2D")),
                                            em(h5("Prediabetes is defined as either impaired glucose tolerance measured during an OGTT or increased insulin 
                                              resistance measured with euglycemic hyperinsulinemic clamps.")),
                                            tags$hr(),
                                            tags$b("Statistics"),
                                            em("Wilcoxon ranked signed test comparing overweight/obesity to lean and Spearman correlation. The p values are not adjusted for multiple testing comparisons."),
                                            tags$hr(),
                                            downloadButton("downloadGeneData", "Download Data")
                                            
                               ),
                               mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                                         plotOutput("boxplot", height="500px") %>% withSpinner(color="#5B768E", type = 8),
                                         tags$hr(),
                                         plotOutput("correlationPlot", height="800px") %>% withSpinner(color="#5B768E", type = 8)
                               )
                             )
                           ),
                           tabPanel(
                             title = "Pathway",
                             
                             fluidRow(style="color:black;background-color:white;padding:0% 2% 0% 2%;",
                                      selectizeInput("inputPathway", "Target Names:", choices=NULL, multiple=T, width=1000),
                                      tags$em("Type the name of multiple genes, proteins and metabolites..."),
                             ),
                             fluidRow(style="color:black;background-color:white;padding:0% 2% 0% 2%;",
                                      column(width = 4,
                                             h3("Transcriptomics"),
                                             plotOutput("transcriptomicPlot", height="800px") %>% withSpinner(color="#5B768E", type = 8),
                                      ),
                                      column(width = 4,
                                             h3("Metabolomics"),
                                             plotOutput("metabolomicPlot", height="800px") %>% withSpinner(color="#5B768E", type = 8)
                                      ),
                                      column(width = 4,
                                             h3("Proteomics"),
                                             tags$em(p("Proteomics data coming soon, stay tuned!")),
                                             #plotOutput("proteomicPlot", height="800px") %>% withSpinner(color="#5B768E", type = 8),
                                      )
                             )
                           ),
                           
                           # Tab with datasets
                           tabPanel(
                             title = "Datasets",
                             tags$p(
                               tags$b("Have we missed a relevant study? Please "),
                               a("let us know", href = "mailto:nicolas.pillon@ki.se", target = "_blank"), "!"
                             ),
                             dataTableOutput("references"),
                             tags$br()
                             
                           ),
                           
                           # Tab with datasets
                           tabPanel(
                             title = "About",
                             h3("Methods"),
                             
                             p("Publicly available gene expression datasets from human skeletal muscle biopsies were collected from the Gene Expression Omnibus (GEO), including both RNA-seq and microarray platforms. RNA-seq raw count data were processed using the edgeR package, with gene filtering, TMM normalization, and voom transformation. Microarray CEL files were normalized using the RMA method in the oligo package and annotated with platform-specific databases."),
                             
                             p("Sample metadata were manually curated. Sex was validated based on the expression of the XIST and RPS4Y1 genes. Samples with inconsistent metadata or low average Spearman correlation were excluded as outliers."),
                             
                             p("All datasets were merged by gene symbol, log2-transformed, and centered on the median expression per gene. Batch effects between studies were corrected using the removeBatchEffect function from the limma package. Genes with over 10% missing values were excluded. Data quality was assessed using principal component analysis (PCA), heatmaps, boxplots, and expression histograms."),
                             
                             tags$hr(),
                             h3("Citation"),
                             "Preliminary analyses for beta-testing. Unpublished data")
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
