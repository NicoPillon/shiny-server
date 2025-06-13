#-----------------------------------------------------------------------
#
# Human Aging in skeletal muscle
#
#----------------------------------------------------------------------

ui <- fluidPage(title="muscleAging",
                style="padding:0%",
                
                # CSS for style
                tags$head(includeCSS("../../www/style.css")),
                
                # main page
                navbarPage("muscleAging",
                           
                           tabPanel(
                             title = "Plots",
                             
                             # main page
                             sidebarLayout(
                               sidebarPanel(width = 3,
                                            selectizeInput("inputGeneSymbol", "Gene Symbol:", choices=NULL, multiple=F, width=1000),
                                            checkboxGroupInput("diagnosis_sarcopenia", 
                                                               label = "Sarcopenia diagnosis", 
                                                               selected = c("Healthy", "Sarcopenia"),
                                                               choices = c("Healthy", "Sarcopenia")),
                                            tags$hr(),
                                            tags$b("Statistics:"),
                                            em(("Spearman correlation and Wilcoxon ranked signed test comparing ages to the middle group. The p values are not adjusted for multiple testing comparisons.")),
                                            tags$hr(),
                                            downloadButton("downloadGeneData", "Download Data")
                                            
                               ),
                               mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                                         plotOutput("geneBoxPlotAging", height="500px") %>% withSpinner(color="#5B768E", type = 8),
                                         tags$hr(),
                                         plotOutput("geneCorrelationPlotAging", height="700px") %>% withSpinner(color="#5B768E", type = 8)
                               )
                             ),
                             
                           ),
                           
                           
                           # description of methods
                           tabPanel("Description",
                                    
                                    h3("Citation"),
                                    "Unpublished data.",
                                    tags$hr(),
                                    
                                    h3("Methods"),
                                    
                                    p("Publicly available transcriptomic datasets profiling human skeletal muscle were retrieved from the Gene Expression Omnibus (GEO), including both RNA sequencing (RNA-seq) and Affymetrix microarray platforms (Human Exon 1.0 ST and HGU133Plus2). Metadata for each dataset, including age, sex, diagnosis, and platform, were manually curated. Only datasets that reported individual age data were retained."),
                                    
                                    p("Raw count data from GTEx and multiple GEO datasets were aggregated by gene and annotated using the Homo.sapiens Bioconductor package. Gene expression filtering and normalization were performed using edgeR, followed by log2 transformation and median centering. For microarrays, raw CEL files were processed with the oligo package, normalized using RMA, filtered for low expression, annotated, and collapsed by gene symbol."),
                                    
                                    p("All expression matrices were merged by gene symbol and batch-corrected using limma. Principal component analysis was used to assess integration and sample clustering. Sex assignments were verified using expression of sex-specific markers (RPS4Y1 and XIST) and adjusted where discrepancies were observed. Outliers were detected based on pairwise sample correlations, and samples below the 1st percentile of average correlation were excluded. Final metadata and normalized expression matrices were saved for visualization and analysis."),
                                    tags$hr(),
                                    
                                    h3("Datasets Included in the Analysis"),
                                    tags$p(
                                      tags$b("Are we missing a relevant study? Please "),
                                      a("let us know!", href = "mailto:nicolas.pillon@ki.se", target = "_blank")
                                    ),
                                    dataTableOutput("references")
                           )
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