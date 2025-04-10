#-----------------------------------------------------------------------
#
# Human Aging in skeletal muscle
#
#----------------------------------------------------------------------

ui <- fluidPage(# Google analytics
  tags$head(includeScript("../../google-analytics.html")),
  
  # CSS for style
  tags$head(includeCSS("../../html/custom_css.css")),
  
  # HTML header
  fluidRow(includeHTML("../../html/header.html")),
  
  # title ribbon
  fluidRow(
    style = "color:black; padding:0% 2% 0% 2%; text-align:left;",
    h3("Gene expression in skeletal muscle accross the human lifespan"),
    h5("Last update 2025-04-11"),
    tags$hr()
  ),
  
  # main page
  fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
           sidebarLayout(
             sidebarPanel(width = 3,
                          selectizeInput("inputGeneSymbol", "Gene Symbol:", choices=NULL, multiple=F, width=1000),
                          checkboxGroupInput("diagnosis_sarcopenia", 
                                             label = "Sarcopenia diagnosis", 
                                             selected = c("Healthy", "Sarcopenia"),
                                             choices = c("Healthy", "Sarcopenia")),
                          tags$b("Statistics"),
                          em(h5("Spearman correlation and Wilcoxon ranked signed test comparing ages to the middle group. The p values are not adjusted for multiple testing comparisons.")),
                          downloadButton("downloadGeneData", "Download Data")
                          
             ),
             mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                       plotOutput("genePlotAging", height="600px") %>% withSpinner(color="#5B768E")
             )
           ),
           
  ),
  
  # description of methods
  fluidRow(
    style = "color:black; background-color:white; padding:0% 2% 0% 2%;",
    tags$hr(),
    h3("Methods"),
    
    p("Publicly available transcriptomic datasets profiling human skeletal muscle were retrieved from the Gene Expression Omnibus (GEO), including both RNA sequencing (RNA-seq) and Affymetrix microarray platforms (Human Exon 1.0 ST and HGU133Plus2). Only datasets that reported individual age data, with a span of more than 15 years were retained. Metadata for each dataset, including age, sex, diagnosis, and platform, were manually curated and standardized across sources."),
    
    p("Raw count data from GTEx and multiple GEO datasets were aggregated by gene and annotated using the Homo.sapiens Bioconductor package. Gene expression filtering and normalization were performed using edgeR, followed by log2 transformation and median centering. For microarrays, raw CEL files were processed with the oligo package, normalized using RMA, filtered for low expression, annotated, and collapsed by gene symbol."),
    
    p("All expression matrices were merged by gene symbol and batch-corrected using limma. Principal component analysis was used to assess integration and sample clustering. Sex assignments were verified using expression of sex-specific markers (RPS4Y1 and XIST) and adjusted where discrepancies were observed."),
    
    p("Outliers were detected based on pairwise sample correlations, and samples below the 1st percentile of average correlation were excluded. Final metadata and normalized expression matrices were saved for visualization and analysis.")
    
    ),
  
  # Table with datasets
  fluidRow(style="color:black;background-color:white;padding:0% 2% 0% 2%;",
           tags$hr(),
           h3("Datasets Included in the Analysis"),
           dataTableOutput("datasets")
  ),
  
  # Citation
  fluidRow(style="color:black;background-color:white;padding:0% 2% 2% 2%;",
           tags$hr(),
           h3("Citation"),
           "Unpublished data."
  ),
  
  # HTML footer
  fluidRow(includeHTML("../../html/footer.html"))
)