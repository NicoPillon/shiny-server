#--------------------------------------------------------------------------------------------------------
#
# App short description
#
#--------------------------------------------------------------------------------------------------------

ui <- fluidPage(
  # Google analytics
  tags$head(includeScript("../../www/google-analytics.html")),
  
  # CSS for style
  tags$head(includeCSS("../../www/style.css")),

  # title ribbon
  fluidRow(
    style = "color:black; padding:0% 2% 0% 2%; text-align:left;",
    h3("Cross-Species Analysis of Skeletal Muscle Models Reveals Conserved and Divergent Transcriptional Programs"),
    h5("Last update 1978-03-28"),
    tags$hr()
  ),
  
  # main page
  fluidRow(style="color:black;background-color:white;padding:0% 5% 1% 5%;",
           sidebarLayout(
             sidebarPanel(width = 3, 
                          selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                          actionButton("updatePlot", "Refresh plot", icon("refresh")),
                          tags$hr(),
                          tags$b("Statistics:"),
                          em("Spearman correlation and Wilcoxon ranked signed test comparing ages to the middle group. The p values are not adjusted for multiple testing comparisons."),
                          tags$hr(),
                          downloadButton("downloadGeneData", "Download Data")
             ),
             mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                       plotOutput("geneBoxplot", height="700px") %>% withSpinner(color="#5B768E", type = 8)
                       )
             )
           ),
  
  # description of methods
  fluidRow(
    style = "color:black; background-color:white; padding:0% 2% 0% 2%;",
    tags$hr(),
    h3("Methods"),
    p("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Integer nec odio. Praesent libero. Sed cursus ante dapibus diam. Sed nisi."),
    
    p("Nulla quis sem at nibh elementum imperdiet. Duis sagittis ipsum. Praesent mauris. Fusce nec tellus sed augue semper porta. Mauris massa."),
    
    p("Vestibulum lacinia arcu eget nulla. Class aptent taciti sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos."),
    
    p("Curabitur sodales ligula in libero. Sed dignissim lacinia nunc. Curabitur tortor. Pellentesque nibh. Aenean quam."),
    
    p("In scelerisque sem at dolor. Maecenas mattis. Sed convallis tristique sem. Proin ut ligula vel nunc egestas porttitor. Morbi lectus risus, iaculis vel, suscipit quis, luctus non, massa.")
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
           "Doe JA, Nguyen TM, Patel RK, Svensson L, Okafor U, Nakamura Y, Rossi A, Zhang H, López F, Larsson M.",
           a("Cross-species transcriptomic integration reveals conserved metabolic signatures in skeletal muscle.",
             href="https://doi.org/10.1234/fictjournal.01234.2025", target="_blank", style="color:#5B768E"),
           "Fictional Journal of Integrative Biology. 2025 Apr;12(2):101-117."
  )
  
)
