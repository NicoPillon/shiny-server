#-----------------------------------------------------------------------
#
# Profile of skeletal muscle fibers
#
#----------------------------------------------------------------------

ui <- fluidPage(
  # CSS for style
  tags$head(includeCSS("../../www/style.css")),
  
  # title ribbon
  fluidRow(
    style = "color:black; padding:0% 2% 0% 2%; text-align:left;",
    h3("Proteome and Transcriptome of Skeletal Muscle Fiber Types"),
    h5("Last update 2025-04-13"),
    tags$hr()
  ),
  
  # main page
  fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
           sidebarLayout(
             sidebarPanel(width = 3,
                          selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                          actionButton("updatePlot", "Refresh plot", icon("refresh")),
                          tags$hr(),
                          tags$b("Statistics:"),
                          em("Kruskal–Wallis tests were performed to assess whether at least one fiber type differs significantly from the others. Reported p-values are unadjusted for multiple comparisons."),
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
  fluidRow(
    style = "color:black; background-color:white; padding:0% 2% 0% 2%;",
    tags$hr(),
    h3("Methods"),
    
    p("This analysis integrates proteomic and transcriptomic datasets of human skeletal muscle fibers to identify and compare fiber-type specific signatures."),
    
    p("For proteomics, raw intensity values were retrieved from supplementary Excel files, and metadata were cleaned and harmonized. Protein expression matrices were constructed by aligning sample identifiers and gene symbols. Rows with multiple gene mappings were split and only those with fewer missing values were retained. Fiber types were annotated based on the relative abundance of MYH isoforms, using published criteria: Type I (MYH7 ≥ 80%), Type IIA (MYH2 ≥ 80%), and Type IIX (MYH1 ≥ 20%)."),
    
    p("For transcriptomics, raw UMI counts or CEL files were preprocessed with standard Bioconductor pipelines. Lowly expressed genes were filtered, normalized using TMM or RMA, and batch effects were removed with `limma::removeBatchEffect`. Gene identifiers were mapped to gene symbols. Metadata were extracted and harmonized to enable sample-wise comparison. Single muscle fibers were annotated based on the expression levels of MYH7, MYH1, and MYH2. Normalized expression values were back-transformed to linear scale and used to compute isoform contributions. Fibers were classified into 'Type I', 'Type IIA', or 'Type IIX' based on tertiles of isoform expression, with intermediate profiles labeled as 'Mixed' and excluded.")
    
  ),
  
  
  # Table with datasets
  fluidRow(style="color:black;background-color:white;padding:0% 2% 1% 2%;",
           h3("Datasets Included in the Analysis"),
           #          "Transcriptomics from manually isolated skeletal muscle fibers from", 
           #          a("Raue et al, 2012 (GSE28392)", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28392", 
           #            target="_blank", style="color:#5B768E"), 
           #          "and",
           #          a("Rubenstein et al, 2020 (GSE130977)", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130977", 
           #            target="_blank", style="color:#5B768E"),
           #          "for a total of 88 fibres.",
           #          tags$br(),
           #          "Proteomics from manually isolated skeletal muscle fibers from", 
           #          a("Murgia et al, 2017", href="https://doi.org/10.1016/j.celrep.2017.05.054", 
           #            target="_blank", style="color:#5B768E"), "and", 
           #          a("Murgia et al, 2022", href="https://doi.org/10.1093/pnasnexus/pgac086", 
           #            target="_blank", style="color:#5B768E"),
           #          "for a total of 387 fibres.",
           tags$p(
             tags$b("Are we missing a relevant study? Please "),
             a("let us know!", href = "mailto:nicolas.pillon@ki.se", target = "_blank")
             )
           ),
           
           # Citation
           fluidRow(
             style = "color:black; background-color:white; padding:0% 2% 2% 2%;",
             tags$hr(),
             h3("Citation"),
             tags$ul(
               tags$li(
                 "Raue U, Trappe TA, Estrem ST, Qian H, Helvering LM, Smith RC, Trappe S. ",
                 a("Transcriptome signature of resistance exercise adaptations: mixed muscle and fiber type specific profiles in young and old adults.", 
                   href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28392", 
                   target = "_blank", style = "color:#5B768E"), 
                 " J Appl Physiol (1985) 2012 May;112(10):1625-36."
               ),
               tags$li(
                 "Rubenstein AB, Smith GR, Raue U, Begue G et al. ",
                 a("Single-cell transcriptional profiles in human skeletal muscle.", 
                   href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130977", 
                   target = "_blank", style = "color:#5B768E"), 
                 " Sci Rep 2020 Jan 14;10(1):229."
               ),
               tags$li(
                 "Murgia M, Toniolo L, Nagaraj N, Ciciliot S, Vindigni V, Schiaffino S, Reggiani C, Mann M. ",
                 a("Single Muscle Fiber Proteomics Reveals Fiber-Type-Specific Features of Human Muscle Aging", 
                   href = "https://doi.org/10.1016/j.celrep.2017.05.054", 
                   target = "_blank", style = "color:#5B768E"), 
                 " Cell Rep. 2017 Jun 13;19(11):2396-2409."
               ),
               tags$li(
                 "Murgia M, Ciciliot S, Nagaraj N, Reggiani C, Schiaffino S, Franchi MV, Pišot R, Šimunič B, Toniolo L, Blaauw B, Sandri M, Biolo G, Flück M, Narici MV, Mann M.",
                 a("Signatures of muscle disuse in spaceflight and bed rest revealed by single muscle fiber proteomics", 
                   href = "https://doi.org/10.1093/pnasnexus/pgac086", 
                   target = "_blank", style = "color:#5B768E"), 
                 " PNAS Nexus. 2022 Jun 11;1(3):pgac086."
               )
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
  