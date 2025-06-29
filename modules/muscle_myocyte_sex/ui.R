#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle cells - UI
#
#--------------------------------------------------------------------------------------------------------

ui <- fluidPage(
  title = "MuscleMyocyteSex",         # Title of the browser tab
  style = "padding:0%",              # Remove default padding around the page
  
  #---------------------------------------------------------
  # Load custom CSS for consistent style across the app
  tags$head(includeCSS("../../www/style.css")),
  
  #---------------------------------------------------------
  # Navigation bar layout with multiple tabs
  navbarPage(
    # Custom title: logo image + text
    title = HTML('
      <div style="display: flex; align-items: center;margin: -10px;">
        <img src="../../www/img/snippet/muscle_myocyte_sex.png" style="height: 40px; margin-right: 10px;">
        <span style="font-size: 20px; font-weight: bold; margin-right: 10px;">Muscle myocyte sex</span>
      </div>
    '),
    
    #=====================================================================
    #============================ TAB 1 =================================
    #=====================================================================
    tabPanel("Analysis",  # Main tab for exploring data and results
             
             # Introductory text with link to the methods section
             p(HTML('Sex differences in the skeletal muscle cell transcriptome  \n To understand how the data was generated, read the 
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
             p("This app offers an accessible resource to explore how biologcial sex affect skeletal muscle cells. 
               By examining the effect of biological sex on the muscle cell transcriptome, this tool opens new 
         avenues for research, therapeutic discovery, and prevention strategies targeting metabolic disorders."),
             
             # Data processing and integration methodology
             h3("Methods"),
             p("Publicly available datasets were downloaded from the Gene Expression Omnibus (GEO), including both 
        microarray and RNA-seq platforms from human, mouse, and rat skeletal muscle cell models. Each dataset 
        was processed using platform-appropriate normalization methods — Robust Multi-array Average (RMA) for 
        microarrays and Variance Stabilizing Transformation (VST) for RNA-seq. Lowly expressed genes were filtered 
        out."),
             
             # Statistical tests
             h3("Statistics"),
             p("Statistical analysis is based on the Wilcoxon signed-rank test comparing females to males. 
        The results table presents both unadjusted and Bonferroni-adjusted p-values to correct for multiple testing across 
        all transcripts in the database, ensuring a highly conservative approach. Significance is indicated as follows: 
        * for FDR < 0.05, ** for FDR < 0.01, and *** for FDR < 0.001. 'ns' denotes non-significant results."),
             
             # Citation
             h3("Citation"),
             HTML('
  <div style="margin: 0rem;">
    <div>
      MacGregor K, Ellefsen S, Pillon NJ, Hammarström D Krook A.
      <strong>Sex differences in skeletal muscle metabolism in exercise and type 2 diabetes mellitus.</strong>
      <em>
Nat Rev Endocrinol.</em> 2025 Mar;21(3):166-179.
      <a href="https://doi.org/10.1038/s41574-024-01058-9" target="_blank">https://doi.org/10.1038/s41574-024-01058-9</a>
    </div>
    <div style="margin-top: 1rem;">
      <script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>
      <span class="__dimensions_badge_embed__" data-id="pub.1182757379" data-legend="always"></span>
    </div>
  </div>
'),
             
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
