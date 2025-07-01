#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic muscle response to hyperinsulinemic euglycemic clamp - UI
#
#--------------------------------------------------------------------------------------------------------

ui <- fluidPage(
  title = "MuscleClamp",         # Title of the browser tab
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
        <img src="../../www/img/snippet/muscle_clamp.png" style="height: 40px; margin-right: 10px;">
        <span style="font-size: 20px; font-weight: bold; margin-right: 10px;">Glucose Clamp</span>
      </div>
    '),
    
    #=====================================================================
    #============================ TAB 1 =================================
    #=====================================================================
    tabPanel("Analysis",  # Main tab for exploring data and results
             
             # Introductory text with link to the methods section
             p(HTML('Skeletal muscle tissue response to hyperinsulinemic euglycemic clamp using integrated transcriptomic data. \n To understand how the data was generated, read the 
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
                            sliderTextInput("duration", "Clamp Duration (min)",
                                            choices = c("30", "120", "180", "240"),
                                            selected = c("30", "240"),
                                            grid = TRUE, force_edges = TRUE),
                            checkboxGroupInput("diagnosis_diabetes", 
                                               label = "Diagnosis", 
                                               selected = c("Healthy", "Prediabetes", "T2D"),
                                               choices = c("Healthy", "Prediabetes", "T2D")),
                            checkboxGroupInput("treatment", 
                                               label = "Treatment", 
                                               selected = c("Control", "Thiazolidinedione"),
                                               choices = c("Control", "Thiazolidinedione")),
                            
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
             p("This app offers an accessible resource to explore how skeleletal muscle tissue responds during an hyperinsulinemic euglycemic clamp. By examining its impact on the muscle cell transcriptome, this tool opens new avenues for research, therapeutic discovery, and prevention strategies targeting metabolic disorders."),
             
             # Data processing and integration methodology
             h3("Methods"),
             p("This dataset enables the integrative analysis of transcriptional responses to insulin in human skeletal muscle across multiple independent hyperinsulinemic euglycemic clamp studies."),
             p("Publicly available gene expression datasets from hyperinsulinemic-euglycemic clamp studies in human skeletal muscle (Affymetrix microarrays and RNA-seq) were retrieved and processed. Gene symbols were harmonized across platforms, and only genes present in ≥2 datasets were retained. Expression values were median-centered per dataset and batch effects were adjusted using removeBatchEffect (limma). Outlier samples were identified and excluded based on PCA, pairwise Spearman correlations, and correlation rank thresholds. Sex assignment was verified using expression of sex-linked genes (XIST, EIF1AY, RPS4Y1, KDM5D)."),
             tags$hr(),
             
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
    )
  )
)
