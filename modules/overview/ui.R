#--------------------------------------------------------------------------------------------------------
#
# OVERVIEW - UI
#
#--------------------------------------------------------------------------------------------------------

ui <- fluidPage(
  title = "Overview",         # Title of the browser tab
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
        <img src="../../www/img/snippet/overview.png" style="height: 40px; margin-right: 10px;">
        <span style="font-size: 20px; font-weight: bold; margin-right: 10px;">Overview</span>
      </div>
    '),
    
    #=====================================================================
    #============================ TAB 1 =================================
    #=====================================================================
    tabPanel("Analysis",  # Main tab for exploring data and results
             
             style = "padding:0% 2% 1% 2%;",
             
             # Introductory text with link to the methods section
             p(HTML('Overview of the main experimental responses across all datasets. \n To understand how the data was generated, read the 
             <a href="#" onclick="$(\'.navbar-nav a:contains(\\\'Methods\\\')\').click()">methods</a>.')),

             # Gene selection input
             selectizeInput("inputGeneSymbol", 
                            label = "Choose a target:",
                            choices = character(0),
                            multiple = FALSE,
                            width = "50%"),

             # Main plot
             highchartOutput("overviewBarHighcharter", height = "400px") %>% withSpinner(color = "#5B768E", type = 8),
             
             # Legend
             tags$em(p("The plot presents the log2-transformed fold-change compared to the control group for each dataset. Values above 0 (log2FC > 0) indicate increased expression, while values below 0 (log2FC < 0) indicate decreased expression. Green bars indicate statistical significance at FDR < 0.05.")),
             
             # Spacer
             br(),
             
             # Gene description
             tags$hr(),
             tags$h4("Target Description"),
             uiOutput("gene_description")
             ),
    
    
    #=====================================================================
    #============================ TAB 2 =================================
    #=====================================================================
    tabPanel("Methods",
             
             style = "padding:0% 5% 1% 5%;",
             
             h3("Why this app?"),
             p("The portal MuscleOmics.org compiles a growing number of transcriptomic and proteomic datasets related to skeletal muscle physiology and metabolism. This app provides a unified interface to explore differential expression results across multiple populations, models, and experimental conditions."),
             
             h3("Methods"),
             p("All datasets were downloaded from public repositories or derived from published studies. Each dataset was analyzed individually using limma-based statistics or equivalent models appropriate for the measurement platform. When applicable, lowly expressed genes were filtered out, and normalization and batch correction were applied. For detailed information on a specific dataset, please refer to the corresponding method section within its dedicated module."),
             
             h3("Data Presentation"),
             p("The results include differential expression statistics computed from normalized expression matrices. All relevant populations were retained to provide a comprehensive overview of each dataset. Each comparison displays the log2 fold-change along with adjusted p-values computed using the Benjamini–Hochberg method. Significance is color-coded: green indicates FDR < 0.05, and grey indicates non-significant results."),
             
             h3("Citation"),
             p("Pillon NJ, Rizo Roca D, Macgregor K. Unpublished. Please contact the authors before reuse.")
    )
  )
)
