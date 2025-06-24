ui <- fluidPage(
  title = "Overview",
  style = "padding:0%",  
  
  # Load CSS
  tags$head(includeCSS("../../www/style.css")),
  
  # Page container with padding
  fluidRow(
    style = "color:black;background-color:white;padding:2% 10%;",
    
    # Header section
    column(12,
           tags$h2("Overview"),
           tags$p(style = "font-size:16px;margin-top:-10px;",
                  "Select a target and explore datasets of interest by clicking the bars below.")
    ),
    
    # Gene selection input
    column(12,
           selectizeInput("inputGeneSymbol", 
                          label = "Choose a target:",
                          choices = NULL,
                          multiple = FALSE,
                          width = "50%")
    ),
    
    # Spacer
    br(),
    
    # Chart output
    column(12,
           highchartOutput("overviewBarHighcharter", height = "400px") %>% 
             withSpinner(color = "#5B768E", type = 8)
    ),
    
    # Gene description
    column(12,
           tags$hr(),
           tags$h4("Target Description"),
           uiOutput("gene_description")
    )
  ),
  
  # iFrame resize script
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
