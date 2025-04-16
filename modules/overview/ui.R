#-------------------------------------------------------------------------------------------
#
# Overview
#
#-------------------------------------------------------------------------------------------

ui <- fluidPage(
  # code staggering computation
  useShinyjs(),
  tags$script(HTML("
  Shiny.addCustomMessageHandler('plotRendered', function(module) {
    Shiny.setInputValue(module + '_done', Math.random());
  });
")),
  
  # Google analytics
  tags$head(includeScript("../../www/google-analytics.html")),
  
  # CSS for style
  tags$head(includeCSS("../../www/style.css")),
  
  # title ribbon
  fluidRow(
    style = "color:black; padding:0% 2% 0% 2%; text-align:left;",
    h3("Skeletal Muscle Database - Overview"),
    tags$hr()
  ),
  
  # main page
  fluidRow(style="color:black;background-color:white;padding:0% 5% 1% 5%;",
           selectizeInput("inputGeneSymbol", "Gene Symbol:", choices=NULL, multiple=F, width=600),
           uiOutput("gene_warning"),
           textOutput("gene_description"),
           tags$hr(),
           fluidRow(
             style = "color:black;background-color:white;",
             div(class = "col-sm-12 col-md-4",
                 tags$a(style="font-weight: bold;",
                   href = "../human_muscle_aging.html",
                   target = "_blank",  # optional: opens in a new tab
                   "Muscle Aging Module"
                 ),
                 plotOutput("plot_human_muscle_aging", height = "300px") %>%
                   withSpinner(color = "#5B768E", type = 8)
             ),
             div(class = "col-sm-12 col-md-4",
                 tags$a(style="font-weight: bold;",
                   href = "../human_muscle_obesity.html",
                   target = "_blank",  # optional: opens in a new tab
                   "Muscle Obesity Module"
                   ),
                 plotOutput("plot_human_muscle_obesity", height = "300px") %>%
                   withSpinner(color = "#5B768E", type = 8)
             ),
             div(class = "col-sm-12 col-md-3",
                 tags$a(style="font-weight: bold;",
                        href = "../fiber_types.html",
                        target = "_blank",  # optional: opens in a new tab
                        "Muscle Fiber Types Module"
                 ),
                 plotOutput("plot_fiber_types", height = "300px") %>%
                   withSpinner(color = "#5B768E", type = 8)
             ),
             tags$hr()
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
