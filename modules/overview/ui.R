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
           column(8,
                  selectizeInput("inputGeneSymbol", "Gene Symbol:", 
                                 choices=NULL, multiple=F),
                  textOutput("gene_description")
           ),
           column(4,
                  r3dmolOutput("structureViewer", height = "200px"))
  ),
  tags$hr(),
  fluidRow(style="color:black;background-color:white;padding:0% 5% 1% 5%;",
           div(class = "col-sm-12 col-md-4",
               style = "padding-top: 10px; padding-bottom: 10px;",
               tags$a(style="font-weight: bold; font-size: 1.75rem;",
                      href = "../human_muscle_aging.html", # check if relative path works on server
                      tags$img(
                        style = "width: 40px; height: auto;",
                        src = "/www/img/snippet/human_muscle_aging.png"
                      ),
                      "Aging",
                      plotOutput("plot_human_muscle_aging", height = "300px") %>%
                        withSpinner(color = "#5B768E", type = 8)
               )
           ),
           div(class = "col-sm-12 col-md-4",
               style = "padding-top: 10px; padding-bottom: 10px;",
               tags$a(style="font-weight: bold; font-size: 1.75rem;",
                      href = "../human_muscle_obesity.html",
                      tags$img(
                        style = "width: 40px; height: auto;",
                        src = "/www/img/snippet/human_muscle_obesity.png" # check if absolute path works on server
                      ),
                      "Obesity",
                      plotOutput("plot_human_muscle_obesity", height = "300px") %>%
                        withSpinner(color = "#5B768E", type = 8)
               )
           ),
           div(class = "col-sm-12 col-md-3",
               style = "padding-top: 10px; padding-bottom: 10px;",
               tags$a(style="font-weight: bold; font-size: 1.75rem;",
                      href = "../fiber_types.html",
                      tags$img(
                        style = "width: 40px; height: auto;",
                        src = "/www/img/snippet/fiber_types.png"
                      ),
                      "Fiber Types",
                      plotOutput("plot_fiber_types", height = "300px") %>%
                        withSpinner(color = "#5B768E", type = 8)
               )
           )
  ),
  fluidRow(style="color:black;background-color:white;padding:0% 5% 1% 5%;",
           div(class = "col-sm-12 col-md-3",
               style = "padding-top: 10px; padding-bottom: 10px;",
               tags$a(style="font-weight: bold; font-size: 1.75rem;",
                      href = "../muscle_models.html",
                      tags$img(
                        style = "width: 40px; height: auto;",
                        src = "/www/img/snippet/muscle_models.png"
                      ),
                      "Muscle Models",
                      plotOutput("plot_muscle_models", height = "400px") %>%
                        withSpinner(color = "#5B768E", type = 8)
               )
           ),
           div(class = "col-sm-12 col-md-3",
               style = "padding-top: 10px; padding-bottom: 10px;",
               tags$a(style="font-weight: bold; font-size: 1.75rem;",
                      href = "../MetaMEx.html",
                      tags$img(
                        style = "width: 40px; height: auto;",
                        src = "/www/img/snippet/MetaMEx.png"
                      ),
                      "MetaMEx",
                      plotOutput("plot_metamex", height = "400px") %>%
                        withSpinner(color = "#5B768E", type = 8)
               )
           ),
           tags$hr()
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
