#--------------------------------------------------------------------------------------------------------
#
# App short description
#
#--------------------------------------------------------------------------------------------------------

ui <- fluidPage(
  # CSS for style
  tags$head(includeCSS("../../www/style.css")),
  
  # title ribbon
  fluidRow(
    style = "color:black; padding:0% 2% 0% 2%; text-align:left;",
    h3("This could be your own app hosted on the portal"),
    h5("Last update XXXX-YY-ZZ"),
    tags$hr()
  ),
  
  navbarPage("",
             tabPanel(
               title = "Transcriptomics",
               
               # Main layout with plot and controls
               fluidRow(
                 style = "color:black;background-color:white;padding:0% 5% 1% 5%;",
                 sidebarLayout(
                   sidebarPanel(
                     width = 3, 
                     selectizeInput("inputGeneSymbol", "Gene Symbols:", choices = NULL, multiple = TRUE, width = 600),
                     actionButton("updatePlot", "Refresh plot", icon("refresh")),
                     tags$hr(),
                     tags$b("Statistics:"),
                     em("Describe briefly what statistics were used in your plots."),
                     tags$hr(),
                     downloadButton("downloadGeneData", "Download Data")
                   ),
                   mainPanel(
                     width = 9,
                     style = "padding:0% 4% 1% 4%;",
                     plotOutput("geneBoxplot", height = "700px") %>% 
                       withSpinner(color = "#5B768E", type = 8)
                   )
                 )
               ),
               
               # Description of methods
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
               
               # Table with references
               fluidRow(
                 style = "color:black;background-color:white;padding:0% 2% 0% 2%;",
                 tags$hr(),
                 h3("Datasets Included in the Analysis"),
                 dataTableOutput("references"),
                 tags$p(
                   tags$b("Are we missing a relevant study? Please "),
                   a("let us know!", href = "mailto:nicolas.pillon@ki.se", target = "_blank")
                 )
               ),
               
               # Citation
               fluidRow(
                 style = "color:black;background-color:white;padding:0% 2% 2% 2%;",
                 tags$hr(),
                 h3("Citation"),
                 p(
                   "Doe JA, Nguyen TM, Patel RK, Svensson L, Okafor U, Nakamura Y, Rossi A, Zhang H, López F, Larsson M. ",
                   a("Your publication showing how the data was used.",
                     href = "https://doi.org/10.1234/fictjournal.01234.2025", 
                     target = "_blank", style = "color:#5B768E"),
                   " Fictional Journal of Skeletal Biology. 2025 Apr;12(2):101-117."
                 )
               )
             ),
             
             # Other tabs
             tabPanel(
               title = "Metabolomics",
               fluidRow(
                 style = "color:black;background-color:white;padding:0% 5% 1% 5%;",
                 "Coming soon..."
               )
             ),
             
             tabPanel(
               title = "Proteomics",
               fluidRow(
                 style = "color:black;background-color:white;padding:0% 5% 1% 5%;",
                 "Coming soon..."
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