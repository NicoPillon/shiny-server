#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic myotube response to electrical pulse stimulation
#
#--------------------------------------------------------------------------------------------------------

ui <- fluidPage(title="MyotubeEPS",
                style="padding:0%",
                
                # CSS for style
                tags$head(includeCSS("../../www/style.css")),

                # main page
                navbarPage("MyotubeEPS",

                           tabPanel(
                             title = "Plots",

                             selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                             actionButton("updatePlot", "Refresh plot", icon("refresh")),
                             downloadButton("downloadGeneData", "Download Data"),
                             tags$hr(),
                             sidebarLayout(
                               sidebarPanel(width = 3,
                                            checkboxGroupInput("cell_type", 
                                                               label = "Cell Type", 
                                                               selected = c("C2C12", "primary"),
                                                               choices = c("C2C12", "primary")),
                                            checkboxGroupInput("duration", 
                                                               label = "Duration (hours)", 
                                                               selected = c("3", "3+Rest3h", "6", "24", "48", "72", "72(1h+3rest)"),
                                                               choices = c("3", "3+Rest3h", "6", "24", "48", "72", "72(1h+3rest)")),
                                            checkboxGroupInput("pulse_param", 
                                                               label = "Pulse Parameters", 
                                                               selected = c("14V, 5Hz, 2ms", "3V, 1Hz, 6ms", "13V, 66Hz, 2ms",
                                                                            "10V, 1Hz, 2ms", "30V, 1Hz, 2ms", "10V, 1Hz, 10ms", 
                                                                            "11.5V, 1Hz, 2ms"),
                                                               choices = c("14V, 5Hz, 2ms", "3V, 1Hz, 6ms", "13V, 66Hz, 2ms",
                                                                           "10V, 1Hz, 2ms", "30V, 1Hz, 2ms", "10V, 1Hz, 10ms", 
                                                                           "11.5V, 1Hz, 2ms")),
                                            tags$hr(),
                                            tags$b("Statistics"),
                                            em("Wilcoxon ranked signed test comparing post to pre. The p values are not adjusted for multiple testing comparisons.")
                                            
                               ),
                               mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                                         plotOutput("geneBoxplot", height="700px") %>% withSpinner(color="#5B768E", type = 8))
                               )
                           ),

                           # description of methods
                           tabPanel("Description",
                                    h3("Citation"),
                                    HTML('Unpublished'),
                                    
                                    h3("Methods"),
                                    p("This application enables the integrated analysis and visualization of transcriptomic responses to electrical pulse stimulation in skeletal myotubes across multiple species and platforms."),
                                    p("Transcriptomic datasets from human, mouse, and rat skeletal myotubes exposed to electrical pulse stimulation (EPS) were compiled and integrated. Orthologous gene annotations were generated using biomaRt to harmonize gene identifiers across species. Raw data from RNA-sequencing and microarray platforms were normalized by median centering. The datasets were then merged based on ortholog mappings. Outlier samples were identified and excluded following quality control based on pairwise correlation matrices. Genes with consistently low expression or detected in fewer than 50% of datasets were filtered out. Batch effects across studies were corrected using the removeBatchEffect function from the limma package. Differential gene expression between EPS and control conditions was assessed using linear modeling with a paired design, followed by moderated t-tests and Bonferroni adjustment for multiple testing. Gene-level fold changes and sample-level pairwise correlations were computed for visualization and secondary quality control. Study-specific metadata were extracted from the GEO accession records and summarized."),
                                    tags$hr(),
                                    
                                    h3("Datasets"),
                                    tags$p(
                                      tags$b("Are we missing a relevant study? Please "),
                                      a("let us know!", href = "mailto:nicolas.pillon@ki.se", target = "_blank")
                                    ),    
                                    dataTableOutput("references")

                                    
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
)

                           
                