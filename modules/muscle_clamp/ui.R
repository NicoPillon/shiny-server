#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic muscle response to hyperinsulinemic euglycemic clamp
#
#--------------------------------------------------------------------------------------------------------

ui <- fluidPage(title="MuscleClamp",
                style="padding:0%",
                
                # CSS for style
                tags$head(includeCSS("../../www/style.css")),

                # main page
                navbarPage("MuscleClamp",

                           tabPanel(
                             title = "Plots",

                             selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                             actionButton("updatePlot", "Refresh plot", icon("refresh")),
                             downloadButton("downloadGeneData", "Download Data"),
                             tags$hr(),
                             sidebarLayout(
                               sidebarPanel(width = 3,
                                            sliderInput("duration", tags$b("Clamp Duration (min)"),
                                                        min = 30, max = 240, value = c(30,240), step = 30, sep = ""),
                                            checkboxGroupInput("diagnosis_diabetes", 
                                                               label = "Diagnosis", 
                                                               selected = c("Healthy", "Prediabetes", "T2D"),
                                                               choices = c("Healthy", "Prediabetes", "T2D")),
                                            checkboxGroupInput("treatment", 
                                                               label = "Treatment", 
                                                               selected = c("Control", "Thiazolidinedione"),
                                                               choices = c("Control", "Thiazolidinedione")),
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
                                    p("This dataset enables the integrative analysis of transcriptional responses to insulin in human skeletal muscle across multiple independent hyperinsulinemic euglycemic clamp studies."),
                                    p("Publicly available gene expression datasets from hyperinsulinemic-euglycemic clamp studies in human skeletal muscle (Affymetrix microarrays and RNA-seq) were retrieved and processed. Gene symbols were harmonized across platforms, and only genes present in ≥2 datasets were retained. Expression values were median-centered per dataset and batch effects were adjusted using removeBatchEffect (limma). Outlier samples were identified and excluded based on PCA, pairwise Spearman correlations, and correlation rank thresholds. Sex assignment was verified using expression of sex-linked genes (XIST, EIF1AY, RPS4Y1, KDM5D)."),
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

                           
                