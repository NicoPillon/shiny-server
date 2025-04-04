#-----------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle response to hypoxia
#
#----------------------------------------------------------------------
# Load data and libraries
library(feather)
library(shinycssloaders)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(cowplot)
library(DT)

#--------------------------------------------------------------------------------------------------------
# AMPK
#--------------------------------------------------------------------------------------------------------
AMPK_datamatrix <- readRDS("data/AMPK_datamatrix.Rds")
AMPK_metadata <- readRDS("data/AMPK_metadata.Rds")


#--------------------------------------------------------------------------------------------------------
# List of genes
#--------------------------------------------------------------------------------------------------------
genelist <- unique(rownames(AMPK_datamatrix))

#function to gene name
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                
                # Google analytics
                tags$head(includeScript("google-analytics.html")),
                
                # Custom CSS to change checkbox tick color
                tags$style(HTML("
                  input[type='checkbox'] {
                    accent-color: #c93f1e; /* Change the checkbox tick color here */
                  }
                ")),
                
                # title ribbon
                fluidRow(style="color:white;background-color:#5B768E;padding:0% 1% 1% 1%;text-align:center; display:flex; justify-content:center; align-items:center;",
                         column(1, 
                                tags$a(href = "https://shiny.nicopillon.com", 
                                       icon("home", class = "fa-2x"), 
                                       style = "color:white; text-decoration:none;")
                         ),
                         column(10,
                                h3("Transcriptomic response of skeletal muscle to AMPK modulation"),
                                h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                           target="_blank", style="color:#D9DADB"), 
                                   "/ last update 2025-04-04")
                         ),
                         column(1,
                         )
                ),
                
                # main page
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         sidebarLayout(
                           sidebarPanel(width = 3,
                                        selectizeInput("inputGeneSymbol", "Gene Symbol:", choices=NULL, multiple=F, width=1000)
                           ),
                           mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                                     plotOutput("ampkPlot", height="600px") %>% withSpinner(color="#5b768e")
                           )
                         ),
                         
                ),
                
                tags$hr(),
                
                # Table with datasets
                fluidRow(style="color:black;background-color:white;padding:0% 2% 1% 2%;",
                         h3("Datasets Included in the Analysis"),
                         #dataTableOutput("datasets")
                ),

                
                # Author section at the bottom
                fluidRow(style="color:white;background-color:#5B768E;padding:2% 1% 2% 1%;display: flex; align-items: top; ",
                         # column(2, align="right", 
                         #        tags$img(src = "https://ki.se/profile-image/nicpil", height = "120px", width = "120px")  # Insert image here
                         # ),
                         column(4, align="left", 
                                tags$b("About the author:"), tags$br(),
                                "Nicolas J. Pillon, PhD", tags$br(),
                                "Associate Professor, Karolinska Institutet", tags$br(),
                                icon("globe"), a("/inflammation-and-metabolism", href="https://ki.se/en/research/research-areas-centres-and-networks/research-groups/inflammation-and-metabolism-nicolas-pillons-research-group",
                                                 target="_blank", style="color:white"), tags$br(),
                                icon("linkedin"), a("/nicopillon", href="https://www.linkedin.com/in/nicopillon/",
                                                    target="_blank", style="color:white"), tags$br(),
                                tags$br(),
                                "Feel free to write to me with feedback or questions:", tags$br(),
                                icon("envelope"), a("nicolas.pillon@ki.se", href="mailto:nicolas.pillon@ki.se",
                                                    target="_blank", style="color:white"), tags$br(),
                                
                         ),
                         column(4, align="center",
                                #tags$b("© 2024 Nicolas Pillon"), tags$br(),
                                
                         ),
                         column(4, align="right",
                                tags$b("Disclaimer:"), tags$br(),
                                em("The authors disclaim any responsibility for the use or interpretation of the data 
                                   presented in this application. Users are solely responsible for ensuring the appropriate 
                                   use of any data they choose to re-use."), tags$br(),
                                tags$br(),
                                tags$b("© 2024 Nicolas Pillon"),
                         ),
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected=c("Cd36"), 
                       options=NULL)
  
  # AMPK plot
  output$ampkPlot <- renderPlot({
    validate(need(input$inputGeneSymbol, " "))
    genename <- c("Cd36")
    genename <- input$inputGeneSymbol
    
    plotdata <- data.frame(AMPK_metadata,
                           data = as.numeric(AMPK_datamatrix[genename,]))
    
    plotdata <- plotdata %>%
      group_by(geo_dataset, Strain) %>%
      mutate(control_mean = median(data[grepl("control", genotype_treatment)], na.rm = TRUE),
             data_normalized = data - control_mean) %>%
      ungroup()
    
    # order groups
    plotdata$genotype_treatment <- factor(plotdata$genotype_treatment,
                                          levels = c("control", "KOα1", "KOα2", "KOα1α2", 
                                                     "KOγ3", "TGγ3",
                                                     "LA2, 3mpk, 6h", "SA2, 3mpk, 6h",
                                                     "LA2, 30mpk, 6h", "SA2, 20mpk, 6h",
                                                     "AICAR, 250 mg/kg/day, 6 days", 
                                                     "AICAR, 500 mg/kg/day, 3 days",
                                                     "AICAR, 500 mg/kg/day, 7 days",
                                                     "AICAR, 500 mg/kg/day, 14 days",
                                                     ""))
    # plot
    ggplot(plotdata, aes(x = genotype_treatment, y = data_normalized)) +
      theme_bw(16) +
      theme(legend.position = "none") +
      labs(y = paste0(genename, " mRNA (relative to control)")) +
      geom_boxplot(size = 0.2, fill = "gray90") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.title.x = element_blank()) +
      geom_point(aes(colour = geo_dataset), size = 2) +
      facet_wrap(.~Experiment, scale = "free_x", ncol = 9)
    
  })
  
  ##################################################################################################################
  #Dataset tables
  output$datasets <- renderDataTable(options=list(signif = 3),{
    DT::datatable(
      human_references, 
      escape = FALSE, 
      rownames = FALSE,
      options = list(
        columnDefs = list(
          list(className = 'dt-center', targets = 1),  # Align column to the center
          list(className = 'dt-center', targets = 2)  # Align column to the center
          # Add more lines for additional columns if needed
        )
      )
    )
  })
  
}


# Run the app ----
shinyApp(ui = ui, server = server)