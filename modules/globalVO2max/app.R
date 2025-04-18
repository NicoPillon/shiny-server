#-----------------------------------------------------------------------
#
# Reference values for VO2 max
#
#----------------------------------------------------------------------
# Load data and libraries
library(shinycssloaders)
library(ggplot2)
library(ggrepel)
library(DT)
library(cowplot)
library(rlang)
library(dplyr)

#--------------------------------------------------------------------------------------------------------
# Palmitate
#--------------------------------------------------------------------------------------------------------
VO2_data <- readRDS("data/data.Rds")

# lsit of studies
studies_included <- readRDS("data/studies_included.Rds")
studies_included$PMID <- as.character(studies_included$PMID)
studies_included <- studies_included[order(studies_included$First.Author),]
sort(unique(studies_included$Country))

# function for plots
fct_nomogram <- function(dat){
  ggplot(dat, aes(x = Age, y = VO2max,
                       color = Percentile, linetype = Percentile)) +
    theme_bw(16) +
    theme(panel.grid.major = element_line(color = "gray80",
                                          linewidth = 0.5,
                                          linetype = 1),
          legend.title = element_blank()) +
    scale_x_continuous(breaks = seq(0, 100, 10),
                       minor_breaks = seq(0, 100, 1)) +
    scale_y_continuous(breaks = seq(0, 100, 10),
                       minor_breaks = seq(0, 100, 1)) +
    labs(x = "Age (Years)",
         y = expression(dot("V")["O"[2]] ~ "peak (mL/min/kg)")) +
    #facet_wrap(.~Sex, scales = "free_y", ncol = 2) +
    geom_smooth(se = FALSE, method = "lm", linewidth = 0.5, alpha = 0.9) +
    scale_colour_manual(values = c("#097910", "#1d7c0f", "#38800d", "#53840c",
                                   "black",
                                   "#8c8d08", "#a99106", "#c99604", "#e69a02")) +
    scale_linetype_manual(values = c(5,4,3,2,1,2,3,4,5))
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
                                style = "height:8vh; display:flex; justify-content:center; align-items:center;",
                                tags$a(href = "https://shiny.nicopillon.com", 
                                       icon("home", class = "fa-2x"), 
                                       style = "color:white; text-decoration:none;")  # Ensuring icon is white and no underline
                         ),
                         column(10,
                                h3("Worldwide reference values for maximal oxygen uptake"),
                                h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                           target="_blank", style="color:#D9DADB"), 
                                   "/ last update 2024-11-22")
                                ),
                         column(1,
                         )
                         ),
                
                # main page
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         sidebarLayout(
                           sidebarPanel(width = 3,
                               selectInput("sex",
                                           label = "Sex", 
                                           selected = "Female", 
                                           choices = c("Female", "Male")),
                               numericInput("age",
                                            label = "Age (years)",
                                            value = 32),
                               numericInput("vo2max",
                                            label = "VO2max (mL/min/kg)",
                                            value = 33),
                               checkboxGroupInput("modality",
                                                  label = "Modality", 
                                                  selected = c("Cycle", "Treadmill"), 
                                                  choices = c("Cycle", "Treadmill")),
                               checkboxGroupInput("country",
                                           label = "Country", 
                                           selected = c("Brazil", "Canada", "China", "Czechia", "Denmark", "Germany", 
                                                        "Greece", "Japan", "Lithuania", "Netherlands", "Norway", "Spain", 
                                                        "Sweden", "Switzerland", "United Kingdom", "United States"), 
                                           choices = c("Brazil", "Canada", "China", "Czechia", "Denmark", "Germany", 
                                                       "Greece", "Japan", "Lithuania", "Netherlands", "Norway", "Spain", 
                                                       "Sweden", "Switzerland", "United Kingdom", "United States")),
                               sliderInput("date", "Publication date",
                                           min = 1990, max = 2023, value = c(1990,2023), step = 1, sep = "")
                               ),
                  mainPanel(width = 9,
                            plotOutput("VO2plotMale", height = "370px", width = "500px"),
                            plotOutput("VO2plotFemale", height = "370px", width = "500px")
                            )
                  )
                  ),

                tags$hr(),
                
                # Table with datasets
                fluidRow(style="color:black;background-color:white;padding:0% 2% 1% 2%;",
                         h3("Datasets Selected in the Analysis"),
                         DT::dataTableOutput("studies_included")
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
  

  plotData <- reactive({
    
    # input values
    input_age = 30
    input_vo2max = 40
    input_sex = "Male"
    input_modality = c("Cycle", "Treadmill")
    input_country = c("Brazil", "Canada", "China", "Czechia", "Denmark", "Germany", "Greece", "Japan", "Lithuania", "Netherlands", "Norway",
                      "Spain", "Switzerland", "United Kingdom", "United States")
    input_date_min = 2000
    input_date_max = 2022
    
    input_modality = input$modality
    input_country = input$country
    input_date_min = input$date[1]
    input_date_max = input$date[2]
    
    # error messages
    validate(need(input_age, "Please provide a numeric value age"))
    validate(need(input_vo2max, "Please provide a numeric value for VO2max"))

    # subset modality
    VO2_data <- subset(VO2_data, Modality %in% input_modality)
    
    # subset country
    VO2_data <- subset(VO2_data, Country %in% input_country)
    
    # subset publication year
    VO2_data <- subset(VO2_data, Publication.year %in% seq(input_date_min, input_date_max))
    VO2_data
    })
  
    
    output$VO2plotMale <- renderPlot({

    # subset selected sex
    VO2_data_M <- subset(plotData(), Sex == "Male")

    # male plot
    nomogram_M <- fct_nomogram(VO2_data_M) +
      labs(subtitle = "Male") +
      theme(legend.position = "right")

    nomogram_M <- if(input$sex == "Male")
      nomogram_M +
      annotate("label",
               x = input$age,
               y = input$vo2max -2,
               label = "Your VO2max",
               color = "red") +
      annotate("pointrange",
               x = input$age,
               y = input$vo2max,
               ymin = input$vo2max, ymax = input$vo2max,
               color = "red",
               size = 0.5) else nomogram_M

    nomogram_M
   })
  
    output$VO2plotFemale <- renderPlot({
      VO2_data_F <- subset(plotData(), Sex == "Female")

      # female plot
      nomogram_F <- fct_nomogram(VO2_data_F) +
        labs(subtitle = "Female") +
        theme(legend.position = "right")
      
      nomogram_F <- if(input$sex == "Female")
        nomogram_F +
        annotate("label",
                 x = input$age,
                 y = input$vo2max -2,
                 label = "Your VO2max",
                 color = "red") +
        annotate("pointrange",
                 x = input$age,
                 y = input$vo2max,
                 ymin = input$vo2max, ymax = input$vo2max,
                 color = "red",
                 size = 0.5) else nomogram_F
      
      nomogram_F

    })
    
    
    
      output$studies_included <- DT::renderDataTable(escape = FALSE, 
                                                 rownames = FALSE, 
                                                 options=list(paging = FALSE,
                                                              dom = 't'), 
                                                 {
                                                   selecTable <- plotData()
                                                   selecTable <- selecTable$PMID
                                                   selecTable <- subset(studies_included, PMID %in% selecTable)
                                                   })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)