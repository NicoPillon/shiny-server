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

#--------------------------------------------------------------------------------------------------------
# Palmitate
#--------------------------------------------------------------------------------------------------------
VO2_data <- readRDS("data/data.Rds")

# make factor to place mean in between percentiles
VO2_data$Percentile <- factor(VO2_data$Percentile,
                              levels = c("Percentile 0.9", "Percentile 0.75", "Percentile 0.6",
                                         "Median", 
                                         "Percentile 0.4", "Percentile 0.25", "Percentile 0.1"))
# lsit of studies
studies_included <- readRDS("data/studies_included.Rds")
studies_excluded <- readRDS("data/studies_excluded.Rds")

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
         y = "VO2max (mL/min/kg)") +
    #facet_wrap(.~Sex, scales = "free_y", ncol = 2) +
    geom_smooth(se = FALSE, method = "lm", linewidth = 0.5, alpha = 0.9) +
    scale_colour_manual(values = c("#2E7F18", "#45731E", "#675E24",
                                   "black",
                                   "#8D472B", "#B13433", "#C82538")) +
    scale_linetype_manual(values = c(4,3,2,1,2,3,4))
}

# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Reference values for maximal oxygen uptake"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), 
                            "/ last update 2023-11-23")
                ),
                
                tags$br(),

                sidebarLayout(
                  sidebarPanel(width = 2,
                               selectInput("sex",
                                           label = "Sex", 
                                           selected = "Female", 
                                           choices = c("Female", "Male")),
                               numericInput("age",
                                            label = "Age (years)",
                                            value = 32),
                               numericInput("vo2max",
                                            label = "VO2max (mL/min/kg)",
                                            value = 33)
                               ),
                  mainPanel(width = 10,
                            plotOutput("VO2plot")
                            )
                ),
                
                fluidRow(style="padding:1% 2% 1% 2%",
                  tags$hr(),
                  h4("The plot above was generated from data extracted from the following publications:"),
                  plotOutput("studies_plot", height = 500),
                  tags$hr(),
                  h4("The following studies were screened for data:"),
                  tabsetPanel(type = "tabs",
                              tabPanel("Studies included", style="padding:1% 2% 1% 2%;text-align:left",
                                       DT::dataTableOutput("studies_included")
                              ),
                              tabPanel("Studies excluded", style="padding:1% 2% 1% 2%;text-align:left",
                                       dataTableOutput("studies_excluded")
                              )
                  )

                )
)


# Define server logic ----
server <- function(input, output, session) {
  

  output$VO2plot <- renderPlot({
    
    # input values
    input_age = 30
    input_vo2max = 40
    input_sex = "Male"
    input_age = input$age
    input_vo2max = input$vo2max
    input_sex = input$sex

    # error messages
    validate(need(input_age, "Please provide a numeric value age"))
    validate(need(input_vo2max, "Please provide a numeric value for VO2max"))
    
    # subset selected sex
    VO2_data_M <- subset(VO2_data, Sex == "Male")
    VO2_data_F <- subset(VO2_data, Sex == "Female")
    
    # plot
    nomogram_legend <- get_legend(
      fct_nomogram(VO2_data_M) + 
        guides(color = guide_legend(ncol = 1)) 
    )

    # male plot
    nomogram_M <- fct_nomogram(VO2_data_M) + 
      labs(subtitle = "Male") +
      theme(legend.position = "none")
    
    nomogram_M <- if(input_sex == "Male")
      nomogram_M + 
      annotate("label", 
               x = input_age, 
               y = input_vo2max -2, 
               label = "Your VO2max", 
               color = "red") +
      annotate("pointrange", 
               x = input_age, 
               y = input_vo2max, 
               ymin = input_vo2max, ymax = input_vo2max,
               color = "red", 
               size = 0.5) else nomogram_M
    
    # female plot
    nomogram_F <- fct_nomogram(VO2_data_F) + 
      labs(subtitle = "Female") +
      theme(legend.position = "none")
    
    nomogram_F <- if(input_sex == "Female")
      nomogram_F + 
      annotate("label", 
               x = input_age, 
               y = input_vo2max -2, 
               label = "Your VO2max", 
               color = "red") +
      annotate("pointrange", 
               x = input_age, 
               y = input_vo2max, 
               ymin = input_vo2max, ymax = input_vo2max,
               color = "red", 
               size = 0.5) else nomogram_F
    
    
    # merge plots
    cowplot::plot_grid(
      nomogram_F,
      nomogram_M,
      nomogram_legend,
      ncol = 3,
      rel_widths = c(2,2,1)
    )
  })
  
  output$studies_included <- DT::renderDataTable(escape = FALSE, 
                                                 rownames = FALSE, 
                                                 options=list(paging = FALSE,
                                                              dom = 't'), 
                                                 { studies_included })
  
  output$studies_excluded <- DT::renderDataTable(escape = FALSE, 
                                                 rownames = FALSE, 
                                                 options=list(paging = FALSE,
                                                              dom = 't'), 
                                                 { studies_excluded })
  
  output$studies_plot <- renderPlot({
    # plot mean only for all studies
    ggplot(VO2_data[VO2_data$Percentile == "Median",], 
           aes(x = Age, y = VO2max, 
               color = Reference, shape = Reference)) +
      labs(x = "Age (Years)",
           y = "Median VO2max (mL/min/kg)") +
      theme_bw(16) +
      theme(legend.key.size = unit(20, "pt")) +
      facet_wrap(.~Sex) +
      geom_point(size = 1.5) +
      geom_line(linewidth = 0.25) +
      scale_x_continuous(breaks = seq(20, 90, 10)) +
      scale_shape_manual(values = rep(c(15,16,17,18), 20))
  })

}

# Run the app ----
shinyApp(ui = ui, server = server)