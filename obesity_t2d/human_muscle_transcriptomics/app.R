#-----------------------------------------------------------------------
#
# Human Obesity
#
#----------------------------------------------------------------------
# Load libraries
library(shinycssloaders)
library(stringr)
library(DT)
library(plyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(feather)
library(ggforce)

# function to format p values
p_value_formatter <- function(p) {
  sapply(p, function(x) {
    if (x < 0.001) {
      return("italic(p) < 0.001")
    } else {
      return(sprintf("italic(p) == %.3f", x))
    }
  })
}

#load metadata
human_genelist <- readRDS("data/human_genelist.Rds")
human_metadata <- readRDS("data/human_metadata.Rds")
human_references <- readRDS("data/human_references.Rds")
  
# matrix
human_datamatrix_1 <- read_feather('data/human_datamatrix_1.feather')
human_datamatrix_2 <- read_feather('data/human_datamatrix_2.feather')
human_datamatrix <- data.frame(rbind(human_datamatrix_1,
                                     human_datamatrix_2))
rownames(human_datamatrix) <- human_genelist

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
                fluidRow(style="color:white;background-color:#5B768E;padding:0% 1% 1% 1%;text-align:center",
                         column(1, 
                                style = "height:8vh; display:flex; justify-content:center; align-items:center;",
                                tags$a(href = "https://shiny.nicopillon.com", 
                                       icon("home", class = "fa-2x"), 
                                       style = "color:white; text-decoration:none;")  # Ensuring icon is white and no underline
                                ),
                         column(10,
                                h3("Gene expression in skeletal muscle from humans with obesity and type 2 diabetes"),
                                h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                           target="_blank", style="color:#D9DADB"), 
                                   "/ last update 2024-10-17")
                                ),
                         column(1,
                         )
                ),

                # main page
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         sidebarLayout(
                           sidebarPanel(width = 3,
                                        selectizeInput("inputHumanGeneSymbol", "Gene Symbol:", choices=NULL, multiple=F, width=1000),
                                        sliderInput("age", tags$b("Age (years)"),
                                                    min = 18, max = 90, value = c(18,90), step = 1, sep = ""),
                                        checkboxGroupInput("diagnosis_diabetes", 
                                                           label = "Diabetes diagnosis", 
                                                           selected = c("Healthy", "Prediabetes", "T2D"),
                                                           choices = c("Healthy", "Prediabetes", "T2D")),
                                        em(h5("Prediabetes is defined as either impaired glucose tolerance measured during an OGTT or increased insulin 
                                              resistance measured with euglycemic hyperinsulinemic clamps."))
                                        
                           ),
                           mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                                     plotOutput("genePlotObesity", height="600px") %>% withSpinner(color="#5B768E")
                                     )
                         ),
                         
                ),
                
                tags$hr(),
                
                # Table with datasets
                fluidRow(style="color:black;background-color:white;padding:0% 2% 1% 2%;",
                         h3("Datasets Included in the Analysis"),
                         dataTableOutput("datasets")
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
  
  updateSelectizeInput(session, 'inputHumanGeneSymbol', 
                       choices=human_genelist, 
                       server=TRUE, 
                       selected=c("LEP"), 
                       options=NULL)
  
  output$genePlotObesity <- renderPlot({
    validate(need(input$inputHumanGeneSymbol, " "))
    genename <- c("AKR1C3")
    genename <- input$inputHumanGeneSymbol

    # Obesity
    plotdata <- data.frame(human_metadata, 
                           genedata = as.numeric(human_datamatrix[genename,]))
    plotdata$bmi_category <- factor(plotdata$bmi_category, 
                                    levels=c("Lean", "Overweight", "Obesity"))

    #filter according to selected categories
    plotdata <- dplyr::filter(plotdata,
                              diagnosis %in% input$diagnosis_diabetes ,
                              age >= input$age[1] & age <= input$age[2])
    
    # label with n size
    plotdata$sex <- gsub("^male", paste0("Male, n = ", nrow(plotdata[plotdata$sex == "male",])), plotdata$sex)
    plotdata$sex <- gsub("^female", paste0("Female, n = ", nrow(plotdata[plotdata$sex == "female",])), plotdata$sex)
    
    cowplot::plot_grid(
      ggplot(plotdata, aes(x=bmi, y=genedata)) +  
        geom_point(aes(color = diagnosis, shape = diagnosis), size = 3, alpha = 0.25)  + 
        geom_smooth(method = "lm", color = "black", se = FALSE) +
        theme_bw(16) +  
        theme(legend.position = "none") +
        facet_wrap(.~sex, ncol = 1) +
        labs(x="BMI, kg/m2",
             y="mRNA expression, log2") +
        scale_shape_manual(values=rep(c(15,16,17), 20)) +
        scale_color_manual(values = c("darkgreen", "orange", "darkred")) +
        scale_y_continuous(expand = c(0,1.5)) +
        stat_cor(size = 4, 
                 vjust = -1, 
                 label.x = 20),
      
      ggplot(plotdata, aes(x=bmi_category, y=genedata)) +  
        #geom_jitter(aes(color = diagnosis, shape = diagnosis), size = 2, alpha = 0.25, width = 0.1) +
        geom_boxplot(outlier.size = 0.1, fill = "gray80", alpha = 0.5)  + 
        geom_sina(aes(color = diagnosis, shape = diagnosis), size = 1.5, position = position_dodge(0), alpha = 0.25) +
        theme_bw(16) + 
        theme(legend.position = "right") +
        facet_wrap(.~sex, ncol = 1) +
        labs(x="BMI category",
             y="mRNA expression, log2") +
        scale_shape_manual(values=rep(c(15,16,17), 20)) +
        scale_color_manual(values = c("darkgreen", "orange", "darkred")) +
        scale_y_continuous(expand = c(0,1.5)) +
        stat_compare_means(aes(label = after_stat(p_value_formatter(..p..))), 
                           ref.group = "Lean",
                           parse = TRUE,
                           size = 4, 
                           vjust = -1),
      
      ncol = 2
    )
    
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
