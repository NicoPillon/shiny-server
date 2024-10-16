#-----------------------------------------------------------------------
#
# Mouse Obesity
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


#load data
HFD_data <- readRDS('data/data_raw.Rds')
genelist <- rownames(HFD_data)
samples <- readRDS("data/metadata.Rds")
datasets <- readRDS("data/datasets.Rds")

#organize diet duration
samples$diet.duration.cat <- "< 6 weeks"
samples$diet.duration.cat[samples$diet.duration >= 6] <- "6-12 weeks"
samples$diet.duration.cat[samples$diet.duration >= 13] <- "13-18 weeks"
samples$diet.duration.cat[samples$diet.duration > 18] <- "> 18 weeks"

#organize age
samples$age.at.start.cat <- "3-4 weeks"
samples$age.at.start.cat[samples$age.at.start > 4] <- "5-6 weeks"
samples$age.at.start.cat[samples$age.at.start > 6] <- "8-9 weeks"

#organize diet
samples$diet.composition.cat <- "32-45 %"
samples$diet.composition.cat[samples$diet.composition > 45] <- "50-60 %"

# function to format p values
p_value_formatter <- function(p) {
  sapply(p, function(x) {
    if (x < 0.001) {
      return("p < 0.001")
    } else {
      return(sprintf("p == %.3f", x))
    }
  })
}


# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                
                # Custom CSS to change checkbox tick color
                tags$style(HTML("
                  input[type='checkbox'] {
                    accent-color: #c93f1e; /* Change the checkbox tick color here */
                  }
                ")),
                
                # title ribbon
                fluidRow(style="color:white;background-color:#5B768E;padding:0% 1% 1% 1%;text-align:center",
                         h3("Gene expression in tissues from mice fed a high fat diet"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), 
                            "/ last update 2024-10-16")
                ),
                
                # main page
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         sidebarLayout(
                           sidebarPanel(width = 3,
                                        selectizeInput("inputGeneSymbol", "Gene Symbol:", choices=NULL, multiple=F, width=1000),
                                        checkboxGroupInput("diet_duration", 
                                                           label = "Diet duration", 
                                                           selected = c("< 6 weeks", "6-12 weeks", "13-18 weeks", "> 18 weeks"), 
                                                           choices = c("< 6 weeks", "6-12 weeks", "13-18 weeks", "> 18 weeks")),
                                        checkboxGroupInput("age_start", 
                                                           label = "Age at start", 
                                                           selected = c("3-4 weeks", "5-6 weeks", "8-9 weeks"), 
                                                           choices = c("3-4 weeks", "5-6 weeks", "8-9 weeks")),
                                        checkboxGroupInput("diet_composition", 
                                                           label = "Diet fat composition", 
                                                           selected = c("32-45 %", "50-60 %"), 
                                                           choices = c("32-45 %", "50-60 %"))
                           ),
                           mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                                     plotOutput("genePlot", height="500px") %>% withSpinner(color="#5b768e")
                                     )
                         )
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
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected=c("Itgax"), 
                       options=NULL)

  output$genePlot <- renderPlot({
    validate(need(input$inputGeneSymbol, " "))
    genename <- c("Itgax")
    genename <- input$inputGeneSymbol
    
    plotdata <- data.frame()
    for (i in genename){
      asd <- HFD_data[i,]
      asd <- data.frame(samples, as.numeric(asd), Gene=i)
      colnames(asd) <- c(colnames(samples), "data", "Gene")
      plotdata <- rbind(plotdata, asd)
    }
    
    #filter according to selected categories
    plotdata <- plotdata[plotdata$diet.duration.cat %in% input$diet_duration &
                           plotdata$age.at.start.cat %in% input$age_start  &
                           plotdata$diet.composition.cat %in% input$diet_composition,]
    
    # label with n size
    plotdata$sex <- gsub("^male", paste0("Male, n = ", nrow(plotdata[plotdata$sex == "male",])), plotdata$sex)
    plotdata$sex <- gsub("^female", paste0("Female, n = ", nrow(plotdata[plotdata$sex == "female",])), plotdata$sex)
    
    
    
    #plot
    ggplot(plotdata, aes(x=tissue, y=data, fill=diet)) +  
      geom_boxplot(outlier.size = 0.1, position = position_dodge(0.85))  + 
      geom_dotplot(binaxis = "y", 
                   stackdir = "center", 
                   dotsize = 2, 
                   position = position_dodge(0.85), 
                   binwidth = (max(plotdata$data, na.rm = T) - min(plotdata$data, na.rm = T))/100, 
                   alpha = 0.5) +
      theme_bw(16) + 
      facet_wrap(~sex, scales="free_y") +
      labs(x=element_blank(),
           y="mRNA expression, log2") +
      scale_fill_manual(values=c("gray90", "#D55E00"))  +
      stat_compare_means(aes(group = diet, 
                             label = ifelse(p < 0.001, "p < 0.001", sprintf("p = %5.3f", as.numeric(..p..)))), 
                         size = 4, vjust=-2) +
      scale_y_continuous(expand = c(0,1)) 
    
  })
  
  
  ##################################################################################################################
  #Dataset tables
  output$datasets <- renderDataTable(options=list(signif = 3),{
    DT::datatable(datasets, 
                  options = list(lengthMenu = c(10, 50, 100), 
                                 pageLength = 10),
                  rownames= FALSE)
  })
  
}


# Run the app ----
shinyApp(ui = ui, server = server)
