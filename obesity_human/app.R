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
theme <- theme(plot.title  = element_text(face="bold", color="black", size=13, angle=0),
               axis.text.x = element_text(color="black", size=12, angle=90, hjust=1, vjust=0.5),
               axis.text.y = element_text(color="black", size=11, angle=0, vjust=0.3),
               axis.title  = element_text(face="bold", color="black", size=13, angle=0),
               legend.text = element_text(color="black", size=13, angle=0),
               legend.title = element_blank(),
               legend.position="right",
               legend.key.size = unit(30, "pt"),
               strip.text = element_text(face="bold", color="black", size=14, angle=0))

#load data
human_data <- readRDS('data/data_raw.Rds')
genelist <- rownames(human_data)
samples <- readRDS("data/metadata.Rds")
datasets <- readRDS("data/datasets.Rds")


# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Gene expression in tissues from human with obesity"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), 
                            "/ last update 2022-06-03")
                ),
                
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;;text-align:center",
                         em(h5("WORK IN PROGRESS: I will keep adding more datasets in the coming weeks."),
                            h5("If you know other datasets that could be incorporated in this tool, please",
                               a("let me know!", href="mailto:nicolas.pillon@ki.se", 
                                 target="_blank", style="color:#5B768E"))),
                         tags$hr()
                ),
                
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=1000),
                         column(2, checkboxGroupInput("diagnosis_obesity", 
                                                      label = "Obesity diagnosis", 
                                                      selected = c("Lean",
                                                                   "Overweight", 
                                                                   "Obesity class I",
                                                                   "Obesity class II",
                                                                   "Obesity class III"), 
                                                      choices = c("Lean",
                                                                  "Overweight", 
                                                                  "Obesity class I",
                                                                  "Obesity class II",
                                                                  "Obesity class III"))),
                         column(2, checkboxGroupInput("diagnosis_diabetes", 
                                            label = "Diabetes diagnosis", 
                                            selected = c("Healthy",
                                                         "Prediabetes", 
                                                         "T2D"), 
                                            choices = c("Healthy",
                                                        "Prediabetes",
                                                        "T2D"))),

                         actionButton("updatePlot", "Refresh plot", icon("refresh"))
                ),
                fluidRow(style="color:black;background-color:white;padding:2% 8% 1% 8%;",
                         plotOutput("genePlot", height="400px") %>% withSpinner(color="#5b768e")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         tags$hr(),
                         h3("Datasets Included in the Analysis"),
                         dataTableOutput("datasets")
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected=c("ITGAM", "IGFBP1"), 
                       options=NULL)
  
  plotInput <- eventReactive(input$updatePlot, {
    validate(need(input$inputGeneSymbol, " "))
    genename <- c("ITGAX", "LEP")
    genename <- input$inputGeneSymbol

    data <- data.frame()
    for (i in genename){
      asd <- human_data[i,]
      asd <- data.frame(samples, as.numeric(asd), Gene=i)
      colnames(asd) <- c(colnames(samples), "data", "Gene")
      data <- rbind(data, asd)
    }
    
    #filter according to selected categories
    data <- data[data$ObesityClass %in% input$diagnosis_obesity &
                   data$diagnosis.diabetes %in% input$diagnosis_diabetes,]
    
    #order the categories
    data$diagnosis.weight <- factor(data$diagnosis.weight, 
                                    levels=c("Lean", "Overweight", "Obesity"))
    
    #plot
    gg <- ggplot(data, aes(x=tissue, y=data, fill=diagnosis.weight)) +  
      geom_boxplot()  + 
      theme_bw() + theme +
      facet_wrap(~Gene, scales="free_y") +
      labs(x=element_blank(),
           y="mRNA expression, log2") +
      scale_fill_manual(values=c("#56B4E9", "#E69F00", "#D55E00"))  +
      stat_compare_means(aes(group = diagnosis.weight, label = ifelse(
        p < 0.001, "p < 0.001", sprintf("p = %5.3f", as.numeric(..p..)))), 
        size = 3, vjust=-1) +
      scale_y_continuous(expand = c(0,1)) 
    gg
    
  })

  output$genePlot <- renderPlot({
    plotInput()
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
