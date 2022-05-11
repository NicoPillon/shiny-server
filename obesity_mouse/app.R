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
HFD_data <- readRDS('data/data_raw.Rds')
genelist <- rownames(HFD_data)
samples <- data.frame(str_split_fixed(colnames(HFD_data), "\\.|_",4))[,1:3]
colnames(samples) <- c('tissue', 'diet', 'GEO')
datasets <- readRDS("data/datasets.Rds")

# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Gene expression in tissues from mice fed a high fat diet"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), 
                            "/ last update 2021-12-02")
                ),
                tags$br(),
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                         #checkboxGroupInput("diet_duration", 
                        #                    label="Diet duration (weeks)", 
                        #                    selected=c("4", "8", "10", "11", "12", "14", "15"), 
                        #                    c("4", "8", "10", "11", "12", "14", "15")),
                         actionButton("updatePlot", "Refresh plot", icon("refresh"))
                ),
                fluidRow(style="color:black;background-color:white;padding:2% 8% 1% 8%;",
                         plotOutput("genePlot", height="400px") %>% withSpinner(color="#5b768e")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         tags$hr(),
                         h3("Datasets"),
                         dataTableOutput("datasets")
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected=c("Itgax", "Lep"), 
                       options=NULL)
  
  plotInput <- eventReactive(input$updatePlot, {
    validate(need(input$inputGeneSymbol, " "))
    genename <- c("Itgax", "Lep")
    genename <- input$inputGeneSymbol

    data <- data.frame()
    for (i in genename){
      asd <- HFD_data[i,]
      asd <- data.frame(samples, as.numeric(asd), Gene=i)
      colnames(asd) <- c(colnames(samples), "data", "Gene")
      data <- rbind(data, asd)
    }
    
    #filter according to selected categories
    #testdata <- testdata[testdata$Treatment %in% input$diet_duration,]
    
    gg <- ggplot(data, aes(x=tissue, y=data, fill=diet)) +  
      geom_boxplot()  + 
      theme_bw() + theme +
      facet_wrap(~Gene, scales="free_y") +
      labs(x=element_blank(),
           y="mRNA expression, log2") +
      scale_fill_manual(values=c("#56B4E9", "#D55E00"))  +
      stat_compare_means(aes(group = diet, label = ifelse(
        p < 0.001, "p < 0.001", sprintf("p = %5.3f", as.numeric(..p..)))), 
        size = 3, vjust=-2) +
      scale_y_continuous(expand = c(0,3)) 
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
