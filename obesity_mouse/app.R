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
samples <- readRDS("data/metadata.Rds")
datasets <- readRDS("data/datasets.Rds")

#organize diet duration
samples$diet.duration.cat <- "< 6 weeks"
samples$diet.duration.cat[samples$diet.duration > 6] <- "6-12 weeks"
samples$diet.duration.cat[samples$diet.duration > 12] <- "12-18 weeks"
samples$diet.duration.cat[samples$diet.duration > 18] <- "> 18 weeks"

#organize age
samples$age.at.start.cat <- "3-4 weeks"
samples$age.at.start.cat[samples$age.at.start > 4] <- "5-6 weeks"
samples$age.at.start.cat[samples$age.at.start > 6] <- "8-9 weeks"

#organize diet
samples$diet.composition.cat <- "< 50 %"
samples$diet.composition.cat[samples$diet.composition >= 50] <- "50-55 %"
samples$diet.composition.cat[samples$diet.composition >= 55] <- "> 55 %"


# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Gene expression in tissues from mice fed a high fat diet"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), 
                            "/ last update 2022-06-02")
                ),
                
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;;text-align:center",
                         em(h5("If you know other datasets that could be incorporated in this tool, please",
                               a("let me know!", href="mailto:nicolas.pillon@ki.se", 
                                 target="_blank", style="color:#5B768E"))),
                         tags$hr()
                ),
                
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=1000),
                         column(2, checkboxGroupInput("diet_duration", 
                                            label = "Diet duration", 
                                            selected = c("< 6 weeks",
                                                         "6-12 weeks", 
                                                         "12-18 weeks",
                                                         "> 18 weeks"), 
                                            choices = c("< 6 weeks",
                                                        "6-12 weeks", 
                                                        "12-18 weeks",
                                                        "> 18 weeks"))),
                         column(2, checkboxGroupInput("age_start", 
                                            label = "Age at start", 
                                            selected = c("3-4 weeks",
                                                         "5-6 weeks",
                                                         "8-9 weeks"), 
                                            choices = c("3-4 weeks",
                                                        "5-6 weeks",
                                                        "8-9 weeks"))),
                         column(2, checkboxGroupInput("diet_composition", 
                                                      label = "Diet composition", 
                                                      selected = c("< 50 %",
                                                                   "50-55 %",
                                                                   "> 55 %"), 
                                                      choices = c("< 50 %",
                                                                  "50-55 %",
                                                                  "> 55 %"))),
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
    data <- data[data$diet.duration.cat %in% input$diet_duration &
                   data$age.at.start.cat %in% input$age_start  &
                   data$diet.composition.cat %in% input$diet_composition,]
    
    #plot
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
