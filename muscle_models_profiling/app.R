#-----------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle cells
#
#----------------------------------------------------------------------
# Load data and libraries
library(shinycssloaders)
library(ggplot2)
theme <- theme(plot.title  = element_text(face="bold", color="black", size=13, angle=0),
               axis.text.x = element_text(color="black", size=12, angle=90, hjust=1, vjust=0.5),
               axis.text.y = element_text(color="black", size=11, angle=0, vjust=0.3),
               axis.title  = element_text(face="bold", color="black", size=13, angle=0),
               legend.text = element_text(color="black", size=13, angle=0),
               legend.title = element_blank(),
               legend.position="none",
               legend.key.size = unit(30, "pt"),
               strip.text = element_text(face="bold", color="black", size=14, angle=0))
library(stringr)
library(plyr)
library(DT)
library(dplyr)
library(DescTools)
library(ggpubr)
library(rstatix)
library(ggprism)

# Load data ----
data_all <- readRDS('data/Muscle_Models_Profiling_data.Rds')
data_all <- cbind(data_all[grepl('HumanCell', colnames(data_all))],
                  data_all[grepl('MouseC2C12', colnames(data_all))],
                  data_all[grepl('RatL6', colnames(data_all))],
                  data_all[grepl('HumanTissue', colnames(data_all))],
                  data_all[grepl('MouseTissue', colnames(data_all))],
                  data_all[grepl('RatTissue', colnames(data_all))])

samples_all <- gsub("_G.*", "", colnames(data_all))
samples_all <- gsub(".*_", "", samples_all)

#--------------------------------------------------------------------------------------------------------
# List of genes
#--------------------------------------------------------------------------------------------------------
genelist <- unique(rownames(data_all))


# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Gene expression in mouse, rat and human skeletal muscle tissue and cells"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), "/ last update 2022-03-08")
                ),
                
                tags$hr(),
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                         actionButton("updatePlot", "Refresh plot", icon("refresh"))
                ),
                
                fluidRow(style="color:black;background-color:white;padding:2% 8% 1% 8%;",
                         plotOutput("genePlot", height="700px") %>% withSpinner(color="#5b768e")
                  )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected=c("CCL2", "IL6"), 
                       options=NULL)
  
  plotData <- eventReactive(input$updatePlot, {
    validate(need(input$inputGeneSymbol, " "))
    genename <- toupper(input$inputGeneSymbol)
    
    #plot
    data <- data.frame()
    for( i in genename) { 
      asd <- data.frame(x=samples_all, 
                        y=as.numeric(data_all[i,]),
                        Gene=i, 
                        stringsAsFactors=FALSE) #empty dataframe to collect data
      data  <- rbind.data.frame(data, asd)              #bind the new gene data at the bottom of the previous one
    }
    
    gg <- ggplot(data) +  geom_boxplot(aes(x=x, y=y, fill=x)) + theme_bw() +
      scale_x_discrete(breaks=sort(unique(samples_all)),
                       labels=c("Human Primary Myotube", "Human Muscle Tissue",
                                "Mouse C2C12 Myotube", "Mouse Muscle Tissue", 
                                "Rat L6 Myotube", "Rat Muscle Tissue")) +
      labs(x="",
           y="mRNA expression (log2)",
           title=element_blank()) +
      scale_y_continuous(breaks = round(seq(-4, 8, by=2),1)) +
      theme +
      scale_fill_manual(values=c("#56B4E9", "#D3C839", "#CC79A7", "#0072B2", "#E69F00", "#D55E00")) +
      scale_color_manual(values=c("#56B4E9", "#D3C839", "#CC79A7", "#0072B2", "#E69F00", "#D55E00")) +
      geom_hline(aes(yintercept=0), linetype="dashed", show.legend=F, color="gray60") +
      facet_wrap(~Gene)
    return(gg)
    gg
  })
  
  
  output$genePlot <- renderPlot({
    plotData()
  })
  
}


# Run the app ----
shinyApp(ui = ui, server = server)