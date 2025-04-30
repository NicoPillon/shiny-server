#-----------------------------------------------------------------------
#
# Single cell RNAseq before/after exercise
#
#----------------------------------------------------------------------
# Load data and libraries
library(shinycssloaders)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggprism)

#--------------------------------------------------------------------------------------------------------
# Data from Lovric et al
#--------------------------------------------------------------------------------------------------------
scRNAseq_data <- readRDS("data/data.Rds")

#--------------------------------------------------------------------------------------------------------
# List of genes
#--------------------------------------------------------------------------------------------------------
genelist <- unique(scRNAseq_data$Gene_name)


# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Single cell RNA sequencing in skeletal muscle before and after an acute bout of exercise"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), "/ last update 2022-10-27")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 0% 8%;",
                         "From Lovric et al, ",
                         a("Single-cell sequencing deconvolutes cellular responses to exercise in human skeletal muscle", 
                           href="https://www.nature.com/articles/s42003-022-04088-z#Sec26", 
                           target="_blank", style="color:#5B768E"), 


                ),
                tags$hr(),
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                         actionButton("updatePlot", "Refresh plot", icon("refresh"))
                ),
                fluidRow(style="color:black;background-color:white;padding:2% 8% 1% 8%;",
                         plotOutput("PalPlot", height="600px") %>% withSpinner(color="#5b768e"),
                         h5("Statistics presented on the plot are false discovery rate (FDR) and scores: value incorporating the 
                         foldchange difference, statistical (e.g. p-val adjusted of a wilcoxon rank sum) difference, and extent 
                         of expression between the cell type post-exercise against pre-exercise.")
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected=c("PPARGC1A", "HES1", "CX3CL1", "IL1B"), 
                       options=NULL)
  
  # AMPK plot
  plotData <- eventReactive(input$updatePlot, {
    validate(need(input$inputGeneSymbol, " "))
    genename <- c("PPARGC1A", "HES1", "CX3CL1", "IL1B")
    genename <- input$inputGeneSymbol
    
    #plot
    scRNAseq_selected <- scRNAseq_data[scRNAseq_data$Gene_name %in% genename,]
    scRNAseq_selected$stats <- paste0("FDR = ", signif(scRNAseq_selected$pvals_adj, 2))
    scRNAseq_selected$position <- scRNAseq_selected$scores
    scRNAseq_selected$position[scRNAseq_selected$position < 0] <- 0
    ggbarplot(scRNAseq_selected, x="cell.type", y="scores", fill="cell.type",
              position = position_dodge(0.7),
              facet.by = "Gene_name",
              size = 0.2) +
      facet_wrap(.~Gene_name, scale="free") +
      theme_bw(14) + theme(legend.position = "none") +
      labs(x = element_blank(),
           y = "Response to exercise, log(score)") +
      geom_hline(yintercept = 0) +
      add_pvalue(scRNAseq_selected, bracket.size = NA, 
                 xmin = "cell.type",
                 xmax = "cell.type",
                 label = "stats",
                 y.position = "position", 
                 label.size = 3.5) + 
      scale_y_continuous(expand = c(0.1, 0.1))
  })
  
  
  output$PalPlot <- renderPlot({
    plotData()
  })
  
  #Statistic tables
  output$stats <- renderDataTable(options=list(signif = 3),{
    datatable(PAL_stats, 
              options = list(lengthMenu = c(10, 50, 100), 
                             pageLength = 10),
              rownames= TRUE) %>% 
      formatSignif(columns = c(1:8), digits = 3)
  })
  
}


# Run the app ----
shinyApp(ui = ui, server = server)