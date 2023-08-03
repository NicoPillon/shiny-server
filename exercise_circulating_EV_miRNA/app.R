#-----------------------------------------------------------------------
#
# microRNA profile of circulating EVs in response to exercise
#
#----------------------------------------------------------------------
# Load data and libraries
library(shinycssloaders)
library(stringr)
library(metafor)
library(metaviz)


#--------------------------------------------------------------------------------------------------------
# Palmitate
#--------------------------------------------------------------------------------------------------------
EV_data <- readRDS("data/data.Rds")

#get stats from limma
EV_stats <- readRDS("data/stats.Rds")


#--------------------------------------------------------------------------------------------------------
# List of genes
#--------------------------------------------------------------------------------------------------------
genelist <- unique(EV_stats$miR)



# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Meta-analysis of microRNA in circulating extracellular vesicles after acute exercise"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), "/ last update 2023-04-20")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 0% 8%;",
                         "Sequencing data from",
                         a("GSE121874", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121874", 
                           target="_blank", style="color:#5B768E"), 
                         a("GSE136997", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136997", 
                           target="_blank", style="color:#5B768E"), 
                         a("GSE144627", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144627", 
                           target="_blank", style="color:#5B768E"), 
                         a("GSE209880", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE209880", 
                           target="_blank", style="color:#5B768E"), 
                         a("GSE218103", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218103", 
                           target="_blank", style="color:#5B768E")
                ),
                tags$hr(),
                fluidRow(style="color:black;background-color:white;padding:0% 8% 0% 8%;",
                         selectizeInput("inputGeneSymbol", "", choices=NULL, multiple=F)
                ),
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         plotOutput("EvPlot", height="400px") %>% withSpinner(color="#5b768e")
                ),
                tags$hr(),
                fluidRow(style="color:black;background-color:white;padding:0% 0% 1% 8%;display: flex;",
                         column(4, style = "padding:0% 0% 0% 0%",
                                plotOutput("funnelPlot") %>% withSpinner(color="#5b768e")
                         ),
                         column(6, style = "padding:15% 0% 0% 0%",
                                tableOutput("metaStats") %>% withSpinner(color="#5b768e")
                         )
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected=c("hsa-miR-122-5p"), 
                       options=NULL)
  

  metaAnalysis <- function(){
    validate(need(input$inputGeneSymbol, " "))
    miRointerest <- c("hsa-miR-122-5p")
    miRointerest <- input$inputGeneSymbol
    
    plotdata <- data.frame(
      GSE121874 = as.numeric(EV_data[["GSE121874"]][miRointerest, ]),
      GSE136997 = as.numeric(EV_data[["GSE136997"]][miRointerest, ]),
      GSE144627.H00 = as.numeric(EV_data[["GSE144627.H00"]][miRointerest, ]),
      GSE144627.H03 = as.numeric(EV_data[["GSE144627.H03"]][miRointerest, ]),
      GSE218103 = as.numeric(EV_data[["GSE218103"]][miRointerest, ]),
      GSE218103.CebPalsy = as.numeric(EV_data[["GSE218103.CebPalsy"]][miRointerest, ]),
      GSE209880.H00 = as.numeric(EV_data[["GSE209880.H00"]][miRointerest, ]),
      GSE209880.H03 = as.numeric(EV_data[["GSE209880.H03"]][miRointerest, ]),
      GSE209880.H24 = as.numeric(EV_data[["GSE209880.H24"]][miRointerest, ])
    )
    rownames(plotdata) <- gsub(".*_", "", colnames(EV_data[["GSE121874"]]))
    
    selectedata <- data.frame(t(plotdata))
    selectedata <- na.omit(selectedata)
    
    meta   <- rma(data = selectedata, 
                  yi = logFC,
                  sei = SE,
                  method = "REML", measure = "MD")
    meta$slab <- rownames(selectedata)
    return(meta)
    }
    
    
    output$EvPlot <- renderPlot({
      meta <- metaAnalysis()
    #Make tables and forest plot figures from rma
    study_table <- data.frame(
      Reference = c("Shah et al, 2017", 
                    "Just et al, 2020", 
                    "Nair et al, 2020", 
                    "Nair et al, 2020", 
                    "Vechetti et al, 2022",
                    "Vechetti et al, 2022",
                    "Lavin et al, 2023",
                    "Lavin et al, 2023",
                    "Lavin et al, 2023"
      ),
      GEO = meta$slab,
      logFC = format(round(meta$data$logFC, digits=2)),
      P.Val = format(meta$data$P.Value,   scientific=T, digits=2),
      n = c("5", "6", "10", "10", "3", "3", "40", "40", "40"))
    
    
    
    summary_table <- data.frame(
      Reference = input$inputGeneSymbol,
      GEO = "Meta-analysis",
      logFC = format(round(meta$beta, digits=2)),
      adj.P.Val = format(meta$pval,   scientific=F, digits=2),
      n = 14)
    
    #forest plot with tables
    return(viz_forest(meta,
                      text_size = 5,
                      study_table = study_table,
                      summary_table = summary_table,
                      col = "#D55E00",
                      xlab="logFC"))
  })
  
    
    output$funnelPlot <- renderPlot({
      meta <- metaAnalysis()
      funnel(meta)
    })

    output$metaStats <- renderTable(colnames = FALSE,{
      meta <- metaAnalysis()
      metastats <- data.frame( c("tau^2 (estimated amount of total heterogeneity):", signif(meta$tau2, 2)),
                               c("I^2 (total heterogeneity / total variability):", signif(meta$I2, 2)), 
                               c("H^2 (total variability / sampling variability):", signif(meta$H2, 2)),
                               c("Test for Heterogeneity:", paste0("Q = ", signif(meta$QE, 4), ", p = ", signif(meta$QEp, 2)))
      )
      t(metastats)
    })
    
    
}


# Run the app ----
shinyApp(ui = ui, server = server)