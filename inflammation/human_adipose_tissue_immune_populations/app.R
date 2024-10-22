#-----------------------------------------------------------------------
#
# FACS of adipose tissue resident cells
#
#----------------------------------------------------------------------
# Load libraries
library(shinycssloaders)
library(stringr)
library(DT)
library(plyr)
library(ggplot2)

#load data
WATFACS_data <- readRDS('data/data_raw.Rds')

samples <- readRDS('data/samples.Rds')
samples$group <- "Adipose cells"
samples$group[samples$CellType %in% c("Leukocyte")] <- "Leukocytes"
samples$group[samples$CellType %in% c("Total T-cells", "CD4+ T-cells", "CD8+ T-cells")] <- "Lymphocytes"
samples$group[samples$CellType %in% c("Monocyte/Macrophage", "M1 Macrophage", "M2 Macrophage", "CD14+ Myeloid")] <- "Macrophages"
samples$group[samples$CellType %in% c("SVF")] <- "SVF"

genelist <- rownames(WATFACS_data)


# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5B768E;padding:0% 1% 1% 1%;text-align:center; display:flex; justify-content:center; align-items:center;",
                         h3("FACS of cells from human adipose tissue"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), 
                            "/ last update 2022-03-08")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 0% 8%;",
                         "Cells from digested human adipose tissue: stroma-vascular fraction and mature adipocytes. From",
                         a("GSE80654", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80654", 
                           target="_blank", style="color:#5B768E"), 
                         "(HTA-2.0) and", 
                         a("GSE100795", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100795", 
                           target="_blank", style="color:#5B768E"), 
                         "(RNAseq)"
                ),
                tags$hr(),
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                         actionButton("updatePlot", "Refresh plot", icon("refresh"))
                ),
                
                fluidRow(style="color:black;background-color:white;padding:2% 8% 1% 8%;",
                         plotOutput("genePlot", height="600px") %>% withSpinner(color="#5b768e", proxy.height=300)
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected=c("ADIPOQ", "MSR1", "ITGAX"), 
                       options=NULL)
  
  plotData <- eventReactive(input$updatePlot, {
    validate(need(input$inputGeneSymbol, " "))
    genename <- c("ADIPOQ", "MSR1", "VWF")
    genename <- input$inputGeneSymbol
    
    #plot
    data <- data.frame()
    for (i in genename){
      asd <- WATFACS_data[i,]
      asd <- data.frame(samples, as.numeric(asd), Gene=i)
      colnames(asd) <- c(colnames(samples), "data", "Gene")
      data <- rbind(data, asd)
    }
    
    # gg <- ggplot(data, aes(x=CellType, y=data, fill=Weight)) +  
    #   geom_boxplot()  + theme_bw() + theme +
    #   facet_wrap(~Gene, scales="free_y") +
    #   labs(x=element_blank(),
    #        y="mRNA expression, log2") +
    #   scale_fill_manual(values=c("gray", "darkred"))

    gg <- ggplot(data, aes(x = CellType, y = data, fill = group)) +  
      geom_boxplot() + theme_bw(16) + 
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
            legend.title = element_blank()) +
      facet_wrap(~Gene, scales="free_y") +
      labs(x=element_blank(),
           y="mRNA expression, log2")
    gg
    
  })

  output$genePlot <- renderPlot({
    plotData()
  })
  

}


# Run the app ----
shinyApp(ui = ui, server = server)
