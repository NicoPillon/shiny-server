#-----------------------------------------------------------------------
#
# My field of research
#
#----------------------------------------------------------------------
# Load data and libraries
# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                tags$head(
                  tags$style(
                    HTML(".shiny-notification {
              height: 50px;
              width: 800px;
              position:fixed;
              top: calc(50% - 50px);;
              left: calc(50% - 400px);;
            }
           "
                    )
                  )
                ),
                sidebarLayout(
                  sidebarPanel(h4("Keywords"),
                               h6(textInput("input_MainField", "Main field of research", value = '(obesity OR diabetes)')),
                               h6(textInput("input_SubFieldA", "Keyword A", value = 'inflammation')),
                               h6(textInput("input_SubFieldB", "Keyword B", value = '"skeletal muscle"')),
                               h6(textInput("input_SubFieldC", "Keyword C", value = 'exercise')),
                               tags$i(h6("All fields are required")),
                               actionButton("do", "Update Figure", icon("refresh")),
                               tags$hr(),
                               h4("Download vector image"),
                               tags$i(h6("An svg image file can be edited in powerpoint or illustrator. 
                                  In powerpoint, import the image and 'ungroup' it to be able to 
                                  edit colors and text")),
                               downloadButton("downloadPlot", "Download plot (.svg)")
                  ),
                  
                  mainPanel(tags$br(),
                            plotOutput("VennOutput", height = "500px")
                  )
                )
)


# Define server logic ----
server <- function(input, output) {
  #install.packages("RISmed")
  library(RISmed)
  library(eulerr)
  
  observeEvent(input$do, {
    
    VennPlot_function <- function(){
      MainField <- isolate(input$input_MainField)
      SubFieldA <- isolate(input$input_SubFieldA)
      SubFieldB <- isolate(input$input_SubFieldB)
      SubFieldC <- isolate(input$input_SubFieldC)
      
      withProgress(message = 'Collecting data from PubMed', value = 1, max=10, {
        #all publications on pubmed
        count_AllResearch <- 31821468 #EUtilsSummary("all[sb]", type="esearch", db="pubmed")@count
        Sys.sleep(0.2) 
        incProgress(1)
        
        count_MainField <- EUtilsSummary(MainField, type="esearch", db="pubmed")@count
        Sys.sleep(0.2) 
        incProgress(1)
        
        count_SubFieldA <- EUtilsSummary(paste(MainField, SubFieldA, sep=" AND "), type="esearch", db="pubmed")@count
        Sys.sleep(0.2) 
        incProgress(1)
        
        count_SubFieldB <- EUtilsSummary(paste(MainField, SubFieldB, sep=" AND "), type="esearch", db="pubmed")@count
        Sys.sleep(0.2) 
        incProgress(1)
        
        count_SubFieldAB <- EUtilsSummary(paste(MainField, SubFieldA, SubFieldB, sep=" AND "), type="esearch", db="pubmed")@count
        Sys.sleep(0.2) 
        incProgress(1)
        
        count_SubFieldC <- EUtilsSummary(paste(MainField, SubFieldC, sep=" AND "), type="esearch", db="pubmed")@count
        Sys.sleep(0.2) 
        incProgress(1)
        
        count_SubFieldAC <- EUtilsSummary(paste(MainField, SubFieldA, SubFieldC, sep=" AND "), type="esearch", db="pubmed")@count
        Sys.sleep(0.2) 
        incProgress(1)
        
        count_SubFieldBC <- EUtilsSummary(paste(MainField, SubFieldB, SubFieldC, sep=" AND "), type="esearch", db="pubmed")@count
        Sys.sleep(0.2) 
        incProgress(1)
        
        count_SubFieldABC <- EUtilsSummary(paste(MainField, SubFieldA, SubFieldB, SubFieldC, sep=" AND "), type="esearch", db="pubmed")@count
        incProgress(1)
        
        percent_of_total <- count_SubFieldABC*100 / count_AllResearch
        
        VennDiag <- euler(c("Total" = (count_AllResearch)^(1/1.2),
                            "Total&Main" = (count_MainField)^(1/1.2),
                            "Total&Main&A" = (count_SubFieldA)^(1/1.2),
                            "Total&Main&B" = (count_SubFieldB)^(1/1.2), 
                            "Total&Main&C" = (count_SubFieldC)^(1/1.2),
                            "Total&Main&A&B" = (count_SubFieldAB)^(1/1.2), 
                            "Total&Main&A&C" = (count_SubFieldAC)^(1/1.2),
                            "Total&Main&B&C" = (count_SubFieldBC)^(1/1.2),
                            "Total&Main&A&B&C" = (count_SubFieldABC)^(1/1.2)),
                          shape = "ellipse")
        Fields <- c("All publications in biology and medicine", MainField, SubFieldA, SubFieldB, SubFieldC)
        Fields <- gsub("\\(|\\)", "", Fields)
        Fields <- gsub('\\"', "", Fields)
        Fields <- gsub('OR', "&", Fields)
        Fields <- gsub(' ', "\n", Fields)
        
        vennplot <- plot(VennDiag, 
                         main=paste("My field of research represents ", signif(percent_of_total, 2), "% of all publications", sep=""),
                         alpha=0.5,
                         edges = list(col=c("black", "#0072B2", "#E69F00", "#009E73", "black", "grey", "grey", "gray", "yellow"),
                                      lex=1),
                         fill=c("white", "#56B4E9", "#E69F00", "#009E73", "black", "#739f3a", "#d98c54", "#668c8d", "yellow"),
                         labels = list(font = 1,
                                       pos=4,
                                       cex = c(1,1,1,1,1),
                                       cex.main=0.5))
        return(vennplot)
      })
      
    }
    
    output$VennOutput <- renderPlot({
      VennPlot_function()
    })
    
    
    output$downloadPlot <- downloadHandler(
      filename = function() { 'MyResearchField_NicoPillon.com.svg' },
      content = function(file) {
        svg(file, width=8, height=8)
        print(VennPlot_function())
        dev.off()
      })
  })
}


# Run the app ----
shinyApp(ui = ui, server = server)