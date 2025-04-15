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
    
    VennPlot_function <- function() {
      MainField <- isolate(input$input_MainField)
      SubFields <- list(A = isolate(input$input_SubFieldA),
                        B = isolate(input$input_SubFieldB),
                        C = isolate(input$input_SubFieldC))
      SubFields <- SubFields[SubFields != ""] # remove empty fields
      field_names <- names(SubFields)
      
      withProgress(message = 'Collecting data from PubMed', value = 1, max = length(SubFields)^2 + 2, {
        
        count_AllResearch <- 31821468
        incProgress(1)
        
        count_MainField <- EUtilsSummary(MainField, type="esearch", db="pubmed")@count
        incProgress(1)
        
        combo_counts <- list()
        for (field in field_names) {
          query <- paste(MainField, SubFields[[field]], sep = " AND ")
          combo_counts[[field]] <- EUtilsSummary(query, type="esearch", db="pubmed")@count
          incProgress(1)
        }
        
        if (length(field_names) >= 2) {
          combn_pairs <- combn(field_names, 2, simplify = FALSE)
          for (pair in combn_pairs) {
            query <- paste(MainField, SubFields[[pair[1]]], SubFields[[pair[2]]], sep = " AND ")
            combo_counts[[paste(pair, collapse = "")]] <- EUtilsSummary(query, type="esearch", db="pubmed")@count
            incProgress(1)
          }
        }
        
        if (length(field_names) == 3) {
          query <- paste(MainField, SubFields$A, SubFields$B, SubFields$C, sep = " AND ")
          combo_counts[["ABC"]] <- EUtilsSummary(query, type="esearch", db="pubmed")@count
          incProgress(1)
        }
        
        # Convert to Venn diagram counts (artificial scaling for better layout)
        euler_input <- c(`Total` = (count_AllResearch)^(1/1.2),
                         `Total&Main` = (count_MainField)^(1/1.2))
        
        for (key in names(combo_counts)) {
          label <- paste("Total&Main", paste(strsplit(key, "")[[1]], collapse = "&"), sep = "&")
          euler_input[[label]] <- (combo_counts[[key]])^(1/1.2)
        }
        
        percent_of_total <- if ("ABC" %in% names(combo_counts)) combo_counts$ABC * 100 / count_AllResearch else NA
        
        vennplot <- plot(euler(euler_input, shape = "ellipse"),
                         main = if (!is.na(percent_of_total)) paste("Your selected field of research represents ", signif(percent_of_total, 2), "% of all publications", sep = "") else "Field intersections",
                         alpha = 0.5,
                         edges = list(col = "black", lex = 1),
                         fill = RColorBrewer::brewer.pal(length(euler_input), "Set3"))
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