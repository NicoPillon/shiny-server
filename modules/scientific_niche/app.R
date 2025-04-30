#-----------------------------------------------------------------------
#
# My field of research
#
#----------------------------------------------------------------------
# Load data and libraries
library(tidyverse)
library(shinycssloaders)
library(RISmed)
library(eulerr)

# Define UI ----
ui <- fluidPage(

  # CSS for style
  tags$head(includeCSS("../../www/style.css")),
  
  # title ribbon
  fluidRow(
    style = "color:black; padding:0% 2% 0% 2%; text-align:left;",
    h3("A Small Part of Something Big - Explore Your Scientific Niche"),
    h5("Last update 2025-04-15"),
    tags$hr()
  ),
  
  # Main layout with plot and controls
  sidebarLayout(
    sidebarPanel(h6(textInput("input_MainField", "Main field of research", value = '(obesity OR diabetes)')),
                 h6(textInput("input_SubFieldA", "Keyword A", value = 'inflammation')),
                 h6(textInput("input_SubFieldB", "Keyword B", value = '"skeletal muscle"')),
                 h6(textInput("input_SubFieldC", "Keyword C", value = 'exercise')),
                 tags$hr(),
                 actionButton("do", "Update Figure", icon("refresh")),
                 tags$hr(),
                 downloadButton("downloadPlot", "Download plot (.svg)")
    ),
    
    mainPanel(tags$br(),
              plotOutput("VennOutput", height = "500px") %>% withSpinner(color = "#5B768E", type = 8),
              div(textOutput("percentText"), style = "text-align: center;")
    )
  ),
  
  # Description of methods
  fluidRow(
    style = "color:black; background-color:white; padding:0% 2% 0% 2%;",
    tags$hr(),
    h3("Methods"),
    p("This application visualizes how a specific scientific niche fits into the broader biomedical research landscape."),
    p("Using PubMed as a data source, it queries the number of publications for a main research field and up to three user-defined keywords. It then computes intersections between these keywords and displays the results as a Venn diagram."),
    p("Each segment of the diagram represents a subset of the literature, scaled proportionally for readability. The percentage of total publications represented by the selected niche is calculated and displayed as a gentle reminder of the vast scope of scientific contributions."),
    p("The tool is intended to provide perspective — highlighting both the relevance of a specific area and its place within the collective body of research.")
  )
  
)


# Define server logic ----
server <- function(input, output) {
  
  venn_result <- eventReactive(input$do, {
    MainField <- input$input_MainField
    SubFields <- list(A = input$input_SubFieldA,
                      B = input$input_SubFieldB,
                      C = input$input_SubFieldC)
    SubFields <- SubFields[SubFields != ""]
    SubFields_names <- names(SubFields)
    
    withProgress(message = 'Collecting data from PubMed', value = 1, max = length(SubFields)^2 + 2, {
      
      # total number of publications on PubMed
      AllResearch_count <- 31821468
      incProgress(1)
      
      # Total number of publications in the main field of study selected
      MainField_count <- EUtilsSummary(MainField, type="esearch", db="pubmed")@count
      incProgress(1)
      
      # Number of publications corresponding to keywords A-B-C
      combo_counts <- list()
      for (field in SubFields_names) {
        query <- paste(MainField, SubFields[[field]], sep = " AND ")
        combo_counts[[field]] <- EUtilsSummary(query, type="esearch", db="pubmed")@count
        #incProgress(1)
      }
      
      # number of publications intersecting keywords A-B-C
      if (length(SubFields_names) > 1) {
        combn_pairs <- combn(SubFields_names, 2, simplify = FALSE)
        for (pair in combn_pairs) {
          query <- paste(MainField, SubFields[[pair[1]]], SubFields[[pair[2]]], sep = " AND ")
          combo_counts[[paste(pair, collapse = "")]] <- EUtilsSummary(query, type="esearch", db="pubmed")@count
          #incProgress(1)
        }
      }
      
      # number of publications interecting all keywords
      if (length(SubFields_names) > 2) {
        query <- paste(c(MainField, SubFields), collapse = " AND ")
        combo_counts[["ABC"]] <- EUtilsSummary(query, type = "esearch", db = "pubmed")@count
        #incProgress(1)
      }
      
      # scale imput
      euler_input <- c(`Total` = (AllResearch_count)^(1/1.2),
                       `Total&Main` = (MainField_count)^(1/1.2))
      
      # labels for plot
      for (key in names(combo_counts)) {
        label <- paste("Total&Main", paste(strsplit(key, "")[[1]], collapse = "&"), sep = "&")
        euler_input[[label]] <- (combo_counts[[key]])^(1/1.2)
      }
      
      # calculate the percentage of the intersect compared to total number of publications
      # Create a key from all selected subfield names (e.g., "A", "AB", "ABC")
      combo_key <- if (length(SubFields_names) >= 1) paste(SubFields_names, collapse = "") else NULL
            # Calculate the percent if the key exists
      if (!is.null(combo_key) && combo_key %in% names(combo_counts)) {
        percent_of_total <- combo_counts[[combo_key]] * 100 / AllResearch_count
      } else {
        percent_of_total <- NA
      }
      
      # plot
      vennplot <- plot(euler(euler_input, shape = "ellipse"),
                       main = NULL,
                       alpha = 0.5,
                       edges = list(col = "black", lex = 1),
                       fill = RColorBrewer::brewer.pal(length(euler_input), "Set3"))
      
      return(list(plot = vennplot, percent = percent_of_total))
    })
  })
  
  output$VennOutput <- renderPlot({
    venn_result()[["plot"]]
  })
  
  output$percentText <- renderText({
    result <- venn_result()
    if (!is.null(result$percent)) {
      paste0("Your selected field of research represents ", signif(result$percent, 2), "% of all publications.")
    } else {
      ""
    }
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { 'MyResearchField_NicoPillon.com.svg' },
    content = function(file) {
      svg(file, width = 8, height = 8)
      print(venn_result()[["plot"]])
      dev.off()
    }
  )
}



# Run the app ----
shinyApp(ui = ui, server = server)