#-----------------------------------------------------------------------
#
# Invenotry of the TaqMan probes
#
#----------------------------------------------------------------------
library(DT)
Egg_all <- readRDS("full_inventory.Rds")
last_update <- "September 28th, 2023"

# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Liquid Nitrogen stocks"),
                         "Cells and biopsies currently available in the lab as of",
                         tags$br(), tags$b(last_update)
                ),
                fluidRow(style="padding:1% 2% 1% 4%;text-align:left;font-size: 100%",
                         
                         tags$u("General rules:"),
                         tags$li("Personal stocks of cells must be in Egg14. Ask Nico if you need more space."),
                         tags$li("Do not take cells at low passages (P2 or less)"),
                         tags$li("If you take cells at P3 always freeze down your own stocks at P5+"),
                         tags$li("After you take a cryovial from the stocks, delete the row in the excel file on the server."),
                         tags$hr(),
                ),
                fluidRow(style="padding:1% 2% 1% 2%;text-align:center;font-size: 90%",
                         "Search and select the samples you want in the first tab.",
                         "All selected samples are available in the second tab",
                         tabsetPanel(type = "tabs", 
                                     tabPanel("All samples", style="padding:1% 2% 1% 2%;text-align:left",
                                              dataTableOutput("LNtable")
                                     ),
                                     tabPanel("Selected samples", style="padding:1% 2% 1% 2%;text-align:left",
                                              dataTableOutput("selected_samples"),
                                              downloadButton('downloadData', 'Download Sample Location')
                                     )
                         )
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  output$LNtable <- renderDataTable({
    alldata <- Egg_all[,1:11]
    datatable(alldata, 
              options = list(lengthMenu = c(10, 50, 100), 
                             pageLength = 10,
                             autoWidth = T),
              filter = 'bottom',
              rownames= FALSE)
  })
  
  selected_samples <- reactive({
    selected_samples <- Egg_all[input$LNtable_rows_selected,]
    return(selected_samples)
  })
  
  # Make file output with selected rows
  output$selected_samples <- renderDataTable({ 
    selected_samples <- selected_samples()
    datatable(selected_samples[,1:11], 
              options = list(lengthMenu = c(10, 50, 100), 
                             pageLength = 10),
              rownames= FALSE)
  })
  
  
  #button to download file with locations
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('SelectedSamples-', Sys.Date(), '.xlsx', sep='')
    },
    content = function(con) {
      write.xlsx(Egg_all[input$LNtable_rows_selected, c(1,4,2,5:11)], con, row.names=F)
    }
  )
  
}

# Run the app ----
shinyApp(ui = ui, server = server)

