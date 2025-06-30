#-----------------------------------------------------------------------
#
# Invenotry of the TaqMan probes
#
#----------------------------------------------------------------------
library(DT)
library(openxlsx)
Egg_all <- readRDS("full_inventory.Rds")
last_update <- readRDS("last_update.Rds")

# Define UI ----
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      /* Change the colors of hyperlinks */
      a {color: #337AB7}
 
      /* Change the colors of selected rows in data table */
      .table.dataTable tbody td.active, .table.dataTable tbody tr.active td {
      background-color: #337AB7 !important;
      color: white !important;
                  "))
  ),
  
                fluidRow(style="color:white;background-color:#337AB7;padding:0% 1% 1% 1%;text-align:center",
                         h3("Liquid Nitrogen stocks"),
                         "Last update:", last_update
                ),
                fluidRow(style="padding:1% 2% 1% 4%;text-align:left;font-size: 100%",
                         
                         tags$u("General rules:"),
                         tags$li("Personal stocks of cells must be in Egg14. Ask Nico if you need more space."),
                         tags$li("Never take cells at passage P2 or less.", tags$b("Ask Katrin or Nico first!")),
                         tags$li("If you take cells at P3 always freeze down your own stocks at P5+"),
                         tags$li("Always", tags$b("delete the row"), "in the excel file on the server after you take a vial!"),
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
              rownames= FALSE,
              style = "default")
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
      asd <- Egg_all[input$LNtable_rows_selected, c(1,4,2,5:11)]
      asd[nrow(asd)+2,1] <- "REMEMBER TO DELETE THE VIALS YOU PICKED IN THE EXCEL FILE!"
      asd[nrow(asd)+1,1] <- "P:/C3_Integrative_Physiology_Group/Liquid Nitrogen/LN inventory.xlsx"
      write.xlsx(asd, con, row.names=F)
      #write.xlsx(Egg_all[input$LNtable_rows_selected, c(1,4,2,5:11)], con, row.names=F)
    }
  )
  
}

# Run the app ----
shinyApp(ui = ui, server = server)


