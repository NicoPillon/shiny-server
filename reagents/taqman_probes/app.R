#-----------------------------------------------------------------------
#
# Inventory of the TaqMan probes
#
#----------------------------------------------------------------------
library(DT)
library(openxlsx)
full_inventory <- readRDS("full_inventory.Rds")
last_update <- readRDS("last_update.Rds")

#=====================================================================================================================
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
                         h3("TaqMan Probes Inventory"),
                         "Last update:", last_update
                ),
                fluidRow(style="padding:1% 2% 1% 4%;text-align:left;font-size: 100%",
                         
                  tags$u("General rules:"),
                  tags$li("Whenever possible, use SYBRgreen oligos instead of TaqMan assays."),
                  tags$li("Always search the inventory before ordering new probes."),
                  tags$li("Talk to", 
                          a("Mutsumi Katayama", href="mailto:mutsumi.katayama@ki.se"),
                          "if you want your probes to be transfered to the common boxes."),
                  tags$li("If a probe cannot be found anywhere, tell Mutsumi and/or update the inventory."),
                  tags$hr(),
                  
                  tags$u("Probes are stored in 2 main locations:"),
                  tags$li("Common boxes. If a probe is labelled with a box number and location, it is in one of the common boxes."),
                  tags$li("Personal boxes. If a probe only has a person's name, it is likely in a personal box. Ask that person!\n")
                ),
                tags$hr(),
                fluidRow(style="padding:1% 2% 1% 2%;text-align:center;font-size: 90%",
                         "Search and select the probes you want in the first tab.",
                         "Selected probes are available to download in the second tab.",
                         tabsetPanel(type = "tabs",
                                     tabPanel("Probes", style="padding:1% 2% 1% 2%;text-align:left",
                                              dataTableOutput("LNtable")
                                     ),
                                     tabPanel("Download Location", style="padding:1% 2% 1% 2%;text-align:left",
                                              dataTableOutput("selected_samples"),
                                              downloadButton('downloadData', 'Download Probe Locations')
                                     )
                         )
                )
)

# Define server logic ----
server <- function(input, output, session) {
  
  output$LNtable <- renderDataTable({
    alldata <- full_inventory
    datatable(alldata, 
              options = list(lengthMenu = c(10, 50, 100), 
                             pageLength = 10,
                             autoWidth = T),
              filter = 'bottom',
              rownames= FALSE,
              style = "bootstrap")
  })
  
  selected_samples <- reactive({
    selected_samples <- full_inventory[input$LNtable_rows_selected,]
    return(selected_samples)
  })
  
  # Make file output with selected rows
  output$selected_samples <- renderDT({ 
    selected_samples <- selected_samples()
    datatable(selected_samples, 
              options = list(lengthMenu = c(10, 50, 100), 
                             pageLength = 10),
              rownames= FALSE)
  })
  
  
  #button to download file with locations
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('SelectedProbes-', Sys.Date(), '.xlsx', sep='')
    },
    content = function(con) {
      asd <- full_inventory[input$LNtable_rows_selected, ]
      asd[nrow(asd)+2,1] <- "CANNOT FIND A PROBE?"
      asd[nrow(asd)+1,1] <- "Write it down and ask Mutsumi to update the inventory!"
      write.xlsx(asd, con, row.names=F)
    }
  )
  
}

# Run the app ----
shinyApp(ui = ui, server = server)

