#-----------------------------------------------------------------------
#
# Inventory of the TaqMan probes
#
#----------------------------------------------------------------------
library(DT)
full_inventory <- readRDS("full_inventory.Rds")
last_update <- "September 28th, 2023"

#=====================================================================================================================
# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("TaqMan Probes Inventory"),
                         "Search probes currently available in the lab as of",
                         tags$br(), tags$b(last_update)
                ),
                fluidRow(style="padding:1% 2% 1% 4%;text-align:left;font-size: 100%",
                         
                  tags$u("General rules:"),
                  tags$li("Whenever possible, use SYBRgreen oligos instead of TaqMan assays."),
                  tags$li("Always search the inventory before ordering new probes."),
                  tags$li("Talk to", 
                          a("Mutsumi Katayama", href="mailto:mutsumi.katayama@ki.se"),
                          "if you want your probes to be transfered to the common boxes."),
                  tags$hr(),
                  
                  tags$u("Probes are stored in 2 main locations:"),
                  tags$li("Common boxes. If a probe is labelled with a box number and location, it is in one of the common boxes."),
                  tags$li("Personal boxes. If a probe only has a person's name, it is likely in a personal box. Ask that person!\n")
                ),
                tags$hr(),
                fluidRow(style="padding:1% 2% 1% 2%;text-align:center;font-size: 90%",
                         "Search and select the probes you want in the first tab.",
                         "Selected probes are available to donwload in the second tab.",
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
              rownames= FALSE)
  })
  
  selected_samples <- reactive({
    selected_samples <- full_inventory[input$LNtable_rows_selected,]
    return(selected_samples)
  })
  
  # Make file output with selected rows
  output$selected_samples <- renderDataTable({ 
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
      write.xlsx(full_inventory[input$LNtable_rows_selected, ], con, row.names=F)
    }
  )
  
}

# Run the app ----
shinyApp(ui = ui, server = server)

