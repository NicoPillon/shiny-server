########################################################################################
tabDownloads <- tabPanel("Download", 
                               ########################################################################################
                               fluidRow(style="background-color: #5b768e;margin:-18px -15px 0px -15px;padding:0.5% 1.5% 0.5% 1.3%",
                                        h5(style="color:white","Here are all the raw data from MetaMEx. For more details, feel free to contact us.")
                               ),
                               fluidRow(style="padding:0% 5% 10% 5%",
                                 column(4, 
                                        tags$br(),
                                        h2("Human"),
                                        tags$br(),
                                        downloadButton("download_human_metadata", label = "Human - Metadata"),
                                        tags$br(),
                                        tags$br(),
                                        downloadButton("download_human_genelist", label = "Human - Gene list"),
                                        tags$br(),
                                        tags$br(),
                                        downloadButton("download_human_matrix", label = "Human - Gene matrix")
                                        ),
                                 column(4,
                                        tags$br(),
                                        h2("Mouse"),
                                        tags$br(),
                                        downloadButton("download_mouse_metadata", label = "Mouse - Metadata"),
                                        tags$br(),
                                        tags$br(),
                                        downloadButton("download_mouse_genelist", label = "Mouse - Gene list"),
                                        tags$br(),
                                        tags$br(),
                                        downloadButton("download_mouse_matrix", label = "Mouse - Gene matrix")
                                        )
                               )
)