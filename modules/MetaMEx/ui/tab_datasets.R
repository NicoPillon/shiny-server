tabDatasets <- navbarMenu("Datasets",
                          ########################################################################################
                          tabPanel("Human Studies", 
                                   fluidRow(style="background-color: #5b768e;margin:-18px -15px 0px -15px;padding:0.5% 1.5% 0.5% 1.3%",
                                            h5(style="color:white","Description of the studies.")
                                   ),
                                   fluidRow(style="font-size: 80%;padding:0% 1% 10% 1%",
                                            DT::dataTableOutput("reftable_human"))

                          ),
                          tabPanel("Mouse Studies", 
                                   fluidRow(style="background-color: #5b768e;margin:-18px -15px 0px -15px;padding:0.5% 1.5% 0.5% 1.3%",
                                            h5(style="color:white","Description of the studies.")
                                   ),
                                   fluidRow(style="font-size: 80%;padding:0% 1% 10% 1%",
                                            DT::dataTableOutput("reftable_mouse"))
                          )
)