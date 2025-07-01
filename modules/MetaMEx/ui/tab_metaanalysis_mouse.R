tabMetaanalysis_mouse <- tabPanel("Mouse", value="panelAppMouse", 
                                  
                                  # Row introduction ####################################################################################
                                  fluidRow(style="margin:-18px -15px 0px -15px;padding:0.5% 1.5% 0% 1.3%", # need to offset the margins
                                           h5("Use the meta-analysis to obtain an overview of the mRNA expression of your gene of interest in all available exercise and inactivity studies."
                                           ),
                                           selectizeInput("genename_metaanalysis_mouse", label=NULL, 
                                                          choices=NULL, 
                                                          selected=NULL, options=NULL)
                                  ),
                                  
                                  # Row App and plots ###################################################################################
                                  fluidRow(style="padding:0% 2% 1% 2%;",
                                           tabsetPanel(type="pills", id="inTabsetMeta", 
                                                       
                                                       #= Tab overview =====================================================        
                                                       tabPanel("Overview", 
                                                                value="panelAppOverviewMouse",
                                                                fluidRow(style="background-color: white", 
                                                                         column(1),
                                                                         column(8, style="padding:2% 5% 2% 2%;margin:0",
                                                                                "This plot provides an overview of the meta-analysis for the selected gene.
                                                                                       The plot is reactive and updated based on the selection criteria used in all individual forest plots.
                                                                                       The individual forest plots uderlying this figure can be found in the different tabs.",
                                                                                tags$br(), tags$br(),
                                                                                plotOutput("mouse_overviewPlot") %>% withSpinner(color="#5b768e"),
                                                                                tags$br(),
                                                                                tags$i("The meta-analysis statistics for the selected gene and criteria can be downloaded here:"),
                                                                                tags$br(),
                                                                                downloadButton("mouse_downloadOverview", label = "Download selected data")
                                                                         )
                                                                )
                                                       ),
                                                       
                                                       #= Tab Acute Exercise =====================================================        
                                                       tabPanel("Acute Exercise", value="panelAppAcute",
                                                                fluidRow(style="background-color:white;padding:1% 1% 1% 1%",
                                                                         column(10, style="padding:0% 5% 2% 2%;margin:0",
                                                                                uiOutput("mouse_forestPlotDynamic_AA"),
                                                                                hr(),
                                                                                checkboxGroupInput("mouse_studies_AA", 
                                                                                                   tags$b("Include/exclude specific datasets"), 
                                                                                                   selected=mouse_input[['GEO_AA']],
                                                                                                   mouse_input[['GEO_AA']],
                                                                                                   inline = TRUE)
                                                                         ),
                                                                         column(2, style="background-color:#eaf1f7;margin-top:60px;padding:0 1% 1% 2%",
                                                                                h3("Filter by", style="color:#011b2a"),
                                                                                sliderInput("mouse_acute_biopsy", tags$b("Time after exercise (h)"),
                                                                                            min = 0, max = 48, value = c(0,48), step = 1, sep = ""),
                                                                                tags$br(),
                                                                                checkboxGroupInput("mouse_exercise_type", tags$b("Exercise type"), 
                                                                                                   choices=mouse_input[['exercise_type_choice']],
                                                                                                   selected=mouse_input[['exercise_type_choice']]),
                                                                                checkboxInput('mouse_exercise_type_allnone', 'All/None', value=T)

                                                                         ),
                                                                ),
                                                       ),
                                                       
                                                       #= Tab Inactivity =====================================================          
                                                       tabPanel("Inactivity", value="panelAppInactivity",
                                                                fluidRow(style="background-color:white;padding:1% 1% 1% 1%",
                                                                         column(10, style="padding:0% 5% 2% 2%;margin:0",
                                                                                uiOutput("mouse_forestPlotDynamic_IN"),
                                                                                hr(), 
                                                                                checkboxGroupInput("mouse_studies_IN", 
                                                                                                   tags$b("Include/exclude specific datasets"), 
                                                                                                   selected=mouse_input[['GEO_IN']], 
                                                                                                   mouse_input[['GEO_IN']],
                                                                                                   inline = TRUE)
                                                                         ),
                                                                         column(2, style="background-color:#eaf1f7;margin-top:60px;padding:0 1% 1% 2%",
                                                                                h3("Filter by", style="color:#011b2a"),
                                                                                sliderInput("mouse_inactivity_duration", tags$b("Duration (days)"),
                                                                                            min = 0, max = 30, value = c(0,30), step = 1, sep = ""),
                                                                                tags$br(),
                                                                                tags$br(),
                                                                                checkboxGroupInput("mouse_inactivity_protocol", 
                                                                                                   tags$b("Protocol"), 
                                                                                                   choices=mouse_input[['inactivity_protocol_choice']],
                                                                                                   selected=mouse_input[['inactivity_protocol_choice']]),
                                                                                checkboxInput('mouse_inactivity_protocol_allnone', 'All/None', value=T)

                                                                         )
                                                                )
                                                       ),
                                                       
                                                       #= Tab Exercise Training =====================================================       
                                                       tabPanel("Exercise Training", value="panelAppTraining",
                                                                fluidRow(style="background-color:white;padding:1% 1% 1% 1%",
                                                                         column(10, style="padding:0% 5% 2% 2%;margin:0",
                                                                                uiOutput("mouse_forestPlotDynamic_TA"),
                                                                                hr(),
                                                                                checkboxGroupInput("mouse_studies_TA", 
                                                                                                   tags$b("Include/exclude specific datasets"),
                                                                                                   selected=mouse_input[['GEO_TA']],
                                                                                                   mouse_input[['GEO_TA']],
                                                                                                   inline = T)
                                                                         ),
                                                                         column(2, style="background-color:#eaf1f7;margin-top:60px;padding:0 1% 1% 2%",
                                                                                h3("Filter by", style="color:#011b2a"),
                                                                                sliderInput("mouse_training_duration", tags$b("Duration (weeks)"),
                                                                                            min = 0, max = 24, value = c(0,24), step = 1, sep = ""),
                                                                                tags$br(),
                                                                                checkboxGroupInput("mouse_training_protocol", 
                                                                                                   tags$b("Protocol"), 
                                                                                                   choices=mouse_input[['training_protocol_choice']],
                                                                                                   selected=mouse_input[['training_protocol_choice']]),
                                                                                checkboxInput('mouse_training_protocol_allnone', 'All/None', value=T)
                                                                         )
                                                                )
                                                       )
                                                       

                                           )
                                           
                                  ),
                                  
                                  # Row Select Population ###################################################################################
                                  tags$br(),
                                  fluidRow(style="position:center;background-color:#f4eae7;padding:0 0 1% 2%;margin:0 6% 5% 5%",
                                           h3("Customize your mouse population of interest", style="color:#c93f1e"),
                                           column(3, checkboxGroupInput("mouse_sex", tags$b("Sex"),
                                                                        choices=mouse_input[['sex_choice']],
                                                                        selected=mouse_input[['sex_choice']]), #checkbox to select category
                                                  # checkboxInput('mouse_sex_allnone', 'All/None', value=T),
                                                  hr(),
                                                  sliderInput("mouse_age", tags$b("Age (weeks)"),
                                                              min = 4, max = 24, value = c(4,24), step = 1, sep = "")
                                           ), 
                                           column(3, checkboxGroupInput("mouse_muscle", tags$b("Muscle"), 
                                                                        choices=mouse_input[['muscle_choice']],
                                                                        selected=mouse_input[['muscle_choice']]
                                           ),
                                           checkboxInput('mouse_muscle_allnone', 'All/None', value=T)
                                           ), 
                                           column(4, checkboxGroupInput("mouse_disease", tags$b("Mouse model"), 
                                                                        choices=mouse_input[['disease_choice']],
                                                                        selected=mouse_input[['disease_choice']]), #checkbox to select category
                                                  checkboxInput('mouse_disease_allnone', 'All/None', value=T)
                                           )
                                  )
)
