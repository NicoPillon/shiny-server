tabMetaanalysis_human <- tabPanel("Human", value="panelAppHuman", 
                                  
                                  tags$script(HTML("
  $('a[data-value=\"panelAppAcute\"]').on('shown.bs.tab', function (e) {
    $('#human_acute_biopsy').data('ionRangeSlider').reset();
  });
")),
                                  tags$script(HTML("
  $('a[data-value=\"panelAppAcute\"]').on('shown.bs.tab', function (e) {
    $('#human_acute_biopsy').data('ionRangeSlider').update();
  });
")),
                                  
                                  # Row introduction ####################################################################################
                                  fluidRow(style="margin:-18px -15px 0px -15px;padding:0.5% 1.5% 0% 1.3%", # need to offset the margins
                                           h5("Use the meta-analysis to obtain an overview of the mRNA expression of your gene of interest in all available exercise and inactivity studies."
                                           ),
                                           selectizeInput("genename_metaanalysis_human", label=NULL, 
                                                          choices=NULL, 
                                                          selected=NULL, options=NULL)
                                  ),
                                  
                                  # Row App and plots ###################################################################################
                                  fluidRow(style="padding:0% 2% 1% 2%;",
                                           tabsetPanel(type="pills", id="inTabsetMeta", 
                                                       
                                                       #= Tab overview =====================================================        
                                                       tabPanel("Overview", 
                                                                value="panelAppOverviewHuman",
                                                                fluidRow(style="background-color: white", 
                                                                         column(1),
                                                                         column(10, style="padding:2% 5% 2% 2%;margin:0",
                                                                                "This plot provides an overview of the meta-analysis for the selected gene.
                                                                                       The plot is reactive and updated based on the selection criteria used in all individual forest plots.
                                                                                       The individual forest plots uderlying this figure can be found in the different tabs.",
                                                                                tags$br(), tags$br(), 
                                                                                plotOutput("human_overviewPlot") %>% withSpinner(color="#5b768e"),
                                                                                tags$br(),
                                                                                tags$i("The meta-analysis statistics for the selected gene and criteria can be downloaded here:"),
                                                                                tags$br(),
                                                                                downloadButton("human_downloadOverview", label = "Download selected data")
                                                                         )
                                                                )
                                                       ),
                                                       
                                                       #= Tab Acute Exercise =====================================================        
                                                       navbarMenu("Acute Exercise",
                                                                  
                                                                  tabPanel("Aerobic",
                                                                           fluidRow(style="background-color:white;padding:1% 1% 1% 1%",
                                                                                    column(10, style="padding:0% 5% 2% 2%;margin:0",
                                                                                           uiOutput("human_forestPlotDynamic_AA") %>% withSpinner(color="#5b768e"),
                                                                                           hr(),
                                                                                           checkboxGroupInput("human_studies_AA", 
                                                                                                              tags$b("Include/exclude specific datasets"), 
                                                                                                              selected=human_input[['GEO_AA']],
                                                                                                              choices=human_input[['GEO_AA']],
                                                                                                              inline = TRUE)
                                                                                    ),
                                                                                    column(2, style="background-color:#eaf1f7;margin-top:60px;padding:0 1% 1% 2%",
                                                                                           h3("Filter by", style="color:#011b2a"),
                                                                                           sliderInput("human_acute_biopsy", tags$b("Time after exercise (h)"),
                                                                                                       min = 0, max = 100, value = c(0,4), step = 1, sep = ""),
                                                                                           tags$br(),
                                                                                           checkboxGroupInput("human_exercise_type", tags$b("Exercise type"),
                                                                                                              choices=human_input[['exercise_type_choice']],
                                                                                                              selected=human_input[['exercise_type_choice']]),
                                                                                           checkboxInput('human_exercise_type_allnone', 'All/None', value=TRUE)
                                                                                    )
                                                                           )
                                                                  ),
                                                                  
                                                                  tabPanel("HIT",
                                                                           fluidRow(style="background-color:white;padding:1% 1% 1% 1%",
                                                                                    column(10, style="padding:0% 5% 2% 2%;margin:0",
                                                                                           uiOutput("human_forestPlotDynamic_AH"),
                                                                                           hr(),
                                                                                           checkboxGroupInput("human_studies_AH", 
                                                                                                              tags$b("Include/exclude specific datasets"), 
                                                                                                              selected=human_input[['GEO_AH']],
                                                                                                              choices=human_input[['GEO_AH']],
                                                                                                              inline = TRUE)
                                                                                    ),
                                                                                    column(2, style="background-color:#eaf1f7;margin-top:60px;padding:0 1% 1% 2%",
                                                                                           h3("Filter by", style="color:#011b2a"),
                                                                                           sliderInput("human_acute_biopsy", tags$b("Time after exercise (h)"),
                                                                                                       min = 0, max = 100, value = c(0,4), step = 1, sep = ""),
                                                                                           tags$br(),
                                                                                           checkboxGroupInput("human_exercise_type", tags$b("Exercise type"),
                                                                                                              choices=human_input[['exercise_type_choice']],
                                                                                                              selected=human_input[['exercise_type_choice']]),
                                                                                           checkboxInput('human_exercise_type_allnone', 'All/None', value=TRUE)
                                                                                    )
                                                                           )
                                                                  ),
                                                                  
                                                                  tabPanel("Resistance",
                                                                           fluidRow(style="background-color:white;padding:1% 1% 1% 1%",
                                                                                    column(10, style="padding:0% 5% 2% 2%;margin:0",
                                                                                           uiOutput("human_forestPlotDynamic_AR"),
                                                                                           hr(),
                                                                                           checkboxGroupInput("human_studies_AR", 
                                                                                                              tags$b("Include/exclude specific datasets"), 
                                                                                                              selected=human_input[['GEO_AR']],
                                                                                                              choices=human_input[['GEO_AR']],
                                                                                                              inline = TRUE)
                                                                                    ),
                                                                                    column(2, style="background-color:#eaf1f7;margin-top:60px;padding:0 1% 1% 2%",
                                                                                           h3("Filter by", style="color:#011b2a"),
                                                                                           sliderInput("human_acute_biopsy", tags$b("Time after exercise (h)"),
                                                                                                       min = 0, max = 100, value = c(0,4), step = 1, sep = ""),
                                                                                           tags$br(),
                                                                                           checkboxGroupInput("human_exercise_type", tags$b("Exercise type"),
                                                                                                              choices=human_input[['exercise_type_choice']],
                                                                                                              selected=human_input[['exercise_type_choice']]),
                                                                                           checkboxInput('human_exercise_type_allnone', 'All/None', value=TRUE)
                                                                                    )
                                                                           )
                                                                  ),
                                                                  
                                                                  tabPanel("Timeline",
                                                                           fluidRow(style="background-color:white;padding:1% 1% 1% 1%",
                                                                                    column(12, 
                                                                                           p("The timeline plot is not reactive. All studies in healthy individuals are included by default, and differences (sex, age, protocol) are blocked for in the statistical model."),
                                                                                           hr(),
                                                                                           plotOutput("human_timelinePlot_acute", height="400px"),
                                                                                           hr(),
                                                                                           tableOutput("human_timelineTable_acute")
                                                                                    )
                                                                           )
                                                                  )
                                                       )
                                                       
                                                       ,
                                                       
                                                       #= Tab Inactivity =====================================================          
                                                       navbarMenu("Inactivity",
                                                                  
                                                                  tabPanel("Meta-analysis",
                                                                           fluidRow(style = "background-color:white;padding:1% 1% 1% 1%",
                                                                                    column(10, style = "padding:0% 5% 2% 2%;margin:0",
                                                                                           uiOutput("human_forestPlotDynamic_IN"),
                                                                                           hr(), 
                                                                                           checkboxGroupInput("human_studies_IN", 
                                                                                                              tags$b("Include/exclude specific datasets"), 
                                                                                                              selected = human_input[['GEO_IN']], 
                                                                                                              choices = human_input[['GEO_IN']],
                                                                                                              inline = TRUE)
                                                                                    ),
                                                                                    column(2, style = "background-color:#eaf1f7;margin-top:60px;padding:0 1% 1% 2%",
                                                                                           h3("Filter by", style = "color:#011b2a"),
                                                                                           sliderInput("human_inactivity_duration", tags$b("Duration (days)"),
                                                                                                       min = 0, max = 100, value = c(0, 100), step = 1, sep = ""),
                                                                                           tags$br(),
                                                                                           checkboxGroupInput("human_inactivity_protocol", 
                                                                                                              tags$b("Protocol"), 
                                                                                                              choices = human_input[['inactivity_protocol_choice']],
                                                                                                              selected = human_input[['inactivity_protocol_choice']]),
                                                                                           checkboxInput('human_inactivity_protocol_allnone', 'All/None', value = TRUE)
                                                                                    )
                                                                           )
                                                                  ),
                                                                  
                                                                  tabPanel("Timeline",
                                                                           fluidRow(style = "background-color:white;padding:1% 1% 1% 1%",
                                                                                    column(12,
                                                                                           p("The timeline plot is not reactive. All studies in healthy individuals are included by default, and differences (sex, age, protocol) are blocked for in the statistical model."),
                                                                                           hr(),
                                                                                           plotOutput("human_timelinePlot_inactivity", height = "500px"),
                                                                                           hr(),
                                                                                           tableOutput("human_timelineTable_inactivity")
                                                                                    )
                                                                           )
                                                                  )
                                                       )
                                                       ,
                                                       
                                                       #= Tab Exercise Training =====================================================       
                                                       navbarMenu("Exercise Training",
                                                                  
                                                                  tabPanel("Aerobic",
                                                                           fluidRow(style = "background-color:white;padding:1% 1% 1% 1%",
                                                                                    column(10, style = "padding:0% 5% 2% 2%",
                                                                                           uiOutput("human_forestPlotDynamic_TA"),
                                                                                           hr(),
                                                                                           checkboxGroupInput("human_studies_TA", 
                                                                                                              tags$b("Include/exclude specific datasets"), 
                                                                                                              selected = human_input[['GEO_TA']],
                                                                                                              choices = human_input[['GEO_TA']],
                                                                                                              inline = TRUE)
                                                                                    ),
                                                                                    column(2, style = "background-color:#eaf1f7;margin-top:60px;padding:0 1% 1% 2%",
                                                                                           h3("Filter by", style = "color:#011b2a"),
                                                                                           sliderInput("human_training_duration", tags$b("Duration (weeks)"),
                                                                                                       min = 0, max = 1000, value = c(0,1000), step = 5, sep = ""),
                                                                                           tags$br(),
                                                                                           sliderInput("human_training_biopsy", tags$b("Biopsy time after exercise (h)"),
                                                                                                       min = 0, max = 96, value = c(0,96), step = 4, sep = "")
                                                                                    )
                                                                           )
                                                                  ),
                                                                  
                                                                  tabPanel("Combined",
                                                                           fluidRow(style = "background-color:white;padding:1% 1% 1% 1%",
                                                                                    column(10, style = "padding:0% 5% 2% 2%",
                                                                                           uiOutput("human_forestPlotDynamic_TC"),
                                                                                           hr(),
                                                                                           checkboxGroupInput("human_studies_TC", 
                                                                                                              tags$b("Include/exclude specific datasets"), 
                                                                                                              selected = human_input[['GEO_TC']],
                                                                                                              choices = human_input[['GEO_TC']],
                                                                                                              inline = TRUE)
                                                                                    ),
                                                                                    column(2, style = "background-color:#eaf1f7;margin-top:60px;padding:0 1% 1% 2%",
                                                                                           h3("Filter by", style = "color:#011b2a"),
                                                                                           sliderInput("human_training_duration", tags$b("Duration (weeks)"),
                                                                                                       min = 0, max = 1000, value = c(0,1000), step = 5, sep = ""),
                                                                                           tags$br(),
                                                                                           sliderInput("human_training_biopsy", tags$b("Biopsy time after exercise (h)"),
                                                                                                       min = 0, max = 96, value = c(0,96), step = 4, sep = "")
                                                                                    )
                                                                           )
                                                                  ),
                                                                  
                                                                  tabPanel("HIT",
                                                                           fluidRow(style = "background-color:white;padding:1% 1% 1% 1%",
                                                                                    column(10, style = "padding:0% 5% 2% 2%",
                                                                                           uiOutput("human_forestPlotDynamic_TH"),
                                                                                           hr(),
                                                                                           checkboxGroupInput("human_studies_TH", 
                                                                                                              tags$b("Include/exclude specific datasets"), 
                                                                                                              selected = human_input[['GEO_TH']],
                                                                                                              choices = human_input[['GEO_TH']],
                                                                                                              inline = TRUE)
                                                                                    ),
                                                                                    column(2, style = "background-color:#eaf1f7;margin-top:60px;padding:0 1% 1% 2%",
                                                                                           h3("Filter by", style = "color:#011b2a"),
                                                                                           sliderInput("human_training_duration", tags$b("Duration (weeks)"),
                                                                                                       min = 0, max = 1000, value = c(0,1000), step = 5, sep = ""),
                                                                                           tags$br(),
                                                                                           sliderInput("human_training_biopsy", tags$b("Biopsy time after exercise (h)"),
                                                                                                       min = 0, max = 96, value = c(0,96), step = 4, sep = "")
                                                                                    )
                                                                           )
                                                                  ),
                                                                  
                                                                  tabPanel("Resistance",
                                                                           fluidRow(style = "background-color:white;padding:1% 1% 1% 1%",
                                                                                    column(10, style = "padding:0% 5% 2% 2%",
                                                                                           uiOutput("human_forestPlotDynamic_TR"),
                                                                                           hr(),
                                                                                           checkboxGroupInput("human_studies_TR", 
                                                                                                              tags$b("Include/exclude specific datasets"), 
                                                                                                              selected = human_input[['GEO_TR']],
                                                                                                              choices = human_input[['GEO_TR']],
                                                                                                              inline = TRUE)
                                                                                    ),
                                                                                    column(2, style = "background-color:#eaf1f7;margin-top:60px;padding:0 1% 1% 2%",
                                                                                           h3("Filter by", style = "color:#011b2a"),
                                                                                           sliderInput("human_training_duration", tags$b("Duration (weeks)"),
                                                                                                       min = 0, max = 1000, value = c(0,1000), step = 5, sep = ""),
                                                                                           tags$br(),
                                                                                           sliderInput("human_training_biopsy", tags$b("Biopsy time after exercise (h)"),
                                                                                                       min = 0, max = 96, value = c(0,96), step = 4, sep = "")
                                                                                    )
                                                                           )
                                                                  )
                                                       )
                                                       ,
                                                       
                                                       
                                                       #= Tab correlations ===================================================== 
                                                       tabPanel("Co-regulation", 
                                                                value="panelAppInactivity",
                                                                fluidRow(style="background-color: white", 
                                                                         column(12, style="padding:2% 5% 2% 2%;margin:0",
                                                                                "This analyses identifies genes co-regulated during exercise.",
                                                                                tags$em("Click 'calculate' after changing selection criteria to update your analysis."),
                                                                                tags$br(), tags$br(),
                                                                                actionButton("updateCorrHuman", "Calculate", icon("sync")),
                                                                                checkboxGroupInput("human_corr_protocol",
                                                                                                   " ",
                                                                                                   choices=human_input[['protocol_choice']],
                                                                                                   selected=human_input[['protocol_choice']],
                                                                                                   inline = TRUE),
                                                                                
                                                                                tags$br(), tags$br()
                                                                         )
                                                                ),
                                                                fluidRow(style="background-color: white;padding:0% 5% 5% 2%",
                                                                         column(5, style="padding:0% 1% 2% 1%;margin:0",
                                                                                DT::dataTableOutput("human_corrTable", height="500px")
                                                                         ),
                                                                         
                                                                         column(7, style="padding:0% 1% 2% 1%;margin:0",
                                                                                plotOutput("human_corrPlot", height="500px"),
                                                                                textOutput("human_corrDescription"),
                                                                                uiOutput ("human_corrGeneCardsLink")
                                                                         )
                                                                )
                                                       )
                                           )
                                           
                                           
                                  ),
                                  
                                  # Row Select Population ###################################################################################
                                  tags$br(),
                                  fluidRow(style="position:center;background-color:#f4eae7;padding:0 0 1% 2%;margin:0 6% 5% 5%",
                                           h3("Customize your human population of interest", style="color:#c93f1e"),
                                           column(3, style="margin:0% 2% 0% 0%",
                                                  checkboxGroupInput("human_sex", tags$b("Sex"), 
                                                                     choices=human_input[['sex_choice']],
                                                                     selected=human_input[['sex_choice']]), #checkbox to select category
                                                  # checkboxInput('human_sex_allnone', 'All/None', value=T),
                                                  hr(),
                                                  checkboxGroupInput("human_age", tags$b("Age group"),
                                                                     choices=human_input[['age_choice']],
                                                                     selected=human_input[['age_choice']]), #checkbox to select category
                                                  # checkboxInput('human_age_allnone', 'All/None', value=T)
                                                  hr(),
                                                  checkboxGroupInput("human_fitness", tags$b("Fitness"),
                                                                     choices=human_input[['training_choice']],
                                                                     selected=human_input[['training_choice']]), #checkbox to select category
                                                  # checkboxInput('human_fitness_allnone', 'All/None', value=T),
                                           ), 
                                           column(3, style="margin:0% 2% 0% 0%",
                                                  checkboxGroupInput("human_weight", tags$b("BMI group"), 
                                                                     choices=human_input[['obesity_choice']],
                                                                     selected=human_input[['obesity_choice']]), #checkbox to select category
                                                  # checkboxInput('human_weight_allnone', 'All/None', value=T)
                                                  hr(),
                                                  checkboxGroupInput("human_muscle", tags$b("Muscle"), 
                                                                     choices=human_input[['muscle_choice']],
                                                                     selected=human_input[['muscle_choice']]),
                                                  checkboxInput('human_muscle_allnone', 'All/None', value=T)
                                           ), 
                                           column(4, style="margin:0% 2% 0% 0%",
                                                  checkboxGroupInput("human_disease", tags$b("Health status"), 
                                                                     choices=human_input[['disease_choice']],
                                                                     selected=human_input[['disease_choice']]), #checkbox to select category
                                                  checkboxInput('human_disease_allnone', 'All/None', value=T)
                                           )
                                  )
)
