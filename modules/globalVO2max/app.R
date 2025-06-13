#-----------------------------------------------------------------------
#
# Reference values for VO2 max
#
#----------------------------------------------------------------------
# Load data and libraries
library(shinycssloaders)
library(tidyverse)
library(ggplot2)
library(DT)
library(cowplot)
library(shinyWidgets)


VO2_data <- readRDS("data/data.Rds")

# lsit of studies
studies_included <- readRDS("data/studies_included.Rds")
countries <- sort(unique(VO2_data$Country))

# Define UI ----
ui <- fluidPage(
  # CSS for style
  tags$head(includeCSS("../../www/style.css")),
  
  # HTML code to display the countries in line
  tags$style(HTML("
  .multicol {
    -webkit-column-count: 2;  /* Chrome, Safari, Opera */
    -moz-column-count: 2;     /* Firefox */
    column-count: 2;
  }

  .bootstrap-select .dropdown-toggle .filter-option {
    white-space: normal !important;
  }
")),
  
  # main page
  navbarPage("GlobalVO2max",
             
             tabPanel(
               title = "Nomogram",
               
               sidebarPanel(width = 3,
                            selectInput("sex",
                                        label = "Your Sex", 
                                        selected = "Female", 
                                        choices = c("Female", "Male")),
                            numericInput("age",
                                         label = "Your Age (years)",
                                         value = 32),
                            numericInput("vo2max",
                                         label = "Your VO2max (mL/min/kg)",
                                         value = 33),
                            checkboxGroupInput("modality",
                                               label = "Modality", 
                                               selected = c("Cycle", "Treadmill"), 
                                               choices = c("Cycle", "Treadmill")),
                            # tags$b("Countries"),                             
                            # tags$br(),
                            # actionButton("select_all_countries", "Select All"),
                            # actionButton("deselect_all_countries", "Deselect All"),
                            # div(class = "multicol",
                            #     checkboxGroupInput("country",
                            #                        label = NULL,
                            #                        selected = countries,
                            #                        choices = countries)
                            # ),
                            shinyWidgets::pickerInput(
                              inputId = "country",
                              label = "Countries",
                              choices = countries,
                              selected = countries,
                              multiple = TRUE,
                              options = list(
                                `actions-box` = TRUE,
                                `live-search` = TRUE
                              )
                            ),
                            # tags$br(),
                            sliderInput("date", "Year",
                                        min = 1968, max = 2024, value = c(1968,2024), step = 1, sep = "")
               ),
               mainPanel(width = 9,
                         plotOutput("vo2Plot", height = "600px") %>% withSpinner(color="#5B768E", type = 8),
                         tags$br(),tags$br(),
                         textOutput("percentileText"),
                         tags$br(),
                         downloadButton("download_nomogram", "Download Nomogram (.pdf)"),
                         downloadButton("download_data", "Download Data (.csv)")
               )
             ),

             tabPanel(
               title = "Description",   
               
               h3("Citation"),
               p("Pillon NJ, Ortiz de Zevallos J, Zierath JR, LaMonte MK, Ainsworth BE. In Press, 2025. Stay tuned!"),
               tags$hr(),
               
               h3("Methods"),
               p("Our review employed a systematic methodology to collate and map existing evidence on oxygen uptake during maximal exercise testing by sex and age. This review followed PRISMA guidelines and was prospectively registered in PROSPERO (CRD42025641493). Studies were included if they directly measured oxygen uptake (VO2peak or VO2max) via incremental exercise testing with direct gas analysis in healthy adults (≥18 years). We excluded studies involving children, adolescents, athletes, or individuals with chronic diseases, as well as studies using indirect VO2peak estimates or mixing exercise modalities."),
               p("We systematically searched MEDLINE/PubMed, Embase, CINAHL, and Web of Science from inception to December 31, 2024. Reference lists of reviews and meta-analyses were also screened. Search terms included combinations of 'cardiopulmonary exercise testing', 'VO2max', 'VO2peak', and related synonyms, combined with terms for 'reference values' and 'age', while excluding animal studies, pediatric populations, athletes, and secondary sources (systematic reviews, meta-analyses). Titles and abstracts were screened independently by two reviewers in Covidence. Full-text reviews were conducted for eligible articles, with discrepancies resolved by consensus. One reviewer extracted study details (author, year, country), participant characteristics (age, sex, height, weight, BMI, sample size), testing modality (treadmill or cycle), and VO2peak statistics, verified by a second reviewer."),
               
               tags$hr(),

               h3("Datasets"),
               tags$p(
                 tags$b("Have we missed a relevant study? Do you want to contribute? Please "),
                 a("contact us", href = "mailto:nicolas.pillon@ki.se", target = "_blank"), "!"
               ),
               DT::dataTableOutput("studies_included")
             )
  ),
  
  
  # Code to send height to resizing iframe
  tags$head(
    tags$script(HTML("
    Shiny.addCustomMessageHandler('resizeFrame', function(message) {
      const height = document.documentElement.scrollHeight;
      parent.postMessage({ frameHeight: height }, '*');
    });

    window.addEventListener('resize', function() {
      const height = document.documentElement.scrollHeight;
      parent.postMessage({ frameHeight: height }, '*');
    });
  "))
  )
  
)


# Define server logic ----
server <- function(input, output, session) {
  
  # observeEvent(input$select_all_countries, {
  #   updateCheckboxGroupInput(session, "country", selected = countries)
  # })
  
  # observeEvent(input$deselect_all_countries, {
  #   updateCheckboxGroupInput(session, "country", selected = character(0))
  # })
  
  observe({
    available_countries <- VO2_data %>%
      filter(Sex == input$sex,
             Modality %in% input$modality,
             Publication_year >= input$date[1],
             Publication_year <= input$date[2]) %>%
      pull(Country) %>%
      unique()
    
    shinyWidgets::updatePickerInput(
      session,
      inputId = "country",
      choices = countries,
      selected = intersect(input$country, countries),
      choicesOpt = list(
        disabled = !(countries %in% available_countries)
      )
    )
  })
  
  # Reactive filtered dataset based on country, sex
  filtered_data <- reactive({
    # input values for testing
    input_sex = "Male"
    input_age = 40
    input_vo2max = 40
    input_modality = c("Cycle", "Treadmill")
    input_country = unique(VO2_data$Country)
    input_date_min = 2000
    input_date_max = 2022
    
    # capture inputs
    input_sex <- input$sex
    input_age <- input$age
    input_vo2max <- input$vo2max
    input_modality <- input$modality
    input_country <- input$country
    input_date_min <- input$date[1]
    input_date_max <- input$date[2]
    
    # validate numeric inputs
    validate(need(is.numeric(input_age) && !is.na(input_age), "Please provide a numeric value for age."))
    validate(need(is.numeric(input_vo2max) && !is.na(input_vo2max), "Please provide a numeric value for VO2max."))
    validate(need(length(input_modality) > 0, "Please select at least one modality (cycle/treadmill)."))
    validate(need(length(input_country) > 0, "Please select at least one country."))
    
    # filter data
    VO2_filtered <- VO2_data %>%
      filter(Sex == input_sex,
             Modality %in% input_modality,
             Country %in% input_country,
             Publication_year >= input_date_min,
             Publication_year <= input_date_max)
    VO2_filtered
  })
  
  
  output$percentileText <- renderPrint({
    df <- VO2_data
    df <- filtered_data()
    
    # check data is available
    validate(need(nrow(df) > 0, ""))
    
    
    # Fit a model per Percentile
    models <- df %>%
      group_by(Percentile) %>%
      do(model = lm(VO2max ~ Age_mean, data = .))
    
    # Predict VO2max at user's age for each percentile
    pred_df <- models %>%
      mutate(VO2_pred = predict(model, newdata = data.frame(Age_mean = input$age)))
    
    # Convert Percentile label to numeric
    percentile_to_number <- function(p) {
      if (p == "Median") return(50)
      as.numeric(gsub("[^0-9]", "", p))
    }
    pred_df$PercentileNum <- sapply(pred_df$Percentile, percentile_to_number)
    
    # Interpolate user's percentile
    approx_result <- approx(x = pred_df$VO2_pred,
                            y = pred_df$PercentileNum,
                            xout = input$vo2max,
                            rule = 2)
    
    user_percentile <- approx_result$y
    top_percent <- 100 - user_percentile
    
    # Access additional input variables
    user_sex <- tolower(input$sex)
    
    # Get unique modalities from filtered data
    modality_list <- unique(df$Modality)
    modality_str <- tolower(paste(modality_list, collapse = "/"))
    
    # Count number of unique studies in filtered data
    n_studies <- length(unique(paste(df$First_author, df$Publication_year)))  # replace Study_ID with the actual column name
    n_participants <- sum(df$N_size) / 9
    n_participants <- format(n_participants, big.mark = " ", scientific = FALSE)
    
    # Prepare message
    output_message <- sprintf(
      "Based on the criteria you selected, the nomogram was calculated from %d stud%s (%s participants). Your VO2max is at the %.0fth percentile of VO2max. This means that your fitness is in the top %.0f%% of the %s population of similar age, tested using %s ergometry across the selected countries.",
      n_studies,
      ifelse(n_studies == 1, "y", "ies"),
      n_participants,
      user_percentile,
      top_percent,
      user_sex,
      modality_str
    )
    
    if (user_percentile <= 40) {
      feedback <- "Try to move more — regular physical activity can help you improve your fitness, live longer, and stay healthier!"
    } else if (user_percentile <= 70) {
      feedback <- "Good job! Keep maintaining your fitness to support your health and well-being."
    } else {
      feedback <- "Fantastic! You’re in excellent shape — congratulations on your high fitness level!"
    }
    
    cat(paste(output_message, feedback))
    
  })
  
  
  

  output$vo2Plot <- renderPlot({
    df <- filtered_data()
    
    # check data is available
    validate(need(nrow(df) > 0, "No data available for the selected filters. Please broaden your selection (e.g., include more countries, modalities, or change sex)."))
    
    ggplot(df, aes(x = Age_mean, y = VO2max,
                   color = Percentile, linetype = Percentile)) +
      theme_bw(16) +
      theme(panel.grid.major = element_line(color = "gray80",
                                            linewidth = 0.5,
                                            linetype = 1),
            legend.title = element_blank()) +
      scale_x_continuous(breaks = seq(0, 100, 10),
                         minor_breaks = seq(0, 100, 1)) +
      scale_y_continuous(breaks = seq(0, 100, 10),
                         minor_breaks = seq(0, 100, 1)) +
      labs(x = "Age (Years)",
           y = expression(dot("V")["O"[2]] ~ "peak (mL/min/kg)")) +
      geom_smooth(se = FALSE, method = "lm", linewidth = 0.5, alpha = 0.9) +
      scale_colour_manual(values = c("#097910", "#1d7c0f", "#38800d", "#53840c",
                                     "black",
                                     "#8c8d08", "#a99106", "#c99604", "#e69a02")) +
      scale_linetype_manual(values = c(5,4,3,2,1,2,3,4,5)) +
      annotate("label", 
               x = input$age, 
               y = input$vo2max * 0.94, 
               alpha = 0.5,
               color = "red",
               label = "Your VO2max") +
      annotate("pointrange",
               x = input$age,
               y = input$vo2max,
               ymin = input$vo2max, ymax = input$vo2max,
               color = "red",
               size = 0.5) +
      coord_fixed(1.5)
  })
  
  
  output$studies_included <- DT::renderDataTable(escape = FALSE, 
                                                 rownames = FALSE, 
                                                 options=list(paging = FALSE,
                                                              dom = 'ft'), 
                                                 {
                                                   studies_included
                                                 })
  
  output$download_nomogram <- downloadHandler(
    filename = function() {
      paste("nomogram_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      df <- filtered_data()
      

      # build plot
      x_range <- diff(range(df$Age_mean, na.rm = TRUE))
      y_range <- diff(range(df$VO2max, na.rm = TRUE))
      ratio <- y_range / x_range
      
      p <- ggplot(df, aes(x = Age_mean, y = VO2max,
                          color = Percentile, linetype = Percentile)) +
        theme_bw(10) +
        theme(panel.grid.major = element_line(color = "gray80",
                                              linewidth = 0.5,
                                              linetype = 1),
              legend.title = element_blank()) +
        scale_x_continuous(breaks = seq(0, 100, 10),
                           minor_breaks = seq(0, 100, 1)) +
        scale_y_continuous(breaks = seq(0, 100, 10),
                           minor_breaks = seq(0, 100, 1)) +
        labs(x = "Age (Years)",
             y = expression(dot("V")["O"[2]] ~ "peak (mL/min/kg)")) +
        geom_smooth(se = FALSE, method = "lm", linewidth = 0.5, alpha = 0.9) +
        scale_colour_manual(values = c("#097910", "#1d7c0f", "#38800d", "#53840c",
                                       "black",
                                       "#8c8d08", "#a99106", "#c99604", "#e69a02")) +
        scale_linetype_manual(values = c(5,4,3,2,1,2,3,4,5)) +
        coord_fixed(ratio = 1.5)
      
      
      # build caption
      selected_sex <- unique(df$Sex)
      selected_modality <- paste(unique(df$Modality), collapse = "/")
      selected_countries <- paste(unique(df$Country), collapse = ", ")
      selected_years <- paste(min(df$Publication_year), "-", max(df$Publication_year))
      
      
      caption_text <- stringr::str_wrap(
        paste0(
          "\n",
        "Source: Pillon NJ, Ortiz de Zevallos J, Zierath JR, LaMonte MK, Ainsworth BE. Submitted manuscript, coming soon.\n\n",
        "The nomogram was created with the following parameters:\n",
        "Sex: ", selected_sex, ", \n",
        "Modality: ", selected_modality, ", \n",
        "Years: ", selected_years, ", \n",
        "Countries: ", selected_countries, "."
        )
        ,
        width = 135
        )
      
      caption_plot <- ggdraw() + 
        draw_label(caption_text, x = 0, y = 1, hjust = 0, vjust = 1, size = 9, fontface = "italic")
      
      final_plot <- plot_grid(p, caption_plot, ncol = 1, rel_heights = c(1, 0.2))
      
      final_plot_with_margin <- ggdraw(final_plot) + 
        theme(plot.margin = margin(t = 20, r = 20, b = 200, l = 20))
      
      # Save as PDF in A4 size
      ggsave(file, plot = final_plot_with_margin, device = "pdf", width = 8.27, height = 11.69, units = "in")
    }
  )
  
  output$download_data <- downloadHandler(
    filename = function() {
      paste("VO2_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      df <- filtered_data()
      
      # pivot wider: Percentile values become columns
      df_wide <- tidyr::pivot_wider(
        df,
        names_from = Percentile,
        values_from = VO2max
      )
      
      # Round all numeric columns to 1 decimal
      df_wide <- df_wide %>%
        mutate(
          Height_cm = round(Height_cm, 0),  # round height to nearest integer
          across(where(is.numeric) & !matches("Height_cm"), ~ round(.x, 1))  # round all other numeric vars to 1 decimal
        )
      
      write.csv(df_wide, file, row.names = FALSE)
    }
  )
  
}

# Run the app ----
shinyApp(ui = ui, server = server)