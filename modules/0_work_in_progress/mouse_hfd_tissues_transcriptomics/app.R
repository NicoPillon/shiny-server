#-----------------------------------------------------------------------
#
# Mice on HFD
#
#----------------------------------------------------------------------
# Load libraries
library(feather)
library(shinycssloaders)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(cowplot)
library(DT)

#load metadata
genelist <- readRDS("data/genelist.Rds")
metadata <- readRDS("data/metadata.Rds")
references <- readRDS("data/references.Rds")

# load matrix
datamatrix_1 <- read_feather('data/datamatrix_1.feather')
datamatrix_2 <- read_feather('data/datamatrix_2.feather')
datamatrix_3 <- read_feather('data/datamatrix_3.feather')
datamatrix <- data.frame(rbind(datamatrix_1,
                               datamatrix_2,
                               datamatrix_3))
rownames(datamatrix) <- genelist

# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                
                # Google analytics
                tags$head(includeScript("google-analytics.html")),
                
                # Custom CSS to change checkbox tick color
                tags$style(HTML("
                  input[type='checkbox'] {
                    accent-color: #c93f1e; /* Change the checkbox tick color here */
                  }
                ")),
                
                # title ribbon
                fluidRow(style="color:white;background-color:#5B768E;padding:0% 1% 1% 1%;text-align:center; display:flex; justify-content:center; align-items:center;",
                         column(1, 
                                tags$a(href = "https://shiny.nicopillon.com", 
                                       icon("home", class = "fa-2x"), 
                                       style = "color:white; text-decoration:none;")  # Ensuring icon is white and no underline
                         ),
                         column(10,
                                h3("Gene expression in tissues from mice fed a high fat diet"),
                                h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                           target="_blank", style="color:#D9DADB"), 
                                   "/ last update 2024-10-16")
                                ),
                         column(1,
                         )
                         ),
                
                # main page
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         sidebarLayout(
                           sidebarPanel(width = 3,
                                        selectizeInput("inputGeneSymbol", "Gene Symbol:", choices=NULL, multiple=F, width=1000),
                                        # sliderInput("age_start", tags$b("Age at start (weeks)"),
                                        #             min = 3, max = 14, value = c(3,14), step = 1, sep = ""),
                                        sliderInput("diet_duration", tags$b("Diet duration (weeks)"),
                                                    min = 0.5, max = 25, value = c(0.5,25), step = 1, sep = ""),
                                        sliderInput("diet_composition", tags$b("Diet fat content (%)"),
                                                    min = 15, max = 60, value = c(15,60), step = 1, sep = ""),
                                        tags$b("Statistics"),
                                        em(h5("Wilcoxon ranked signed test comparing HFD to control. The p values are not adjusted for multiple testing comparisons.")),
                                        em(h5("In this analysis of more than 700 animals, less than 20 are female, making it impossible to analyse sex differences."))
                                        
                           ),
                           mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                                     plotOutput("genePlot", height="500px") %>% withSpinner(color="#5b768e")
                                     )
                         )
                ),
                
                tags$hr(),
                
                # Table with references
                fluidRow(style="color:black;background-color:white;padding:0% 2% 1% 2%;",
                         h3("Datasets Included in the Analysis"),
                         dataTableOutput("references")
                ),
                
                
                # Author section at the bottom
                fluidRow(style="color:white;background-color:#5B768E;padding:2% 1% 2% 1%;display: flex; align-items: top; ",
                         # column(2, align="right", 
                         #        tags$img(src = "https://ki.se/profile-image/nicpil", height = "120px", width = "120px")  # Insert image here
                         # ),
                         column(4, align="left", 
                                tags$b("About the author:"), tags$br(),
                                "Nicolas J. Pillon, PhD", tags$br(),
                                "Associate Professor, Karolinska Institutet", tags$br(),
                                icon("globe"), a("/inflammation-and-metabolism", href="https://ki.se/en/research/research-areas-centres-and-networks/research-groups/inflammation-and-metabolism-nicolas-pillons-research-group",
                                                 target="_blank", style="color:white"), tags$br(),
                                icon("linkedin"), a("/nicopillon", href="https://www.linkedin.com/in/nicopillon/",
                                                    target="_blank", style="color:white"), tags$br(),
                                tags$br(),
                                "Feel free to write to me with feedback or questions:", tags$br(),
                                icon("envelope"), a("nicolas.pillon@ki.se", href="mailto:nicolas.pillon@ki.se",
                                                    target="_blank", style="color:white"), tags$br(),
                                
                         ),
                         column(4, align="center",
                                #tags$b("© 2024 Nicolas Pillon"), tags$br(),
                                
                         ),
                         column(4, align="right",
                                tags$b("Disclaimer:"), tags$br(),
                                em("The authors disclaim any responsibility for the use or interpretation of the data 
                                   presented in this application. Users are solely responsible for ensuring the appropriate 
                                   use of any data they choose to re-use."), tags$br(),
                                tags$br(),
                                tags$b("© 2024 Nicolas Pillon"),
                         ),
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected=c("Lep"), 
                       options=NULL)

  output$genePlot <- renderPlot({
    validate(need(input$inputGeneSymbol, " "))
    genename <- c("Itgax")
    genename <- input$inputGeneSymbol
    
    # collect data
    plotdata <- data.frame(metadata,
                           genedata = as.numeric(datamatrix[genename,]))
    
    #filter according to selected categories
    plotdata <- dplyr::filter(plotdata,
                              diet.duration >= input$diet_duration[1] & diet.duration <= input$diet_duration[2],
                              #age.at.start >= input$age_start[1] & age.at.start <= input$age_start[2],
                              diet.fat.percent >= input$diet_composition[1] & diet.fat.percent <= input$diet_composition[2])

    # label tissues with n sizes
    plotdata$tissue <- gsub("BAT", paste0("BAT\nn = ", nrow(plotdata[plotdata$tissue == "BAT",])), plotdata$tissue)
    plotdata$tissue <- gsub("Liver", paste0("Liver\nn = ", nrow(plotdata[plotdata$tissue == "Liver",])), plotdata$tissue)
    plotdata$tissue <- gsub("Muscle", paste0("Muscle\nn = ", nrow(plotdata[plotdata$tissue == "Muscle",])), plotdata$tissue)
    plotdata$tissue <- gsub("subcutWAT", paste0("subcutWAT\nn = ", nrow(plotdata[plotdata$tissue == "subcutWAT",])), plotdata$tissue)
    plotdata$tissue <- gsub("visceralWAT", paste0("visceralWAT\nn = ", nrow(plotdata[plotdata$tissue == "visceralWAT",])), plotdata$tissue)
    
    # express relative to each control
    plotdata_normalized <- plotdata %>%
      group_by(GEO, tissue) %>%
      mutate(
        control_median = median(genedata[diet == "chow"], na.rm = TRUE), # Calculate median for "Control" in each group
        normalized_data = genedata - control_median # Normalize `data` to the control median
      ) %>%
      ungroup() %>%
      dplyr::select(-control_median) # Remove the temporary `control_median` column if no longer needed
    
    # Compute p-values **for each tissue**
    p_values <- plotdata_normalized %>%
      group_by(tissue) %>%
      summarise(
        p = wilcox.test(normalized_data ~ diet, data = cur_data(), exact = FALSE)$p.value
      ) %>%
      ungroup()
    
    # Define function to format p-values
    p_value_formatter <- function(p) {
      sapply(p, function(x) {
        if (x < 0.001) {
          return("italic(p) < 0.001")
        } else {
          return(sprintf("italic(p) == %.3f", x))
        }
      })
    }
    
    # Create formatted p-values
    p_values$formatted_p <- p_value_formatter(p_values$p)
    
    # Define y-position for annotation (above max value for each tissue)
    y_pos <- plotdata_normalized %>%
      group_by(tissue) %>%
      summarise(y_position = max(normalized_data, na.rm = TRUE) * 1.1)
    
    # Merge y_position into p_values dataframe
    p_values <- left_join(p_values, y_pos, by = "tissue")
    
    # Assign `diet` column to place p-values at `HFD`
    p_values$diet <- "HFD"
    
    # plot
    ggplot(plotdata_normalized, aes(x=diet, y=normalized_data, fill = diet)) +
      facet_wrap_paginate(.~tissue, ncol = 5, scales = "free_y") +
      geom_boxplot(position = position_dodge(0.8), outlier.size = 0) +
      geom_sina(size = 0.5, position = position_dodge(0.8)) +
      theme_bw(16) + 
      theme(legend.position = "none") +
      labs(x=element_blank(),
           y = paste(genename, "mRNA\nlog2(relative to control)")) +
      scale_fill_manual(values=c("gray90", "#FDE725"))  +
      scale_y_continuous(expand = expansion(mult = c(0, .1)))  +
      geom_text(data = p_values, aes(x = 1.5, y = y_position, label = formatted_p), 
                parse = TRUE, size = 4, vjust = -1)  # Manually add p-values
  })
  
  
  ##################################################################################################################
  #Dataset tables
  output$references <- renderDataTable(options=list(signif = 3),{
    DT::datatable(
      references, 
      escape = FALSE, 
      rownames = FALSE,
      options = list(
        columnDefs = list(
          list(className = 'dt-center', targets = 2),  # Align column to the center
          list(className = 'dt-center', targets = 3),  # Align column to the center
          list(className = 'dt-center', targets = 4)  # Align column to the right
          # Add more lines for additional columns if needed
        )
      )
    )
  })
  
}


# Run the app ----
shinyApp(ui = ui, server = server)
