#-----------------------------------------------------------------------
#
# Mice exercise
#
#----------------------------------------------------------------------
# Load libraries
library(shinycssloaders)
library(stringr)
library(DT)
library(plyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(feather)
library(ggforce)

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
                                style = "height:8vh; display:flex; justify-content:center; align-items:center;",
                                tags$a(href = "https://shiny.nicopillon.com", 
                                       icon("home", class = "fa-2x"), 
                                       style = "color:white; text-decoration:none;")  # Ensuring icon is white and no underline
                         ),
                         column(10,
                                h3("Gene expression in Tissues from mice after acute exercise"),
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
                                        tags$b("Statistics"),
                                        em(h5("ANOVA with GEO studies included as a covariate for each tissue. 
                                              P-values are not adjusted for multiple testing comparisons."))
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
                       selected=c("Nr4a3"), 
                       options=NULL)

  output$genePlot <- renderPlot({
    validate(need(input$inputGeneSymbol, " "))
    genename <- c("Ptgs2")
    genename <- input$inputGeneSymbol
    
    # collect data
    plotdata <- data.frame(metadata,
                           genedata = as.numeric(datamatrix[genename,]))
    
    #filter according to selected categories

    # label Tissues with n sizes
    plotdata$Tissue <- gsub("BAT", paste0("BAT\nn = ", nrow(plotdata[plotdata$Tissue == "BAT",])), plotdata$Tissue)
    plotdata$Tissue <- gsub("Liver", paste0("Liver\nn = ", nrow(plotdata[plotdata$Tissue == "Liver",])), plotdata$Tissue)
    plotdata$Tissue <- gsub("Skeletal Muscle", paste0("Skeletal Muscle\nn = ", nrow(plotdata[plotdata$Tissue == "Skeletal Muscle",])), plotdata$Tissue)
    plotdata$Tissue <- gsub("iWAT", paste0("iWAT\nn = ", nrow(plotdata[plotdata$Tissue == "iWAT",])), plotdata$Tissue)
    plotdata$Tissue <- gsub("Heart", paste0("Heart\nn = ", nrow(plotdata[plotdata$Tissue == "Heart",])), plotdata$Tissue)
    
    # express relative to each control
    plotdata_normalized <- plotdata %>%
      group_by(GEO, Tissue) %>%
      mutate(
        control_median = median(genedata[Time == "Control"], na.rm = TRUE), # Calculate median for "Control" in each group
        normalized_data = genedata - control_median # Normalize `data` to the control median
      ) %>%
      ungroup() %>%
      dplyr::select(-control_median) # Remove the temporary `control_median` column if no longer needed

    ##################################################################################
    # Statistics
    ##################################################################################
    
    # Split data based on number of GEO levels within each Tissue
    multi_geo <- plotdata_normalized %>%
      group_by(Tissue) %>%
      filter(n_distinct(GEO) > 1) %>%
      ungroup()
    
    single_geo <- plotdata_normalized %>%
      group_by(Tissue) %>%
      filter(n_distinct(GEO) == 1) %>%
      ungroup()
    
    # Run ANOVA with GEO as a covariate for tissues with multiple GEO levels
    results_multi_geo <- multi_geo %>%
      group_by(Tissue) %>%
      do({
        model <- aov(normalized_data ~ Time + GEO, data = .)
        anova_summary <- summary(model)
        p_value <- anova_summary[[1]]["Time", "Pr(>F)"]
        data.frame(Tissue = unique(.$Tissue), p_value = p_value)
      }) %>%
      ungroup()
    
    # Run ANOVA without GEO for tissues with a single GEO level
    results_single_geo <- single_geo %>%
      group_by(Tissue) %>%
      do({
        model <- aov(normalized_data ~ Time, data = .)
        anova_summary <- summary(model)
        p_value <- anova_summary[[1]]["Time", "Pr(>F)"]
        data.frame(Tissue = unique(.$Tissue), p_value = p_value)
      }) %>%
      ungroup()
    
    # Combine results
    results <- bind_rows(results_multi_geo, results_single_geo)
    
    # function to format p values
    p_value_formatter <- function(p) {
      sapply(p, function(x) {
        if (x < 0.001) {
          return("italic(p) < 0.001")
        } else {
          return(sprintf("italic(p) == %.3f", x))
        }
      })
    }
    results$p.label <- p_value_formatter(results$p_value)
    
    # plot
    ggplot(plotdata_normalized, aes(x=Time, y=normalized_data, fill = Time)) +
      facet_wrap_paginate(.~Tissue, ncol = 5, scales = "free_y") +
      geom_boxplot(position = position_dodge(0.8), outlier.size = 0) +
      geom_sina(size = 0.5, position = position_dodge(0.8)) +
      theme_bw(16) + 
      theme(legend.position = "none") +
      labs(x=element_blank(),
           y = paste(genename, "mRNA\nlog2(relative to control)")) +
      scale_fill_manual(values=c("gray90", "#B8DE29"))  +
      # stat_compare_means(aes(group = Time, 
      #                        label = after_stat(p_value_formatter(..p..))), 
      #                    size = 4, 
      #                    label.x = 1.5, hjust = 0.5,
      #                    vjust = 0, # Adjusted closer to help visibility
      #                    parse = TRUE) +
      geom_text(data = results, 
                aes(label = p.label, x = 1.5, y = Inf, vjust = -1), 
                parse = TRUE,
                vjust = 1.5, 
                size = 4, 
                inherit.aes = FALSE) +
      coord_cartesian(clip = "off") # Enables room for labels without affecting y-axis
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
          list(className = 'dt-center', targets = 3)  # Align column to the center
          # Add more lines for additional columns if needed
        )
      )
    )
  })
  
}


# Run the app ----
shinyApp(ui = ui, server = server)
