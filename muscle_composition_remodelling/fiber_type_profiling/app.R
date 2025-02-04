#-----------------------------------------------------------------------
#
# Profile of skeletal muscle fibers
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
library(tidyverse)
library(ggforce)


#--------------------------------------------------------------------------------------------------------
# Palmitate
#--------------------------------------------------------------------------------------------------------
data_proteome <- readRDS("data/data_proteome.Rds")
data_transcriptome <- readRDS("data/data_transcriptome.Rds")

# metadata
metadata_proteome <- readRDS("data/metadata_proteome.Rds")
metadata_transcriptome <- readRDS("data/metadata_transcriptome")

#--------------------------------------------------------------------------------------------------------
# List of genes
#--------------------------------------------------------------------------------------------------------
genelist <- unique(rownames(data_transcriptome))


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
                                h3("Proteome and transcriptome of individual skeletal muscle fibers"),
                                h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                           target="_blank", style="color:#D9DADB"), 
                                   "/ last update 2024-10-26")
                         ),
                         column(1,
                         )
                ),
                
                # main page
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         sidebarLayout(
                           sidebarPanel(width = 3,
                                        selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                                        actionButton("updatePlot", "Refresh plot", icon("refresh"))
                           ),
                           mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                                     plotOutput("MainPlot", height="800px") %>% withSpinner(color="#5b768e")
                                     )
                           )
                         ),
                         
                         tags$hr(),
                                  
                         # Table with references
                         fluidRow(style="color:black;background-color:white;padding:0% 2% 1% 2%;",
                                  h3("Datasets Included in the Analysis"),
                                  "Transcriptomics from manually isolated skeletal muscle fibers from", 
                                  a("Raue et al, 2012 (GSE28392)", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28392", 
                                    target="_blank", style="color:#5B768E"), 
                                  "and",
                                  a("Rubenstein et al, 2020 (GSE130977)", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130977", 
                                    target="_blank", style="color:#5B768E"),
                                  "for a total of 88 fibres.",
                                  tags$br(),
                                  "Proteomics from manually isolated skeletal muscle fibers from", 
                                  a("Murgia et al, 2017", href="https://doi.org/10.1016/j.celrep.2017.05.054", 
                                    target="_blank", style="color:#5B768E"), "and", 
                                  a("Murgia et al, 2022", href="https://doi.org/10.1093/pnasnexus/pgac086", 
                                    target="_blank", style="color:#5B768E"),
                                  "for a total of 387 fibres.",
                                  tags$br(),tags$br(),
                                  # "Raw data was downloaded, processed in unison, and corrected for batch effects to ensure consistency across samples. 
                                  # Fibers were identified based on the expression of the MYH isoforms: MYH1, MYH2, MYH4, and MYH7. For each fiber, the 
                                  # total intensity of these MYH isoforms was calculated, and the percentage contribution of each isoform to the total 
                                  # intensity was determined. Tertiles were then used to categorize fibers based on the relative abundance of each isoform. 
                                  # Fibers were annotated as 'high', 'mid', or 'low' for each isoform, depending on their placement within these tertiles. 
                                  # Fibers with high MYH7 expression were categorized as Type I, fibers with high MYH2 expression were categorized as Type IIA,
                                  # fibers with high MYH1 expression were categorized as Type IIX. Fibers that did not show a predominant expression of any single isoform, 
                                  # or showed intermediate levels across multiple isoforms, were classified as 'Mixed'."
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
                                       selected=c("LDHA", "LDHB", "MYH7", "MYH1", "MYH2"), 
                                       options=NULL)
                  
                  # AMPK plot
                  plotData <- eventReactive(input$updatePlot, {
                    validate(need(input$inputGeneSymbol, " "))
                    genes_of_interest <- c("LDHA", "LDHB", "MYH7", "MYH1", "MYH2")
                    genes_of_interest <- input$inputGeneSymbol
                    
                    # Extract and reshape data for plotting
                    mRNA <- data.frame(FiberType = metadata_transcriptome$FiberType,
                                       method = "transcriptome",
                                       t(data_transcriptome[genes_of_interest, ]))
                    protein <- data.frame(FiberType = metadata_proteome$Fiber,
                                          method = "proteome",
                                          t(data_proteome[genes_of_interest, ]))
                    
                    
                    plotdata <- full_join(mRNA, protein)
                    plotdata <- pivot_longer(plotdata, cols = all_of(genes_of_interest), names_to = "gene", values_to = "data")
                    plotdata$FiberType <- factor(plotdata$FiberType,
                                                 levels = c("Type I", "Mixed Type I/II", "Type IIA", "Mixed Type IIA/IIX", "Type IIX"))
                    
                    # exclude mixed fibers
                    plotdata <- plotdata[!grepl("Mixed", plotdata$FiberType),]
                    
                    # Compute ANOVA p-values for each gene and method
                    p_values <- plotdata %>%
                      group_by(gene, method) %>%
                      summarise(
                        p_value = ifelse(length(unique(FiberType)) > 1, 
                                         tryCatch(
                                           summary(aov(data ~ FiberType, data = cur_data()))[[1]][["Pr(>F)"]][1], 
                                           error = function(e) NA
                                         ), 
                                         NA),
                        .groups = "drop"
                      )
                    
                    # Function to format p-values
                    p_value_formatter <- function(p) {
                      sapply(p, function(x) {
                        if (!is.na(x) && x < 0.001) {
                          return("italic(p) < 0.001")
                        } else if (!is.na(x)) {
                          return(sprintf("italic(p) == %.3f", x))
                        } else {
                          return(NA)
                        }
                      })
                    }
                    
                    
                    # Add formatted p-values
                    p_values$formatted_p <- p_value_formatter(p_values$p_value)
                    
                    # Compute max data for each method
                    method_max_values <- plotdata %>%
                      group_by(method) %>%
                      summarise(y_position = max(data, na.rm = TRUE) * 1.05, .groups = "drop")
                    
                    # Merge to assign y_position based on method
                    p_values <- p_values %>%
                      left_join(method_max_values, by = "method")
                    
                    # Assign `FiberType` column to place p-values above **PAL** (not BSA)
                    p_values$FiberType <- "Type IIA"
                    
                    # Plot
                    ggplot(plotdata, aes(x = gene, y = data, fill = FiberType)) +
                      geom_boxplot(position = position_dodge(0.8), outlier.size = 0) +
                      geom_sina(size = 0.5, position = position_dodge(0.8), alpha = 0.5) +
                      facet_wrap(.~method, scales = "free", ncol = 1) +
                      theme_bw(base_size = 16) + 
                      labs(x = NULL, y = "Relative expression, log2") +
                      scale_y_continuous(expand = c(0, 4)) +
                      scale_fill_manual(values = c("Type I" = "#8B0000", "Type IIA" = "#F5DEB3", "Type IIX"= "#D3D3D3", "Mixed" = "#A0522D"))+
                      geom_text(data = p_values, aes(x = gene, y = y_position, label = formatted_p), 
                                parse = TRUE, size = 4, vjust = -1)  # Manually add p-values
                  })
                  
                  
                  output$MainPlot <- renderPlot({
                    plotData()
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