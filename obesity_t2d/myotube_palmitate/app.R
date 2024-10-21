#-----------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle response to palmitate
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


#--------------------------------------------------------------------------------------------------------
# Palmitate
#--------------------------------------------------------------------------------------------------------
PAL_data <- readRDS("data/data.Rds")

#get stats from limma
PAL_stats <- readRDS("data/stats.Rds")


#--------------------------------------------------------------------------------------------------------
# List of genes
#--------------------------------------------------------------------------------------------------------
genelist <- unique(rownames(PAL_data))


#--------------------------------------------------------------------------------------------------------
# Sample
#--------------------------------------------------------------------------------------------------------
#find sample list
Sample_norm <- data.frame(sample=colnames(PAL_data),
                          str_split_fixed(colnames(PAL_data), "_|\\.", 4))
colnames(Sample_norm) <- c("sample", "GEO", "treatment", "time", "repeats")

Sample_norm$concentration <- gsub("PAL", "", Sample_norm$treatment)
Sample_norm$concentration <- gsub("BSA", "0", Sample_norm$concentration)
Sample_norm$concentration <- as.numeric(Sample_norm$concentration)

Sample_norm$time <- as.numeric(gsub("H", "", Sample_norm$time))

Sample_norm$treatment <- gsub("[0-9].*", "", Sample_norm$treatment)

Sample_norm$model
Sample_norm$model[Sample_norm$GEO %in% "GSE6766"] <- "C2C12"
Sample_norm$model[Sample_norm$GEO %in% "GSE18589"] <- "HSMC"
Sample_norm$model[Sample_norm$GEO %in% "GSE38590"] <- "C2C12"
Sample_norm$model[Sample_norm$GEO %in% "GSE53116"] <- "LHCN-M2"
Sample_norm$model[Sample_norm$GEO %in% "GSE126101"] <- "HSMC"
Sample_norm$model[Sample_norm$GEO %in% "GSE205677"] <- "HSMC"

Sample_norm$repeats <- paste(Sample_norm$GEO,
                             Sample_norm$time,
                             Sample_norm$repeats)

#add description of stidues
Sample_norm$description <- paste0(
  Sample_norm$GEO, ", ",
  Sample_norm$model, ", ",
  Sample_norm$concentration, "uM, ",
  Sample_norm$time, "h"
)


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
                fluidRow(style="color:white;background-color:#5B768E;padding:0% 1% 1% 1%;text-align:center",
                         column(1, 
                                style = "height:8vh; display:flex; justify-content:center; align-items:center;",
                                tags$a(href = "https://shiny.nicopillon.com", 
                                       icon("home", class = "fa-2x"), 
                                       style = "color:white; text-decoration:none;")  # Ensuring icon is white and no underline
                         ),
                         column(10,
                                h3("Transcriptomic response of skeletal myotubes to palmitate"),
                                h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                           target="_blank", style="color:#D9DADB"), 
                                   "/ last update 2022-07-07")
                         ),
                         column(1,
                         )
                ),
                
                # main page
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         sidebarLayout(
                           sidebarPanel(width = 3,
                                        selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                                        checkboxGroupInput("concentration", 
                                                                     label = "Concentration (µM)", 
                                                                     selected = c(100,
                                                                                  200,
                                                                                  400,
                                                                                  500), 
                                                                     choices = c(100,
                                                                                 200,
                                                                                 400,
                                                                                 500)),
                                        checkboxGroupInput("time", 
                                                                     label = "Exposure (hours)", 
                                                                     selected = c(12,
                                                                                  16,
                                                                                  18,
                                                                                  20,
                                                                                  24,
                                                                                  30,
                                                                                  36,
                                                                                  42,
                                                                                  48,
                                                                                  54), 
                                                                     choices = c(12,
                                                                                 16,
                                                                                 18,
                                                                                 20,
                                                                                 24,
                                                                                 30,
                                                                                 36,
                                                                                 42,
                                                                                 48,
                                                                                 54)),
                                        checkboxGroupInput("model", 
                                                                     label = "Cell type", 
                                                                     selected = c("C2C12",
                                                                                  "HSMC",
                                                                                  "LHCN-M2"), 
                                                                     choices = c("Mouse C2C12 myotube" = "C2C12",
                                                                                 "Human primary myotube" = "HSMC",
                                                                                 "Human immortalized LHCN-M2" = "LHCN-M2")),
                                        actionButton("updatePlot", "Refresh plot", icon("refresh"))
                           ),
                           mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                                     plotOutput("PalPlot", height="600px") %>% withSpinner(color="#5b768e")
                                     )
                           )
                         ),
                         
                         
                         tags$hr(),
                         
                         # Statistics
                         fluidRow(style="color:black;background-color:white;padding:0% 2% 1% 2%;",
                                  h3("Statistics (calculated with limma on all samples)"),
                                  dataTableOutput("stats")
                         ),
                         
                         tags$hr(),
                                  
                         # Table with references
                         fluidRow(style="color:black;background-color:white;padding:0% 2% 1% 2%;",
                                  h3("Datasets Included in the Analysis"),
                                  "Myotubes exposed to BSA-conjugated palmitate from",
                                  a("GSE6766", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6766", 
                                    target="_blank", style="color:#5B768E"), 
                                  "(C2C12, n=3),",
                                  a("GSE18589", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18589", 
                                    target="_blank", style="color:#5B768E"), 
                                  "(HSMC, n=3),",
                                  a("GSE38590", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38590", 
                                    target="_blank", style="color:#5B768E"), 
                                  "(C2C12, n=1),",
                                  a("GSE53116", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53116", 
                                    target="_blank", style="color:#5B768E"), 
                                  "(LHCN-M2, n=2), ",
                                  a("GSE126101", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126101", 
                                    target="_blank", style="color:#5B768E"), 
                                  "(HSMC, n=4) and", 
                                  a("GSE205677", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205677", 
                                    target="_blank", style="color:#5B768E"), 
                                  "(HSMC, n=7).", 
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
                                       selected=c("PDK4", "ANGPTL4", "EPHA7", "NEFL", "IL6"), 
                                       options=NULL)
                  
                  # AMPK plot
                  plotData <- eventReactive(input$updatePlot, {
                    validate(need(input$inputGeneSymbol, " "))
                    genename <- c("GPSM2", "ANGPTL4", "HMOX1", "EPHA7", "IL6")
                    genename <- input$inputGeneSymbol
                    
                    #plot
                    plotdata <- data.frame()
                    for( i in genename) { 
                      data <- as.numeric(PAL_data[i,])              #collect data for gene name i
                      datay <- cbind.data.frame(Sample_norm[,2:8], data, gene=rep(i))   #create table with x="sample type", y="data", "gene name"
                      plotdata  <- rbind.data.frame(plotdata, datay)              #bind the new gene data at the bottom of the previous one
                    }
                    
                    plotdata$gene <- factor(plotdata$gene,
                                            levels=unique(str_sort(plotdata$gene, numeric = TRUE)))
                    
                    #filter according to selected categories
                    plotdata <- plotdata[plotdata$concentration %in% input$concentration &
                                           plotdata$time %in% input$time &
                                           plotdata$model %in% input$model 
                                         ,]
                    
                    # Apply the count of 'n' size for each unique gene
                    plotdata$gene <- sapply(plotdata$gene, function(g) {
                      gene_count <- nrow(na.omit(plotdata[plotdata$gene == g,]))
                      paste0(g, "\nn = ", gene_count)
                    })
                    
                    
                    # plot
                    ggplot(plotdata, aes(x = gene, y = data, fill = treatment)) +
                      geom_boxplot(position = position_dodge(0.8), outlier.size = 0) +
                      geom_sina(size = 0.5, position = position_dodge(0.8)) +
                      theme_bw(16) + 
                      labs(x=element_blank(),
                           y="mRNA expression, log2") +
                      scale_fill_manual(values=c("gray90", "#D55E00"))  +
                      stat_compare_means(aes(label = after_stat(p_value_formatter(..p..))), 
                                         size = 4, vjust = -2, parse = TRUE) +
                      scale_y_continuous(expand = c(0,1.5)) 
                  })
                  
                  
                  output$PalPlot <- renderPlot({
                    plotData()
                  })
                  
                  #Statistic tables
                  output$stats <- renderDataTable(options=list(signif = 3),{
                    datatable(PAL_stats, 
                              options = list(lengthMenu = c(10, 50, 100), 
                                             pageLength = 10),
                              rownames= TRUE) %>% 
                      formatSignif(columns = c(1:8), digits = 3)
                  })
                  
                }
                
                
                # Run the app ----
                shinyApp(ui = ui, server = server)