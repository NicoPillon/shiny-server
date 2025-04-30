#-----------------------------------------------------------------------
#
# Diabetic Myotube differentiation
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


library(stringr)
library(plyr)
library(dplyr)
library(feather)
library(rstatix)

#load data
BlastDiffT2D_data <- readRDS('data/data_raw.Rds')

#stats, overall effect
BlastDiffT2D_stats_differentiation <- readRDS('data/data_stats_differentiation.Rds')
BlastDiffT2D_stats_interaction <- readRDS('data/data_stats_interaction.Rds')
BlastDiffT2D_stats_T2D <- readRDS('data/data_stats_T2D.Rds')

#stata pair-wise comparisons
BlastDiffT2D_stats_differentiation_NGT <- readRDS('data/data_stats_differentiation_NGT.Rds')
BlastDiffT2D_stats_differentiation_T2D <- readRDS('data/data_stats_differentiation_T2D.Rds')
BlastDiffT2D_stats_T2D_blasts <- readRDS('data/data_stats_T2D_blast.Rds')
BlastDiffT2D_stats_T2D_tubes <- readRDS('data/data_stats_T2D_tubes.Rds')

#list of genes
genelist <- rownames(BlastDiffT2D_data)


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
                                h3("Primary skeletal muscle cells from individuals with or without type 2 diabetes"),
                                h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                           target="_blank", style="color:#D9DADB"), 
                                   "/ last update 2021-09-30")
                         ),
                         column(1,
                         )
                          
                ),

                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         sidebarLayout(
                           sidebarPanel(width = 3,
                                        selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                                        actionButton("updatePlot", "Refresh plot", icon("refresh")),
                                        tags$br(),tags$br(),
                                        tags$b("Statistics"),
                                        em(h5("Linear model with empirical Bayes (eBayes) moderation in limma, adjusting for intra-study correlation. P-values are corrected with the false discovery rate."))
                           ),                           
                           mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                                     plotOutput("genePlot", height="600px") %>% withSpinner(color="#5b768e", proxy.height="200px")
                           )
                         )
                ),
                
                tags$hr(),
                
                # Table with datasets
                fluidRow(style="color:black;background-color:white;padding:0% 2% 1% 2%;",
                         h3("Datasets Included in the Analysis"),
                         "Primary human myoblasts differentiated into myotubes from",
                         a("GSE55650", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55650", target="_blank", style="color:#5B768E"), 
                         "(Affymetrix HG-U133 plus2, n=5-6) and", 
                         a("GSE166502", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166502", target="_blank", style="color:#5B768E"), 
                         "(Illumina HumanHT-12 V4.0, n=13)"
                ),
                
                tags$hr(),
                
                fluidRow(style="color:black;background-color:white;padding:0% 2% 1% 2%;",
                         h3("Statistics"),
                         h4("Overall effect of differentiation"),
                         dataTableOutput("stats_diff"),
                         h4("Overall effect of T2D"),
                         dataTableOutput("stats_t2d"),
                         h4("Interaction effect differentiation:diabetes"),
                         dataTableOutput("stats_int")
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
                       selected=c("MYF5", "CKM", "AIM1", "MYOD", "MYOG", "HMOX1"), 
                       options=NULL)
  
  plotData <- eventReactive(input$updatePlot, {
    validate(need(input$inputGeneSymbol, " "))
    diff_list <- c("MKI67", "PAX7", "MYF5", "MSTN", "CKM", "FBN2", "MYOD", "MYOG", "HMOX1")
    diff_list <- input$inputGeneSymbol
    
    diff_data <- t(BlastDiffT2D_data[diff_list,])
    colnames(diff_data) <- diff_list
    
    #normalize to control
    if (length(diff_list)>1){
    controldata <- colMeans(diff_data[grepl("Myoblast_NGT", rownames(diff_data)),], na.rm=T)
    } else controldata <- mean(diff_data[grepl("Myoblast_NGT", rownames(diff_data)),], na.rm=T)
      
    #calculate RQ by dividing everything by control means
    diff_data <- sweep(diff_data, 2, controldata, FUN="-")
    
    #add sample groups
    samples <- data.frame(str_split_fixed(rownames(diff_data), '_|\\.', 4))
    colnames(samples) <- c('GEO', 'differentiation', 'T2D', 'subject')
    samples$subject <- gsub('\\..*', '', samples$subject)
    
    diff_data <- data.frame(samples, diff_data)
    
    #plot data
    plotdata <- data.frame()
    for(i in 5:ncol(diff_data)){
      asd <- diff_data[,c(1:4, i)]
      asd$gene <- colnames(asd)[5]
      colnames(asd)[5] <- "data"
      plotdata <- rbind(plotdata, asd)  
    }
    
    plotdata$group <- paste(plotdata$T2D, plotdata$differentiation, sep=".")
    plotdata$group <- factor(plotdata$group, levels=c("NGT.Myoblast", "NGT.Myotube", "T2D.Myoblast", "T2D.Myotube"))
    
    #get stats from limma
    stats_diffNGT <- data.frame(BlastDiffT2D_stats_differentiation_NGT[diff_list,], group1="NGT.Myoblast", group2="NGT.Myotube")
    stats_diffT2D <- data.frame(BlastDiffT2D_stats_differentiation_T2D[diff_list,], group1="T2D.Myoblast", group2="T2D.Myotube")
    stats_T2Dblasts <- data.frame(BlastDiffT2D_stats_T2D_blasts[diff_list,], group1="NGT.Myoblast", group2="T2D.Myoblast")
    stats_T2Dtubes <- data.frame(BlastDiffT2D_stats_T2D_tubes[diff_list,], group1="NGT.Myotube", group2="T2D.Myotube")
    stat.limma <- rbind(stats_diffNGT[,c(1,10,11,7,8)],
                        stats_diffT2D[,c(1,10,11,7,8)],
                        #stats_T2Dblasts[,c(1,10,11,7,8)],
                        stats_T2Dtubes[,c(1,10,11,7,8)]
                        )
    colnames(stat.limma)[1] <- "gene"
    
    #get stats from ggpubr
    stat.test <- compare_means(data=plotdata, data ~ group, 
                               group.by="gene")

    #merge with stats from limma
    stat.test <- full_join(stat.test, stat.limma)
    stat.test <- stat.test[stat.test$adj.P.Val < 0.1,]
    stat.test <- na.omit(stat.test)
    
    # function to format p values
    p_value_formatter <- function(p) {
      sapply(p, function(x) {
        if (x < 0.001) {
          return("fdr < 0.001")
        } else {
          return(sprintf("fdr = %.3f", x))
        }
      })
    }
    
    #format p values
    stat.test$FDR_format <- p_value_formatter(stat.test$adj.P.Val)
    
    # Add position for adjusted p-values
    stat.test <- stat.test %>% add_x_position(x = "group")

    #manually get y position for each gene
    max.y <- data.frame(aggregate(plotdata$data, by = list(plotdata$gene), max, na.rm=T))
    colnames(max.y) <- c("gene", "y.position")
    max.y$y.position <- max.y$y.position + 0.5
    stat.test <- full_join(stat.test, max.y)

    #plot
    ggboxplot(plotdata, x = "group", y = "data", facet.by = "gene", outlier.size=0.5, scales="free") +
      geom_hline(yintercept=0, linetype=3, color='grey30', size=0.5) +
      theme_bw(16) + theme(legend.title = element_blank()) +
      labs(x=element_blank(),
           y="mRNA, log2(relative to control)") +
      stat_pvalue_manual(stat.test, label = "FDR_format", hide.ns = T, label.size = 4, bracket.size = 0.5) +
      geom_line(aes(group=subject, color=GEO),
                size=0.3) +
      scale_colour_manual(values=c("#56B4E9", "#D55E00")) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
      scale_x_discrete(labels=c("NGT.Myoblast" = "NGT\nMyoblast", 
                                "NGT.Myotube"  = "NGT\nMyotube",
                                "T2D.Myoblast" = "T2D\nMyoblast",
                                "T2D.Myotube"  = "T2D\nMyotube"))

  })

  output$genePlot <- renderPlot({
    plotData()
  })
  
  
  ##################################################################################################################
  #Statistic tables
  output$stats_diff <- renderDataTable(options=list(signif = 3),{
    datatable(BlastDiffT2D_stats_differentiation[2:nrow(BlastDiffT2D_data),], 
                  options = list(lengthMenu = c(5, 10, 50, 100), 
                                 pageLength = 5),
                  rownames= FALSE) %>% 
      formatSignif(columns = c(2:9), digits = 3)
  })
  
  output$stats_t2d <- renderDataTable({
    datatable(BlastDiffT2D_stats_T2D[2:nrow(BlastDiffT2D_data),], options = list(lengthMenu = c(5, 10, 50, 100), 
                                                    pageLength = 5),
                  rownames= FALSE) %>% 
      formatSignif(columns = c(2:9), digits = 3)
  })
  
  output$stats_int <- renderDataTable({
    datatable(BlastDiffT2D_stats_interaction[2:nrow(BlastDiffT2D_data),], options = list(lengthMenu = c(5, 10, 50, 100), 
                                                        pageLength = 5),
                  rownames= FALSE) %>% 
      formatSignif(columns = c(2:9), digits = 3)
  })
  
  
}


# Run the app ----
shinyApp(ui = ui, server = server)
