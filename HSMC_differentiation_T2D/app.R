#-----------------------------------------------------------------------
#
# Myotube differentiation
#
#----------------------------------------------------------------------
# Load libraries
library(shinycssloaders)
library(stringr)
library(plyr)
library(dplyr)
library(DT)
library(ggpubr)
library(rstatix)
library(ggplot2)
library(ggprism)
theme <- theme(plot.title  = element_text(face="bold", color="black", size=13, angle=0),
               axis.text.x = element_text(face="bold", color="black", size=13, angle=90, hjust=1, vjust=0.5),
               axis.text.y = element_text(color="black", size=11, angle=0, vjust=0.3),
               axis.title  = element_text(face="bold", color="black", size=13, angle=0),
               legend.text = element_text(color="black", size=13, angle=0),
               legend.title = element_blank(),
               legend.position="right",
               legend.key.size = unit(30, "pt"),
               strip.text = element_text(face="bold", color="black", size=14, angle=0))

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
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("In vitro differentiation of primary skeletal muscle cells"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), 
                            "/ last update", 
                            Sys.Date())
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 0% 8%;",
                         "Primary human myoblasts differentiated into myotubes from",
                         a("GSE55650", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55650", target="_blank", style="color:#5B768E"), 
                         "(Affymetrix HG-U133 plus2, n=5-6) and", 
                         a("GSE166502", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166502", target="_blank", style="color:#5B768E"), 
                         "(Illumina HumanHT-12 V4.0, n=13)"
                ),
                tags$hr(),
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                         actionButton("updatePlot", "Refresh plot", icon("refresh"))
                ),
                fluidRow(style="color:black;background-color:white;padding:2% 8% 1% 8%;",
                         plotOutput("genePlot", height="600px") %>% withSpinner(color="#5b768e", proxy.height="200px")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         tags$hr(),
                         h3("Statistics - differentiation"),
                         dataTableOutput("stats_diff"),
                         tags$hr(),
                         h3("Statistics - T2D"),
                         dataTableOutput("stats_t2d"),
                         tags$hr(),
                         h3("Statistics - interaction"),
                         dataTableOutput("stats_int")
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
    diff_list <- c("MKI67", "PAX7", "MYF5", "MSTN", "CKM", "AIM1", "ANOS1", "MYOD", "MYOG", "HMOX1")
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
                        stats_diffT2D[,c(1,10,11,7,8)]
                        #stats_T2Dblasts[,c(1,10,11,7,8)],
                        #stats_T2Dtubes[,c(1,10,11,7,8)]
                        )
    colnames(stat.limma)[1] <- "gene"
    
    #get stats from ggpubr
    stat.test <- compare_means(data=plotdata, data ~ group, 
                               group.by="gene")

    #merge with stats from limma
    stat.test <- full_join(stat.test, stat.limma)
    stat.test <- na.omit(stat.test)
    
    #replace stats by the ones from limma
    stat.test$p.adj.signif <- "ns"
    stat.test$p.adj.signif[stat.test$adj.P.Val<0.05] <- "*"
    stat.test$p.adj.signif[stat.test$adj.P.Val<0.01] <- "**"
    stat.test$p.adj.signif[stat.test$adj.P.Val<0.001] <- "***"
    
    # Add position for adjusted p-values
    stat.test <- stat.test %>% add_x_position(x = "group")

    #manually get y position for each gene
    max.y <- data.frame(aggregate(plotdata$data, by = list(plotdata$gene), max, na.rm=T))
    colnames(max.y) <- c("gene", "y.position")
    max.y$y.position <- max.y$y.position + 0.5
    stat.test <- full_join(stat.test, max.y)
    
    #stat.test <- stat.test %>% add_y_position(data=plotdata, formula=data ~ group,
    #                                          step.increase = 0.15, scales="free")
    
    #plot
    diffplot_mean <- ggboxplot(plotdata, x = "group", y = "data", facet.by = "gene", outlier.size=0.5, scales="free") +
      geom_hline(yintercept=0, linetype=3, color='grey30', size=0.5) +
      theme_bw(14) + theme + theme(legend.title = element_blank()) +
      labs(x=element_blank(),
           y="mRNA, log2(relative to control)") +
      stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = T, label.size = 8, bracket.size = 0.5) +
      geom_line(aes(group=subject, color=GEO),
                size=0.3) +
      scale_colour_manual(values=c("#56B4E9", "#D55E00")) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
      scale_x_discrete(labels=c("NGT.Myoblast" = "NGT\nMyoblast", 
                                "NGT.Myotube"  = "NGT\nMyotube",
                                "T2D.Myoblast" = "T2D\nMyoblast",
                                "T2D.Myotube"  = "T2D\nMyotube"))
    diffplot_mean
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
