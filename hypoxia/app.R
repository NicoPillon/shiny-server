#-----------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle response to hypoxia
#
#----------------------------------------------------------------------
# Load data and libraries
library(shinycssloaders)
library(ggplot2)
theme <- theme(plot.title  = element_text(face="bold", color="black", size=13, angle=0),
               axis.text.x = element_text(face="bold", color="black", size=13, angle=45, hjust=1),
               axis.text.y = element_text(color="black", size=11, angle=0, vjust=0.3),
               axis.title  = element_text(face="bold", color="black", size=13, angle=0),
               legend.text = element_text(color="black", size=11, angle=0),
               legend.title = element_text(face="bold", color="black", size=12, angle=0),
               legend.position="right",
               legend.key.size = unit(18, "pt"))
library(stringr)
library(plyr)
library(DT)
library(dplyr)
library(DescTools)
library(ggpubr)
library(rstatix)


#--------------------------------------------------------------------------------------------------------
# C2C12
#--------------------------------------------------------------------------------------------------------
C2C12_data <- readRDS("data/C2C12_data.Rds")
C2C12_data <- C2C12_data[!grepl("Rik", rownames(C2C12_data)),]
C2C12_stats <- readRDS("data/C2C12_stats.Rds")

#--------------------------------------------------------------------------------------------------------
# HMDM
#--------------------------------------------------------------------------------------------------------
HMDM_data <- readRDS("data/HMDM_data.Rds")
HMDM_stats <- readRDS("data/HMDM_stats.Rds")

#--------------------------------------------------------------------------------------------------------
# Endothelial cells
#--------------------------------------------------------------------------------------------------------
EC_data <- readRDS("data/EC_data.Rds")
EC_stats <- readRDS("data/EC_stats.Rds")


#--------------------------------------------------------------------------------------------------------
# List of genes
#--------------------------------------------------------------------------------------------------------
genelist <- unique(c(rownames(C2C12_data), 
                     rownames(EC_data),
                     rownames(HMDM_stats)))


#function to gene name
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Transcriptomic response to hypoxia"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), "/ last update 2021-09-09")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         selectizeInput("genesymbol", "Gene Symbol:", choices=NULL, multiple=F),
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;text-align:center",
                         column(4,
                                h5(tags$b("Mouse C2C12 myotubes")),
                                plotOutput("C2C12Plot") %>% withSpinner(color="#5b768e"),
                         ),
                         column(4,
                                h5(tags$b("Human monocyte-derived macrophages")),
                                plotOutput("humanMDMPlot") %>% withSpinner(color="#5b768e")
                         ),
                         column(4,
                                h5(tags$b("Human endothelial cells")),
                                plotOutput("EndothelialPlot") %>% withSpinner(color="#5b768e")
                         )
                )
)



# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'genesymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected='Slc2a1', 
                       options=NULL)
  
  # C2C12 plot
  output$C2C12Plot <- renderPlot({
    validate(need(input$genesymbol, " "))
    genetotest <- firstup(input$genesymbol)
    validate(need(genetotest %in% rownames(C2C12_data)==T, "Gene undetectable"))
    
    #load data
    res_norm <- C2C12_data
    
    #find sample list
    Sample_norm <- data.frame(str_split_fixed(colnames(res_norm), "\\.|_", 4))[,1:3]
    colnames(Sample_norm) <- c("GEO", "treatment", "Repeat")
    
    #plot
    genedata <- t(res_norm[genetotest,])
    testdata <- data.frame(Sample_norm, genedata)
    colnames(testdata)[4] <- 'data'
    testdata$treatment <- gsub("[2-5]", "", testdata$treatment)
    
    testdata$time <- "2 days"
    testdata$time[testdata$GEO %in% "GSE132234"] <- "4 days"
    
    testdata$concentration <- 2
    testdata$concentration[testdata$GEO %in% "GSE132234"] <- 5
    testdata$legend <- paste(testdata$GEO, " (", testdata$concentration, "% ,", testdata$time, ")", sep="")
    
    testdata$treatment <- gsub("control", "Normoxia", testdata$treatment)
    testdata$treatment <- gsub("hypoxia", "Hypoxia", testdata$treatment)
    testdata$treatment <- factor(testdata$treatment, levels=c("Normoxia", "Hypoxia")) #for a box plot, x should be a factor
    
    #get stats from ggpubr
    stat.test <- t_test(testdata, data ~ treatment)
    
    #replace stats by the ones from limma
    stat.test$p <- C2C12_stats[genetotest,]$P.Value
    stat.test$p.adj <- C2C12_stats[genetotest,]$adj.P.Val
    stat.test$p.adj.signif <- "ns"
    stat.test$p.adj.signif[stat.test$p.adj<0.05] <- "*"
    stat.test$p.adj.signif[stat.test$p.adj<0.01] <- "**"
    stat.test$p.adj.signif[stat.test$p.adj<0.001] <- "***"
    
    # Add position for adjusted p-values
    stat.test <- stat.test %>% add_xy_position(x = "treatment")
    stat.test$y.position <- stat.test$y.position * 1.1
    
    #plot
    ggplot(testdata, aes(x=treatment, y=data)) +  
      geom_boxplot() + geom_point() +
      labs(x="",
           y=paste(genetotest)) +
      theme_bw() + theme + theme(legend.title=element_blank()) + 
      geom_line(aes(group=Repeat, color=legend),
                size=1) +
      stat_pvalue_manual(stat.test, 
                         label = "p.adj.signif", 
                         hide.ns = T, 
                         remove.bracket=F, 
                         size=7, 
                         bracket.size=0.5)
    
  })
  
  
  # Endothelial plot
  output$EndothelialPlot <- renderPlot({
    validate(need(input$genesymbol, " "))
    genetotest <- toupper(input$genesymbol)
    validate(need(genetotest %in% rownames(EC_data)==T, "Gene undetectable"))
    
    #load data
    res_norm <- EC_data
    
    #find sample list
    Sample_norm <- data.frame(str_split_fixed(colnames(res_norm), "\\.|_", 3))
    colnames(Sample_norm) <- c("GEO", "treatment", "Repeat")
    
    #collect data
    genedata <- t(res_norm[genetotest,])
    testdata <- data.frame(Sample_norm, genedata)
    colnames(testdata)[4] <- 'data'
    
    #add group time
    testdata$time <- gsub("hypoxia", "", testdata$treatment)
    testdata$time <- gsub("normoxia", "", testdata$time)
    
    #make only 2 groups: HYP and control
    testdata$treatment <- gsub("[0-9]", "", testdata$treatment)
    testdata$treatment <- gsub("hypoxia", "Hypoxia", testdata$treatment)
    testdata$treatment <- gsub("normoxia", "Normoxia", testdata$treatment)
    testdata$treatment <- factor(testdata$treatment, levels=c("Normoxia", "Hypoxia"))
    
    #label repeats
    testdata$Repeat <- paste(testdata$GEO, testdata$Repeat, sep=".")
    
    
    #add concentration of oxygen
    testdata$concentration <- 1
    testdata$concentration[testdata$GEO %in% c("GSE163827", "GSE107029")] <- 0.2
    testdata$concentration[testdata$GEO %in% c("GSE151610")] <- 0.5
    testdata$concentration[testdata$GEO %in% c("GSE86793")] <- 2
    
    testdata$legend <- paste(testdata$GEO, " (", testdata$concentration, "%, ", testdata$time, "h)", sep="")
    
    
    #get stats from ggpubr
    stat.test <- t_test(testdata, data ~ treatment)
    
    #replace stats by the ones from limma
    stat.test$p <- EC_stats[genetotest,]$P.Value
    stat.test$p.adj <- EC_stats[genetotest,]$adj.P.Val
    stat.test$p.adj.signif <- "ns"
    stat.test$p.adj.signif[stat.test$p.adj<0.05] <- "*"
    stat.test$p.adj.signif[stat.test$p.adj<0.01] <- "**"
    stat.test$p.adj.signif[stat.test$p.adj<0.001] <- "***"
    
    # Add position for adjusted p-values
    stat.test <- stat.test %>% add_xy_position(x = "treatment")
    stat.test$y.position <- stat.test$y.position * 1.1
    
    #plot
    ggplot(testdata, aes(x=treatment, y=data)) +  
      geom_point() +     
      geom_line(aes(group=Repeat, color=legend), size=1) +
      geom_boxplot(alpha=0) + 
      labs(x="",
           y=paste(genetotest)) +
      theme_bw() + theme + theme(legend.title=element_blank()) +
      stat_pvalue_manual(stat.test, 
                         label = "p.adj.signif", 
                         hide.ns = T, 
                         remove.bracket=F, 
                         size=5, 
                         bracket.size=0.5) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
    
  })

  # Human blood-derived macrophages
  output$humanMDMPlot <- renderPlot({
    validate(need(input$genesymbol, " "))
    genetotest <- toupper(input$genesymbol)
    validate(need(genetotest %in% rownames(HMDM_data)==T, "Gene undetectable"))
    
    #load data
    res_norm <- HMDM_data
    
    #find sample list
    Sample_norm <- data.frame(str_split_fixed(colnames(res_norm), "\\.|_", 4))
    colnames(Sample_norm) <- c("GEO", "treatment", "time", "repeat")
    
    #plot
    genedata <- t(res_norm[genetotest,])
    samples <- Sample_norm
    testdata <- data.frame(samples, genedata)
    colnames(testdata)[5] <- 'data'
    testdata$treatment <- factor(testdata$treatment, levels=c("Normoxia", "Hypoxia")) #for a box plot, x should be a factor
    
    testdata$pair <- gsub('_Hypoxia', '', colnames(res_norm))
    testdata$pair <- gsub('_Normoxia', '', testdata$pair)
    
    testdata$GEO <- gsub('_.*', '', colnames(res_norm))
    
    testdata$concentration <- 2
    testdata$concentration[testdata$GEO %in% "GSE16099"] <- 0.1
    testdata$concentration[testdata$GEO %in% "GSE15949"] <- 0.5
    testdata$concentration[testdata$GEO %in% "GSE4630"] <- 0
    
    testdata$legend <- paste(testdata$GEO, " (", testdata$concentration, "% ,", testdata$time, ")", sep="")
    
    #get stats from ggpubr
    stat.test <- t_test(testdata, data ~ treatment)
    
    #replace stats by the ones from limma
    stat.test$p <- c(HMDM_stats[genetotest,]$P.Value)
    stat.test$p.adj <- c(HMDM_stats[genetotest,]$adj.P.Val)
    stat.test$p.adj.signif <- "ns"
    stat.test$p.adj.signif[stat.test$p.adj<0.05] <- "*"
    stat.test$p.adj.signif[stat.test$p.adj<0.01] <- "**"
    stat.test$p.adj.signif[stat.test$p.adj<0.001] <- "***"
    
    # Add position for adjusted p-values
    stat.test <- stat.test %>% add_xy_position(x = "treatment")
    stat.test$y.position <- stat.test$y.position * 1.1
    
    #plot
    ggplot(testdata, aes(x=treatment, y=data)) +  
      geom_boxplot() + geom_point() +
      labs(x="",
           y=paste(genetotest)) +
      theme_bw() + theme + theme(legend.title=element_blank()) +
      geom_line(aes(group=pair, color=legend),
                size=1) +
      stat_pvalue_manual(stat.test, 
                         label = "p.adj.signif", 
                         hide.ns = T, 
                         remove.bracket=F, 
                         size=7, 
                         bracket.size=0.5)
    
  })
  
}


# Run the app ----
shinyApp(ui = ui, server = server)