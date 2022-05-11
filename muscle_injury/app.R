#-----------------------------------------------------------------------
#
# Transcriptomic profile of human and mouse macrophages
#
#----------------------------------------------------------------------
# Load data and libraries
library(shinycssloaders)
library(ggplot2)
theme <- theme(plot.title  = element_text(face="bold", color="black", size=13, angle=0),
               axis.text.x = element_text(face="bold", color="black", size=13, angle=0, hjust=0.5),
               axis.text.y = element_text(color="black", size=11, angle=0, vjust=0.3),
               axis.title  = element_text(face="bold", color="black", size=13, angle=0),
               legend.text = element_text(color="black", size=11, angle=0),
               legend.title = element_text(face="bold", color="black", size=12, angle=0),
               legend.position="right",
               legend.key.size = unit(22, "pt"))
library(stringr)
library(plyr)
library(DT)
library(dplyr)
library(DescTools)
library(ggpubr)
library(rstatix)
library(ggprism)


#--------------------------------------------------------------------------------------------------------
# Plot
#--------------------------------------------------------------------------------------------------------
data_all <- readRDS("data/MuscleInjury_data.Rds")

data_stats <- readRDS("data/MuscleInjury_stats.Rds")

samples <- readRDS("data/MuscleInjury_samples.Rds")

genelist <- rownames(data_all)

datasets <- readRDS("data/MuscleInjury_datasets.Rds")

#function to gene name
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Transcriptomic response of skeletal muscle to various insults"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", target="_blank", style="color:#D9DADB"), "/ last update 2021-09-14")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         selectizeInput("genesymbol", "Gene Symbol:", choices=NULL, multiple=F),
                         plotOutput("mousePlot") %>% withSpinner(color="#5b768e"),
                         tags$hr(),
                         h3("Statistics"),
                         "Results present mean and Sd for each protocol (A) and area under the curve (B)", tags$br(),
                         "Statistics: F-test performed separately for each protocol (compared to untreated animals). *FDR<0.05",
                         tags$hr(),
                         h3("Datasets included in the analysis"),
                         dataTableOutput("mouseDatasets"),
                         tags$hr(),
                         h3("All Genes included in the analysis"),
                         dataTableOutput("allStats")
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'genesymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected='Ppargc1a', 
                       options=NULL)
  
  
  output$mousePlot <- renderPlot({
    validate(need(input$genesymbol, " "))
    
    genetotest <- "Nsd2"
    genetotest <- input$genesymbol
    genedata <- t(data_all[genetotest,])
    testdata <- data.frame(samples[,2:3], genedata)
    colnames(testdata) <- c('Treatment', 'time', 'data')
    testdata$time <- as.numeric(gsub("T", "", testdata$time))
    testdata$group <- paste(testdata$Treatment, testdata$time, sep=".")
    
    #get stats from ggpubr
    stat.test <- compare_means(data=na.omit(testdata), data ~ time, group.by = "Treatment", ref.group = "0")
    
    # Add position for adjusted p-values
    stat.test <- stat.test %>% add_xy_position(data=testdata, x = "time", formula=data ~ time,  group="Treatment")
    stat.test$group <- paste(stat.test$Treatment, stat.test$group2, sep=".")
    
    #find n size
    nsize <- data.frame(nsize=table(testdata$group))
    testdata <- merge(testdata, nsize, by.x="group", by.y=1)
    
    #summarize data
    cdata <- ddply(testdata, c('Treatment', "time", 'group', 'nsize.Freq'), summarise,
                   mean = median(data, na.rm=TRUE),
                   sd   = sd(data, na.rm=TRUE))
    cdata$se <- cdata$sd / sqrt(cdata$nsize.Freq)
    cdata <- full_join(cdata, stat.test)
    cdata$y.position <- cdata$mean+cdata$se
    cdata$p.signif[cdata$p.signif=="ns"] <- NA
    
    plot_mean <- ggplot(cdata, aes(x=time, y=mean, color=Treatment, shape=Treatment, linetype=Treatment)) +
      theme_bw() + theme(legend.title=element_blank()) + theme +
      geom_hline(yintercept=0, linetype=3, color='grey30', size=0.5) +
      geom_point(size=3) + geom_line(aes(group=Treatment)) +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.4, size=0.5, linetype=1) +
      labs(title="",
           x="Time (hours)",
           y="-log10(p value), scaled") 
      #add_pvalue(cdata, bracket.size = NA,
      #           xmin = "time",
      #           xmax = "time",
      #           label = "p.signif",
      #           y.position = "y.position", label.size=6)
    
    #calculate AUC
    AUC_all <- numeric()
    for (i in unique(cdata$Treatment)){
      asd <- subset(cdata, Treatment %in% i)
      asd <- AUC(asd$time, asd$mean)
      AUC_all <- c(AUC_all, asd)
    }
    
    AUC_all <- data.frame(Treatment=unique(cdata$Treatment),
                          AUC=AUC_all)
    AUC_all$fdr <- as.numeric(data_stats[genetotest, c(2,4,6,8,10)])
    AUC_all$sign <- NA
    AUC_all$sign[AUC_all$fdr<0.05] <- "*"
    AUC_all$sign[AUC_all$fdr<0.01] <- "**"
    AUC_all$sign[AUC_all$fdr<0.001] <- "***"
    AUC_all$y.position <- ifelse(AUC_all$AUC > 0,
                                 AUC_all$AUC*1.1,
                                 AUC_all$AUC*1.5)
    
    #rename Treatments
    cdata$Treatment <- gsub("BLUNT", "Blunt", cdata$Treatment)
    cdata$Treatment <- gsub("ECCEN", "Eccentric\ncontraction", cdata$Treatment)
    cdata$Treatment <- gsub("FREEZ", "Freeze", cdata$Treatment)
    cdata$Treatment <- gsub("INCIS", "Incision", cdata$Treatment)
    cdata$Treatment <- gsub("EXERC", "Exercise", cdata$Treatment)
    AUC_all$Treatment <- gsub("BLUNT", "Blunt", AUC_all$Treatment)
    AUC_all$Treatment <- gsub("ECCEN", "Eccentric\ncontraction", AUC_all$Treatment)
    AUC_all$Treatment <- gsub("FREEZ", "Freeze", AUC_all$Treatment)
    AUC_all$Treatment <- gsub("INCIS", "Incision", AUC_all$Treatment)
    AUC_all$Treatment <- gsub("EXERC", "Exercise", AUC_all$Treatment)
    
    #plot AUC
    AUC_plot <- ggplot(AUC_all, aes(x=Treatment, y=AUC, fill=Treatment)) +
      theme_bw() + theme + theme(legend.position="none",
                                 axis.text.x = element_text(color="black", size=14, angle=45, hjust=1)) +
      geom_bar(stat = "identity", size=0.2, colour="black") +  
      geom_hline(yintercept=0) +
      labs(title=" ",
           x=element_blank(),
           y="Response to Treatment, AUC") +
      add_pvalue(AUC_all, bracket.size = NA,
                 xmin = "Treatment",
                 xmax = "Treatment",
                 label = "sign",
                 y.position = "y.position", label.size=8)
    AUC_plot
    
    cowplot::plot_grid(plot_mean, NA, AUC_plot,
                       ncol=3, rel_widths = c(2.5, 0.1, 1),
                       labels=c("A", "", "B"), label_size = 20)
  })
  
  output$mouseDatasets <- renderDataTable({
    datasets
  })
  
  output$allStats <- renderDataTable({
    datasets <- data.frame(SYMBOL=rownames(data_stats),
                           data_stats[,c(2,4,6,8,10)])
    
    datatable(datasets, options = list(lengthMenu = c(10, 50, 100), 
                                       pageLength = 10),
              rownames= FALSE) %>% 
      formatSignif(columns = c(2:6), digits = 3)
  })
  
  
}


# Run the app ----
shinyApp(ui = ui, server = server)