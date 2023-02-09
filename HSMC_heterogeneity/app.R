#-----------------------------------------------------------------------
#
# Muscle cell heterogeneity
#
#----------------------------------------------------------------------
# Load libraries
library(shinycssloaders)
library(stringr)
library(plyr)
library(DT)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(matrixStats)
library(ggprism)
theme <- theme(plot.title  = element_text(face="bold", color="black", size=13, angle=0),
               axis.text.x = element_text(face="bold", color="black", size=13, angle=90, hjust=0.5, vjust=0.5),
               axis.text.y = element_text(color="black", size=11, angle=0, vjust=0.3),
               axis.title  = element_text(face="bold", color="black", size=13, angle=0),
               legend.text = element_text(color="black", size=13, angle=0),
               legend.title = element_blank(),
               legend.position="right",
               legend.key.size = unit(30, "pt"),
               strip.text = element_text(face="bold", color="black", size=14, angle=0))

#load data
HSMC_data <- readRDS('data/data_corrected.Rds')
HSMC_sex_data <- readRDS('data/sex_data.Rds')
HSMC_sex_stats <- readRDS('data/sex_stats.Rds')

#genes
genelist <- rownames(HSMC_sex_stats)

#normalize to female
sample.data <- data.frame(t(data.frame(str_split(colnames(HSMC_sex_data), "_"))))
rownames(sample.data) <- c()
colnames(sample.data) <- c("donor", "experimenter", "condition", "platform")
sample.data$sex <- "M"
sample.data$sex[sample.data$donor %in% c("DK02", "DK04", "DK06", "DK08", "DK10")] <- "F"
res.F <- rowMedians(as.matrix(HSMC_sex_data[,sample.data$sex == "F"]))
HSMC_sex_data <- sweep(HSMC_sex_data, 1, res.F)


# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Heterogeneity and sex differences in primary human myotubes"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), 
                            "/ last update 2023-02-09")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 0% 8%;",
                         "Primary human myoblasts differentiated into myotubes from published and unpublished transcriptomic analyses
                         from the integrative physiology laboratory at Karolinska Institutet"
                ),
                tags$hr(),
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                         column(5, style="color:black;background-color:white;padding:0% 1% 1% 1%;text-align=center",
                                tags$b("Heterogeneity of primary human mytoubes"),
                                plotOutput("donorPlot", height="600px") %>% withSpinner(color="#5b768e")
                         ),
                         column(7, style="color:black;background-color:white;padding:0% 1% 1% 1%;text-align=center",
                                tags$b("Sex diffences in primary human mytoubes"),
                                plotOutput("sexPlot", height="600px") %>% withSpinner(color="#5b768e")
                         )
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         tags$hr(),
                         h3("Statistics - sex"),
                         dataTableOutput("stats_sex"),
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected=c("MYF5", "CKM", "EIF1AX", "DDX3Y"), 
                       options=NULL)

  output$donorPlot <- renderPlot({
    validate(need(input$inputGeneSymbol, " "))
    diff_list <- c("MKI67", "PAX7", "MYF5", "MSTN", "CKM", "DDX3Y", "MYOD", "EIF1AX", "HMOX1")
    diff_list <- input$inputGeneSymbol
    
    #cells
    sample.data <- data.frame(t(data.frame(str_split(colnames(HSMC_data), "_"))))
    rownames(sample.data) <- c()
    colnames(sample.data) <- c("donor", "experimenter", "condition", "platform")
    sample.data$sex <- "M"
    sample.data$sex[sample.data$donor %in% c("DK02", "DK04", "DK06", "DK08", "DK10")] <- "F"
    
    testdata   <- data.frame()                                #create an empty dataframe to collect data
    for( i in 1:length(diff_list)) { 
      y     <- as.numeric(HSMC_data[diff_list[i],])              #collect data for gene name i
      datay <- cbind.data.frame(sample.data, y, rep(diff_list[i]))   #create table with x="sample type", y="data", "gene name"
      colnames(datay) <- c("donor", "experimenter", "condition", "platform","sex", "data", "symbol") #rename column names for rbind later
      testdata  <- rbind.data.frame(testdata, datay)              #bind the new gene data at the bottom of the previous one
    }
    
    #plot
    plot_donor <- ggplot(testdata, aes(x=symbol, y=data)) + theme_bw(14) +
      geom_boxplot(outlier.colour = NA, fill="gray90") + 
      geom_jitter(aes(x=symbol, y=data, color=donor, shape=donor), width = 0.15) +
      labs(x=element_blank(),
           y="log2(relative mRNA)") +
      theme(legend.title=element_blank()) +
      scale_shape_manual(values=rep(c(15,1,17,5),10)) +
      guides(color=guide_legend(ncol=1), shape=guide_legend(ncol=1))
    plot_donor
    
  })
  
  output$sexPlot <- renderPlot({
    validate(need(input$inputGeneSymbol, " "))
    diff_list <- c("MKI67", "PAX7", "MYF5", "MSTN", "CKM", "DDX3Y", "MYOD", "EIF1AX", "HMOX1")
    diff_list <- input$inputGeneSymbol
    
    #cells
    sample.data <- data.frame(t(data.frame(str_split(colnames(HSMC_sex_data), "_"))))
    rownames(sample.data) <- c()
    colnames(sample.data) <- c("donor", "experimenter", "condition", "platform")
    sample.data$sex <- "M"
    sample.data$sex[sample.data$donor %in% c("DK02", "DK04", "DK06", "DK08", "DK10")] <- "F"

    testdata   <- data.frame()                                #create an empty dataframe to collect data
    for( i in 1:length(diff_list)) { 
      y     <- as.numeric(HSMC_sex_data[diff_list[i],])              #collect data for gene name i
      datay <- cbind.data.frame(sample.data, y, rep(diff_list[i]))   #create table with x="sample type", y="data", "gene name"
      colnames(datay) <- c("donor", "experimenter", "condition", "platform","sex", "data", "symbol") #rename column names for rbind later
      testdata  <- rbind.data.frame(testdata, datay)              #bind the new gene data at the bottom of the previous one
    }
    
    #get stats from limma
    stat.limma <- HSMC_sex_stats[diff_list,]
    stat.limma
    
    #get stats from ggpubr
    stat.test <- compare_means(data ~ sex, testdata, group.by="symbol")
    stat.test
    
    #check that data is in the same order
    stat.test$symbol == stat.limma$SYMBOL
    
    #replace stats by the ones from limma
    stat.test$p <- stat.limma$P.Value
    stat.test$p.adj <- stat.limma$adj.P.Val
    stat.test$p.adj.signif <- "ns"
    stat.test$p.adj.signif[stat.test$p.adj<0.05] <- "*"
    stat.test$p.adj.signif[stat.test$p.adj<0.01] <- "**"
    stat.test$p.adj.signif[stat.test$p.adj<0.001] <- "***"

    # Add position for adjusted p-values
    stat.test <- stat.test %>% add_xy_position(data=testdata, x = "symbol", formula=data ~ sex,  group="sex")
    stat.test$y.position <- max(testdata$data)*1.1
    
    #plot
    plot_sex <- ggboxplot(testdata, x = "symbol", y = "data", fill="sex", group="symbol", outlier.size=0.5) + 
      theme_bw(14) + theme(legend.title = element_blank()) +
      labs(x=element_blank(),
           y="log2(relative mRNA)") +
      scale_fill_manual(values=c("#CC79A7", "#0072B2")) +
      add_pvalue(stat.test, label.size=5,
                 xmin = "xmin", 
                 xmax = "xmax",
                 label = "p.adj.signif",
                 tip.length = 0.01)
    plot_sex
    return(plot_sex)
  })
  
  
  ##################################################################################################################
  #Statistic tables
  output$stats_sex <- renderDataTable(options=list(signif = 3),{
    datatable(HSMC_sex_stats[,c(1,5,2,7,8)], 
              options = list(lengthMenu = c(5, 10, 50, 100), 
                             pageLength = 5),
              rownames= FALSE) %>% 
      formatSignif(columns = c(2:5), digits = 3)
  })
  

}


# Run the app ----
shinyApp(ui = ui, server = server)
