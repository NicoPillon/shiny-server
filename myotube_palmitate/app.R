#-----------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle response to palmitate
#
#----------------------------------------------------------------------
# Load data and libraries
library(shinycssloaders)
library(ggplot2)
theme <- theme(plot.title  = element_text(face="bold", color="black", size=13, angle=0),
               axis.text.x = element_text(face="bold", color="black", size=13, angle=45, hjust=1),
               axis.text.y = element_text(color="black", size=11, angle=0, vjust=0.3),
               axis.title  = element_text(face="bold", color="black", size=13, angle=0),
               legend.text = element_text(color="black", size=12, angle=0),
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
library(ggprism)


#--------------------------------------------------------------------------------------------------------
# AMPK
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
                          str_split_fixed(colnames(PAL_data), "_|\\.", 3))
colnames(Sample_norm) <- c("sample", "GEO", "group", "repeats")

Sample_norm$group <- gsub("[0-9].*", "", Sample_norm$group)

#batch pairs
Sample_norm$repeats <- gsub('BSA', '', Sample_norm$sample)
Sample_norm$repeats <- gsub('PAL200', '', Sample_norm$repeats)
Sample_norm$repeats <- gsub('PAL500', '', Sample_norm$repeats)



# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Transcriptomic response of skeletal myotubes to palmitate"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), "/ last update 2021-10-04")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 0% 8%;",
                         "Myotubes exposed to BSA-conjugated palmitate from",
                         a("GSE6766", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6766", target="_blank", style="color:#5B768E"), 
                         "(C2C12, n=3),",
                         a("GSE18589", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18589", target="_blank", style="color:#5B768E"), 
                         "(HSMC, n=3),",
                         a("GSE38590", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38590", target="_blank", style="color:#5B768E"), 
                         "(C2C12, n=1),",
                         a("GSE53116", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53116", target="_blank", style="color:#5B768E"), 
                         "(LHCN-M2, n=2) and",
                         a("GSE126101", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126101", target="_blank", style="color:#5B768E"), 
                         "(HSMC, n=4).", 


                ),
                tags$hr(),
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                         actionButton("updatePlot", "Refresh plot", icon("refresh"))
                ),
                fluidRow(style="color:black;background-color:white;padding:2% 8% 1% 8%;",
                         plotOutput("PalPlot", height="400px") %>% withSpinner(color="#5b768e")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         tags$hr(),
                         h3("Statistics"),
                         dataTableOutput("stats")
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected=c("GPSM2", "ANGPTL4", "HMOX1", "CCL2", "IL6"), 
                       options=NULL)
  
  # AMPK plot
  plotData <- eventReactive(input$updatePlot, {
    validate(need(input$inputGeneSymbol, " "))
    genename <- c("GPSM2", "ANGPTL4", "HMOX1", "CCL2", "IL6")
    genename <- input$inputGeneSymbol
    
    #plot
    testdata <- data.frame()
    for( i in genename) { 
      y     <- as.numeric(PAL_data[i,])              #collect data for gene name i
      datay <- cbind.data.frame(Sample_norm[,2:4], y, rep(i))   #create table with x="sample type", y="data", "gene name"
      colnames(datay) <- c("GEO", "treatment","Repeat", "data", "gene")                #rename column names to make it possible to rbind later
      testdata  <- rbind.data.frame(testdata, datay)              #bind the new gene data at the bottom of the previous one
    }
    
    testdata$concentration <- ""
    testdata$concentration[testdata$GEO %in% "GSE6766"] <- "C2C12, 500uM, 16h"
    testdata$concentration[testdata$GEO %in% "GSE18589"] <- "HSMC, 100uM, 24h"
    testdata$concentration[testdata$GEO %in% "GSE38590"] <- "C2C12, 200uM, 20h"
    testdata$concentration[testdata$GEO %in% "GSE53116"] <- "LHCN-M2, 500uM, 16h"
    testdata$concentration[testdata$GEO %in% "GSE126101"] <- "HSMC, 500uM, 48h"
    testdata$legend <- paste(testdata$GEO, " (", testdata$concentration, ")", sep="")
    
    testdata$gene <- factor(testdata$gene,
                            levels=unique(str_sort(testdata$gene, numeric = TRUE)))
    
    #get stats from ggpubr
    stat.test <- compare_means(data=na.omit(testdata), data ~ treatment, group.by = "gene")
    
    #replace stats by the ones from limma
    stat.test$p <- na.omit(PAL_stats[genename,])$P.Value
    stat.test$p.adj <- na.omit(PAL_stats[genename,])$adj.P.Val
    stat.test$p.adj.signif <- "ns"
    stat.test$p.adj.signif[stat.test$p.adj<0.05] <- "*"
    stat.test$p.adj.signif[stat.test$p.adj<0.01] <- "**"
    stat.test$p.adj.signif[stat.test$p.adj<0.001] <- "***"
    
    #plot 2 - only treated condition relative to each basal
    foldchange <- testdata[testdata$treatment %in% "PAL",]
    foldchange$data <- testdata[testdata$treatment %in% "PAL",]$data - testdata[testdata$treatment %in% "BSA",]$data 
    
    #manually get y position for each gene
    max.y <- data.frame(aggregate(foldchange$data, by = list(foldchange$gene), max, na.rm=T))
    colnames(max.y) <- c("gene", "y.position")
    max.y$y.position <- max.y$y.position + 0.5
    stat.test <- full_join(stat.test, max.y)
    
    #remove ns
    stat.test$p.adj.signif[stat.test$p.adj.signif=="ns"] <- ""

    #bar plot
    plot_C2C12 <- ggbarplot(foldchange, x="gene", y="data", add="mean_se", fill="gray80", size=0.5) +  
      geom_hline(yintercept = 0, linetype=1) +
      geom_point(aes(color=legend, shape=legend), size=2, position=position_jitter(0.1)) +
      labs(title="Myotubes exposed to palmitate",
           x="",
           y="Palmitate-induced mRNA\nlog2(relative to control)") +
      theme_bw() + theme + theme(legend.title=element_blank()) +
      add_pvalue(stat.test, bracket.size = NA, 
                 xmin = "gene",
                 xmax = "gene",
                 label = "p.adj.signif",
                 y.position = "y.position", label.size=10) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))
    plot_C2C12
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