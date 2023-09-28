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
# EPS
#--------------------------------------------------------------------------------------------------------
EPS_data <- readRDS("data/data.Rds")

#--------------------------------------------------------------------------------------------------------
# List of genes
#--------------------------------------------------------------------------------------------------------
genelist <- unique(rownames(EPS_data))


#--------------------------------------------------------------------------------------------------------
# Sample
#--------------------------------------------------------------------------------------------------------
#find sample list
Sample_norm <- data.frame(sample=colnames(EPS_data),
                          str_split_fixed(colnames(EPS_data), "_|\\.", 6))
colnames(Sample_norm) <- c("sample", "GEO", "model", "time", "conditions", "treatment","repeats")

Sample_norm$time <- as.numeric(gsub("H", "", Sample_norm$time))

Sample_norm$repeats <- paste(Sample_norm$GEO,
                             Sample_norm$time,
                             Sample_norm$repeats)

#add description of stidues
Sample_norm$description <- paste0(
  Sample_norm$GEO, ", ",
  Sample_norm$model, ", ",
  Sample_norm$conditions, ", ",
  Sample_norm$time, "h"
)


# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Transcriptomic response of skeletal myotubes to electrical pulse stimulation"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), "/ last update 2022-07-08")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 0% 8%;",
                         "From",
                         a("GSE44051", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44051", 
                           target="_blank", style="color:#5B768E"), 
                         "(HSMC, n=12),",
                         a("GSE139872", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139872", 
                           target="_blank", style="color:#5B768E"), 
                         "(C2C12, n=4),",
                         a("GSE183159", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183159", 
                           target="_blank", style="color:#5B768E"), 
                         "(C2C12, n=1),",
                         a("GSE185530", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185530", 
                           target="_blank", style="color:#5B768E"), 
                         "(C2C12, n=6),",
                         a("GSE200335", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200335", 
                           target="_blank", style="color:#5B768E"), 
                         "(HSMC, n=7),", 
                         a("GSE201340", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201340", 
                           target="_blank", style="color:#5B768E"), 
                         "(HSMC, n=4)."


                ),
                tags$hr(),
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                         column(2, checkboxGroupInput("time", 
                                                      label = "Duration (h)", 
                                                      selected = c(3,
                                                                   24,
                                                                   48), 
                                                      choices = c(3,
                                                                  24,
                                                                  48))),
                         column(2, checkboxGroupInput("model", 
                                                      label = "Cell type", 
                                                      selected = c("C2C12",
                                                                   "HSMC"), 
                                                      choices = c("C2C12",
                                                                  "HSMC"))),
                         actionButton("updatePlot", "Refresh plot", icon("refresh"))
                ),
                fluidRow(style="color:black;background-color:white;padding:2% 8% 1% 8%;",
                         plotOutput("PalPlot", height="400px") %>% withSpinner(color="#5b768e"),
                         h5("Statistics presented on the plot are paired t-test (EPS vs control)
                            with Bonferroni correction for multiple testing. Statistics are dynamically calculated 
                            based on the selection criteria.")
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected=c("TNN", "PLEKHG3", "HMOX1", "OMD", "IL6", "NR4A3"), 
                       options=NULL)
  
  # AMPK plot
  plotData <- eventReactive(input$updatePlot, {
    validate(need(input$inputGeneSymbol, " "))
    genename <- c("GPSM2", "ANGPTL4", "HMOX1", "CCL2", "IL6")
    genename <- input$inputGeneSymbol
    
    #plot
    testdata <- data.frame()
    for( i in genename) { 
      data <- as.numeric(EPS_data[i,])              #collect data for gene name i
      datay <- cbind.data.frame(Sample_norm[,2:8], data, gene=rep(i))   #create table with x="sample type", y="data", "gene name"
      testdata  <- rbind.data.frame(testdata, datay)              #bind the new gene data at the bottom of the previous one
    }
    
    testdata$gene <- factor(testdata$gene,
                            levels=unique(str_sort(testdata$gene, numeric = TRUE)))
    
    #filter according to selected categories
    testdata <- testdata[testdata$time %in% input$time &
                           testdata$model %in% input$model 
                         ,]
    
    #get stats from ggpubr
    stat.test <- compare_means(data = testdata, 
                               data ~ treatment, 
                               method = "wilcox.test",
                               group.by = "gene",
                               paired = TRUE)
    
    #adjust for multiple testing
    stat.test$p.adj <- p.adjust(stat.test$p, 
                                method="bonferroni",
                                n=nrow(EPS_data)/4)
    
    #format stats
    stat.test$p.format <- paste0("p =\n",
                                     format(signif(stat.test$p, 2), scientific = TRUE)
    )
    stat.test$p.adj.format <- paste0("p.adj =\n",
                  format(signif(stat.test$p.adj, 2), scientific = TRUE)
                  )
    
    #add significance
    stat.test$p.adj.signif <- "ns"
    stat.test$p.adj.signif[stat.test$p.adj<0.05] <- "*"
    stat.test$p.adj.signif[stat.test$p.adj<0.01] <- "**"
    stat.test$p.adj.signif[stat.test$p.adj<0.001] <- "***"
    
    #plot only treated condition relative to each basal
    foldchange <- testdata[testdata$treatment %in% "EPS",]
    foldchange$data <- testdata[testdata$treatment %in% "EPS",]$data - testdata[testdata$treatment %in% "CTRL",]$data 
    
    #manually get y position for each gene
    max.y <- data.frame(aggregate(foldchange$data, by = list(foldchange$gene), max, na.rm=T))
    colnames(max.y) <- c("gene", "y.position")
    max.y$y.position <- max.y$y.position + 0.5
    stat.test <- full_join(stat.test, max.y)
    
    #remove ns
    stat.test$p.adj.signif[stat.test$p.adj.signif=="ns"] <- ""

    #bar plot
    plot_EPS <- ggbarplot(foldchange, x="gene", y="data", add="mean_se", fill="gray80", size=0.5) +  
      geom_hline(yintercept = 0, 
                 linetype=1) +
      geom_point(aes(color = description, shape = description), 
                 size=2, 
                 position=position_jitter(0.1),
                 alpha=0.8) +
      labs(subtitle = "",
           x="",
           y="EPS-induced mRNA\nlog2(relative to control)") +
      theme_bw() + theme + theme(legend.title=element_blank()) +
      add_pvalue(stat.test, bracket.size = NA, 
                 xmin = "gene",
                 xmax = "gene",
                 label = "p.adj.format",
                 y.position = "y.position", 
                 label.size = 3.5) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))
    plot_EPS
  })
  
  
  output$PalPlot <- renderPlot({
    plotData()
  })
  
}


# Run the app ----
shinyApp(ui = ui, server = server)