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
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Transcriptomic response of skeletal myotubes to palmitate"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), "/ last update 2022-07-07")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 0% 8%;",
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
                tags$hr(),
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                         column(2, checkboxGroupInput("concentration", 
                                                      label = "Concentration (µM)", 
                                                      selected = c(100,
                                                                   200,
                                                                   400,
                                                                   500), 
                                                      choices = c(100,
                                                                  200,
                                                                  400,
                                                                  500))),
                         column(2, checkboxGroupInput("time", 
                                                      label = "Duration (h)", 
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
                                                                  54))),
                         column(2, checkboxGroupInput("model", 
                                                      label = "Cell type", 
                                                      selected = c("C2C12",
                                                                   "HSMC",
                                                                   "LHCN-M2"), 
                                                      choices = c("C2C12",
                                                                  "HSMC",
                                                                  "LHCN-M2"))),
                         actionButton("updatePlot", "Refresh plot", icon("refresh"))
                ),
                fluidRow(style="color:black;background-color:white;padding:2% 8% 1% 8%;",
                         plotOutput("PalPlot", height="400px") %>% withSpinner(color="#5b768e"),
                         h5("Statistics presented on the plot are paired t-test (palmitate vs BSA-vehicle)
                            with Bonferroni correction for multiple testing of 12 265 genes in the dataset.
                            Statistics are dynamically calculated based on the selection criteria.")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         tags$hr(),
                         h3("Statistics (calculated with limma on all samples)"),
                         dataTableOutput("stats")
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
    genename <- c("GPSM2", "ANGPTL4", "HMOX1", "CCL2", "IL6")
    genename <- input$inputGeneSymbol
    
    #plot
    testdata <- data.frame()
    for( i in genename) { 
      data <- as.numeric(PAL_data[i,])              #collect data for gene name i
      datay <- cbind.data.frame(Sample_norm[,2:8], data, gene=rep(i))   #create table with x="sample type", y="data", "gene name"
      testdata  <- rbind.data.frame(testdata, datay)              #bind the new gene data at the bottom of the previous one
    }
    
    testdata$gene <- factor(testdata$gene,
                            levels=unique(str_sort(testdata$gene, numeric = TRUE)))
    
    #filter according to selected categories
    testdata <- testdata[testdata$concentration %in% input$concentration &
                           testdata$time %in% input$time &
                           testdata$model %in% input$model 
                         ,]
    
    #get stats from ggpubr
    stat.test <- compare_means(data=na.omit(testdata), 
                               data ~ treatment, 
                               group.by = "gene",
                               paired = T)
    
    #adjust for multiple testing
    stat.test$p.adj <- p.adjust(stat.test$p, 
                                method="bonferroni",
                                n=nrow(PAL_stats))
    
    #format FDR
    stat.test$p.adj.format <- ifelse(stat.test$p.adj > 0.05,
           "ns",
           paste0(
             "p.adj =\n",
             format(signif(stat.test$p.adj, 2), scientific = TRUE)
             )
           )
    
    #add significance
    stat.test$p.adj.signif <- "ns"
    stat.test$p.adj.signif[stat.test$p.adj<0.05] <- "*"
    stat.test$p.adj.signif[stat.test$p.adj<0.01] <- "**"
    stat.test$p.adj.signif[stat.test$p.adj<0.001] <- "***"
    
    #plot only treated condition relative to each basal
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
      geom_hline(yintercept = 0, 
                 linetype=1) +
      geom_point(aes(color=description), 
                 size=2, 
                 position=position_jitter(0.1),
                 alpha=0.2) +
      labs(subtitle = "",
           x="",
           y="Palmitate-induced mRNA\nlog2(relative to control)") +
      theme_bw() + theme + theme(legend.title=element_blank()) +
      add_pvalue(stat.test, bracket.size = NA, 
                 xmin = "gene",
                 xmax = "gene",
                 label = "p.adj.format",
                 y.position = "y.position", 
                 label.size = 3.5) +
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