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
library(forestplot)
library(metafor)
library(stringr)
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
MuscleClamp_stats_MetaAnalysis <- readRDS('data/stats_MetaAnalysis.Rds')
MuscleClamp_stats_TimeCouse <- readRDS('data/stats_TimeCourse.Rds')

GSE22309_data <- readRDS('data/GSE22309_data.Rds')
GSE22309_stats_IS <- readRDS('data/GSE22309_stats_IS.Rds')
colnames(GSE22309_stats_IS) <- gsub("_GSE.*", "", colnames(GSE22309_stats_IS))
GSE22309_stats_IR <- readRDS('data/GSE22309_stats_IR.Rds')
colnames(GSE22309_stats_IR) <- gsub("_GSE.*", "", colnames(GSE22309_stats_IR))
GSE22309_stats_DIA <- readRDS('data/GSE22309_stats_DIA.Rds')
colnames(GSE22309_stats_DIA) <- gsub("_GSE.*", "", colnames(GSE22309_stats_DIA))
GSE22309_stats_intDIA <- readRDS('data/GSE22309_stats_int.DIA.Rds')
GSE22309_stats_intIR <- readRDS('data/GSE22309_stats_int.IR.Rds')

#genes
genelist <- rownames(GSE22309_data)

#function to extract data for a gene
DataForGeneName <- function(x){
  x <- data.frame(t(x[grepl('logFC',    colnames(x))]), # M-value (M) is the log2-fold change
                  t(x[grepl('adj.P.Val',colnames(x))]), # Benjamini and Hochberg's method to control the false discovery rate
                  t(x[grepl('CI.L',     colnames(x))]), # lower limit of the 95% confidence interval
                  t(x[grepl('CI.R',     colnames(x))]), # upper limit of the 95% confidence interval
                  t(x[grepl('mean.pre', colnames(x))]), # mean of control condition
                  t(x[grepl('mean.post', colnames(x))]), # mean of exercise condition
                  t(x[grepl('Sd.pre',   colnames(x))]), # standard deviation of control condition
                  t(x[grepl('Sd.post',   colnames(x))]), # standard deviation of exercise condition
                  t(x[grepl('size',     colnames(x))])) # number of subjects in the study
  x <- cbind(x, str_split_fixed(rownames(x), "_", 11))
  colnames(x) <- c('logFC', 'adj.P.Val', 
                   'CI.L', 'CI.R',
                   'Mean_Ctrl', 'Mean_Ex', 
                   'Sd_Ctrl', 'Sd_Ex', 'size',
                   'Studies', 'GEO', 'Exercisetype', 
                   'Muscle', 'Sex', 'Age', 'Training',
                   'Obesity', 'Disease', 'Biopsy', 'Duration')
  x$Studies <- gsub("logFC_","", rownames(x))
  x
}

#function for meta-analysis calculation
MetaAnalysis <- function(x, nrow){
  #order by time
  x <- x[c("logFC_GSE9105_HLY_30", 
           "logFC_GSE7146A_HLY_120", 
           "logFC_GSE7146B_HLY_180", 
           "logFC_GSE22309_HLY_180", 
           "logFC_GSE22309_IR_180", 
           "logFC_GSE22309_DIA_180", 
           "logFC_GSE9105_HLY_240"),]
  #if only one row, skip calculation
  if(nrow(x)<2){
    x <- rbind(x[,1:10],
               c(as.numeric(x[,1:4]), rep(NA, 4), sum(x$size, na.rm=T)))
    x$Studies <- c(gsub("logFC_", "", x$Studies[1:(nrow(x)-1)]),
                   "Meta-analysis score")
    return(x)
    
  } else {
    #Calculate metascore
    meta <- metafor::rma(m1 = Mean_Ex, m2 = Mean_Ctrl,
                         sd1 = Sd_Ex, sd2 = Sd_Ctrl,
                         n1 = size, n2 = size,
                         method = "REML",
                         measure = "MD",
                         data = x,
                         control=list(maxiter=1000, stepadj=0.5),
                         weighted=T, weights=x$size)
    fdr  <- p.adjust(meta$pval, method='bonferroni', n=nrow)
    x <- rbind(x[,1:10],
               c(meta$beta, fdr, meta$ci.lb, meta$ci.ub, rep(NA, 4), sum(x$size, na.rm=T)))
    x$Studies <- c(gsub("logFC_", "", x$Studies[1:(nrow(x)-1)]),
                   "Meta-analysis score")
    return(x)
  }
}



#function to make plots
ModuleForestPlot <- function(metadata, genename, color, title) {
  metadata <- data.frame(metadata)
  max_n <- max(metadata[1:(nrow(metadata)-1),]$size, na.rm=T)+1 #used to plot weights
  tabledata <- cbind(mean = c(NA , metadata[,1]), 
                     lower= c(NA , metadata[,3]),
                     upper= c(NA , metadata[,4]))
  tabletext <- cbind(c(paste(genename, title, sep=' ') , metadata[,10]),
                     c("logFC" , format(round(metadata[,1], digits=2))), 
                     c("FDR"   , format(metadata[,2],   scientific=T, digits=2)),
                     c("n" , metadata[,9]))
  finalplot <- forestplot(tabletext, grid=T,
                          tabledata, new_page=TRUE,
                          is.summary=c(TRUE,rep(FALSE,(nrow(metadata)-1)),TRUE),
                          xlog=F,  txt_gp = own, xlab="logFC",
                          zero = 0, lwd.zero=1,
                          col=fpColors(box=color, line=color, summary=color),
                          boxsize=c(NA, metadata$size[1:(nrow(metadata)-1)]/max_n, 1)) #Make the size of the symbols relative to sample size
  return(finalplot)
}

#graphical settings for forest plots
own <- fpTxtGp()
own$label$cex <- 0.9
own$summary$cex <- 1.1
own$title$cex <- 1.4
own$xlab$cex <- 0.9
own$ticks$cex <- 0.8

# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Transcriptomic response of skeletal muscle to hyperinsulinemic euglycemic clamp"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), 
                            "/ last update 2021-09-12")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 0% 8%;",
                         "Skeletal muscle data from healthy volonteers from",
                         a("GSE7146", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7146", 
                           target="_blank", style="color:#5B768E"),
                         "(120min and 180min, n=6) and",
                         a("GSE9105", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9105", 
                           target="_blank", style="color:#5B768E"), 
                         "(30min and 240min, n=12). Data of skeletal muscle response to clamp in healthy, 
                         insulin resistant and diabetic individuals from",
                         a("GSE22309", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22309", 
                           target="_blank", style="color:#5B768E"), 
                         "(180min, n=54).",
                         tags$hr(),
                         selectizeInput("inputGeneSymbol", "Gene Symbol:", choices=NULL, multiple=F, width=300),
                ),
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         column(7,
                                plotOutput("metaPlot", height="250px") %>% withSpinner(color="#5b768e")
                         ),
                         column(5,
                                plotOutput("timePlot", height="250px") %>% withSpinner(color="#5b768e")
                         )
                ),
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         column(7,
                                plotOutput("insPlot", height="350px") %>% withSpinner(color="#5b768e")
                         ),
                         column(4,
                                plotOutput("intPlot", height="350px") %>% withSpinner(color="#5b768e")
                         )
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected="HES1", 
                       options=NULL)
  
  output$metaPlot <- renderPlot({
    validate(need(input$inputGeneSymbol, " "))
    genename <- "HES1"
    genename <- input$inputGeneSymbol
    
    selectedata <- DataForGeneName(MuscleClamp_stats_MetaAnalysis[genename,])
    metadata <- MetaAnalysis(selectedata, nrow(MuscleClamp_stats_MetaAnalysis))
    ModuleForestPlot(metadata, genename, "#9D6807", "during clamp")
  })
  
  
  output$insPlot <- renderPlot({
    validate(need(input$inputGeneSymbol, " "))
    genename <- "HES1"
    genename <- input$inputGeneSymbol
    
    #make table for plots
    testdata <- data.frame(t(data.frame(str_split(colnames(GSE22309_data), "_"))),
                           data=as.numeric(GSE22309_data[genename,])  )
    colnames(testdata) <- c("GEO", "diagnosis", "time", "subject", "data")
    rownames(testdata) <- colnames(GSE22309_data)
    
    #order groups
    testdata$time <- factor(testdata$time, levels=c("PRE", "180"))
    
    #get stats from limma
    stat.limma <- rbind(GSE22309_stats_IS[genename,],
                        GSE22309_stats_IR[genename,],
                        GSE22309_stats_DIA[genename,])
    stat.limma
    
    #get stats from ggpubr
    stat.test <- compare_means(data ~ time, testdata, group.by="diagnosis")
    stat.test
    
    #replace stats by the ones from limma
    stat.test$p <- stat.limma$P.Value
    stat.test$p.adj <- stat.limma$adj.P.Val
    stat.test$p.adj.signif <- "ns"
    stat.test$p.adj.signif[stat.test$p.adj<0.1] <- "*"
    stat.test$p.adj.signif[stat.test$p.adj<0.05] <- "**"
    stat.test$p.adj.signif[stat.test$p.adj<0.01] <- "***"
    
    # Add position for adjusted p-values
    stat.test <- stat.test %>% add_xy_position(data=testdata, x = "time", formula=data ~ time,  group="diagnosis")
    stat.test$y.position <- stat.test$y.position * 1.05
    
    #plot
    plot_clamp <- ggboxplot(testdata, x = "diagnosis", y = "data", fill="time", outlier.size=0.5) + 
      theme_bw(13) + theme(legend.title = element_blank()) +
      labs(title="Response to clamp",
           x=element_blank(),
           y=paste(genename, ", log2(relative mRNA)")) +
      scale_fill_manual(values=c("gray80", "gray40")) +
      stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = T, x="diagnosis", size=8) 
    plot_clamp
    
  })
  
  
  output$intPlot <- renderPlot({
    validate(need(input$inputGeneSymbol, " "))
    genename <- "HES1"
    genename <- input$inputGeneSymbol
    
    res <- GSE22309_data
    
    #caculate ratios
    res.PRE <- res[grepl("PRE", colnames(res))]
    res.180 <- res[grepl("180", colnames(res))]
    res <- res.180 - res.PRE
    
    #make table for plots
    testdata <- data.frame(t(data.frame(str_split(colnames(res), "_"))),
                           data=as.numeric(res[genename,])  )
    colnames(testdata) <- c("GEO", "diagnosis", "time", "subject", "data")
    rownames(testdata) <- colnames(res)
    testdata$diagnosis <- factor(testdata$diagnosis, levels=c("HLY", "IR", "DIA"))
    
    #get stats from limma
    stat.limma <- rbind(GSE22309_stats_intIR[genename,],
                        GSE22309_stats_intDIA[genename,])
    stat.limma
    
    #get stats from ggpubr
    stat.test <- t_test(testdata, data ~ diagnosis, ref.group = "HLY")
    stat.test
    
    #replace stats by the ones from limma
    stat.test$p <- stat.limma$P.Value
    stat.test$p.adj <- stat.limma$adj.P.Val
    stat.test$p.adj.signif <- "ns"
    stat.test$p.adj.signif[stat.test$p.adj<0.05] <- "*"
    stat.test$p.adj.signif[stat.test$p.adj<0.01] <- "**"
    stat.test$p.adj.signif[stat.test$p.adj<0.001] <- "***"
    
    # Add position for adjusted p-values
    stat.test <- stat.test %>% add_xy_position(x = "diagnosis", dodge=0.8)
    stat.test$y.position <- stat.test$y.position * 1.05
    
    #plot
    plot_clamp <- ggboxplot(testdata, x = "diagnosis", y = "data", fill="diagnosis", outlier.size=0.5) + 
      theme_bw(13) + theme(legend.title = element_blank()) +
      labs(title="Interaction (health x clamp)",
           x=element_blank(),
           y=paste(genename, "mRNA, relative to basal")) +
      scale_fill_manual(values=c("#009E73", "#E69F00", "#D55E00")) +
      stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = T, size=8) 
    plot_clamp
  })
  
  
  
  
  output$timePlot <- renderPlot({
    validate(need(input$inputGeneSymbol, " "))
    genename <- "HES1"
    genename <- input$inputGeneSymbol
    
    #data for plot
      testdata <- data.frame(time=c(0,30,120,180,240),
                        mean=as.numeric(MuscleClamp_stats_TimeCouse[genename, 
                                                                    grepl("Median", colnames(MuscleClamp_stats_TimeCouse))]),
                        sem=as.numeric(MuscleClamp_stats_TimeCouse[genename, 
                                                                   grepl("SEM", colnames(MuscleClamp_stats_TimeCouse))]),
                        pvalue=c(1, as.numeric(MuscleClamp_stats_TimeCouse[genename, 
                                                                           grepl("P.Value", colnames(MuscleClamp_stats_TimeCouse))])),
                        fdr=c(1, as.numeric(MuscleClamp_stats_TimeCouse[genename, 
                                                                        grepl("adj.P.Val", colnames(MuscleClamp_stats_TimeCouse))])))
    
    #add significance level
    testdata$sign <- ""
    testdata$sign[testdata$fdr<0.05] <- "*"
    testdata$sign[testdata$fdr<0.01] <- "**"
    testdata$sign[testdata$fdr<0.001] <- "***"
    testdata$y.position <- (testdata$mean+testdata$sem)+0.05
    
    #plot
    ggplot(testdata, aes(x=time, y=mean)) +
      theme_bw(13) +
      geom_hline(yintercept=0, linetype=3, color='grey30', size=1) +
      geom_point() +   geom_line() +
      geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1) +
      labs(x="Time (min)",
           y="mRNA expression\nlog2(relative to basal)") +
      add_pvalue(testdata, bracket.size = NA,
                 xmin = "time",
                 xmax = "time",
                 label = "sign",
                 y.position = "y.position", label.size=4)
  })
  
  
  
  
  
  
  ##################################################################################################################
  #Statistic tables
  output$stats_sex <- renderDataTable(options=list(signif = 3),{
    datatable(MuscleClamp_sex_stats[,c(1,5,2,7,8)], 
              options = list(lengthMenu = c(5, 10, 50, 100), 
                             pageLength = 5),
              rownames= FALSE) %>% 
      formatSignif(columns = c(2:5), digits = 3)
  })
  
  
}


# Run the app ----
shinyApp(ui = ui, server = server)
