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
               legend.title = element_blank(),
               legend.position="bottom",
               legend.key.size = unit(18, "pt"))
library(stringr)
library(grid)
library(DT)
library(gridExtra)
library(ggpubr)
library(rstatix)

datasets <- readRDS("data/datasets.Rds")


#Mouse data
mouse_data <- readRDS("data/mouse_data.Rds")
mouse_means <- readRDS("data/mouse_means.Rds")
mouse_M1vsM0 <- readRDS("data/mouse_M1vsM0.Rds")
mouse_M2vsM0 <- readRDS("data/mouse_M2vsM0.Rds")
mouse_M1vsM2 <- readRDS("data/mouse_M1vsM2.Rds")

mouse_stats_table <- data.frame(SYMBOL=rownames(mouse_M1vsM0),
                                M1vsM0_logFC=mouse_M1vsM0$logFC,
                                M1vsM0_FDR=mouse_M1vsM0$adj.P.Val,
                                M2vsM0_logFC=mouse_M2vsM0$logFC,
                                M2vsM0_FDR=mouse_M2vsM0$adj.P.Val,
                                M1vsM2_logFC=mouse_M1vsM2$logFC,
                                M1vsM2_FDR=mouse_M1vsM2$adj.P.Val)
mouse_stats_table <- mouse_stats_table[order(mouse_stats_table$M1vsM2_FDR, decreasing = F),]


#Human data
human_data <- readRDS("data/human_data.Rds")
human_means <- readRDS("data/human_means.Rds")
human_M1vsM0 <- readRDS("data/human_M1vsM0.Rds")
human_M2vsM0 <- readRDS("data/human_M2vsM0.Rds")
human_M1vsM2 <- readRDS("data/human_M1vsM2.Rds")

human_stats_table <- data.frame(SYMBOL=rownames(human_M1vsM0),
                                M1vsM0_logFC=human_M1vsM0$logFC,
                                M1vsM0_FDR=human_M1vsM0$adj.P.Val,
                                M2vsM0_logFC=human_M2vsM0$logFC,
                                M2vsM0_FDR=human_M2vsM0$adj.P.Val,
                                M1vsM2_logFC=human_M1vsM2$logFC,
                                M1vsM2_FDR=human_M1vsM2$adj.P.Val)
human_stats_table <- human_stats_table[order(human_stats_table$M1vsM2_FDR, decreasing = F),]


#Human WAT data
human_WAT_data <- readRDS("data/human_WAT_data.Rds")
human_WAT_means <- readRDS("data/human_WAT_means.Rds")
human_WAT_M1vsM0 <- readRDS("data/human_WAT_M1vsM0.Rds")
human_WAT_M2vsM0 <- readRDS("data/human_WAT_M2vsM0.Rds")
human_WAT_M1vsM2 <- readRDS("data/human_WAT_M1vsM2.Rds")

human_WAT_stats_table <- data.frame(SYMBOL=rownames(human_WAT_M1vsM0),
                                    M1vsM0_logFC=human_WAT_M1vsM0$logFC,
                                    M1vsM0_FDR=human_WAT_M1vsM0$adj.P.Val,
                                    M2vsM0_logFC=human_WAT_M2vsM0$logFC,
                                    M2vsM0_FDR=human_WAT_M2vsM0$adj.P.Val,
                                    M1vsM2_logFC=human_WAT_M1vsM2$logFC,
                                    M1vsM2_FDR=human_WAT_M1vsM2$adj.P.Val)
human_WAT_stats_table <- human_WAT_stats_table[order(human_WAT_stats_table$M1vsM2_FDR, decreasing = F),]


#gene list
genelist <- unique(c(rownames(mouse_data), rownames(human_data), rownames(human_WAT_data)))

#function to gene name
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Transcriptomic profile of human and mouse macrophages"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), "/ last update 2022-04-07")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         selectizeInput("inputGeneSymbol", "Gene Symbol:", choices=NULL, multiple=F),
                         sidebarLayout(
                           sidebarPanel(width = 3, style=("font-size:75%"),
                                        
                                        #checkboxGroupInput("treatment_M1", 
                                        #                   label="M1 treatment", 
                                        #                   selected=c("LPS", "IFNg", "LPS+IFNg"), 
                                        #                   c("LPS", "IFNg", "LPS+IFNg")),
                                        #checkboxGroupInput("treatment_M2", 
                                        #                   label="M2 treatment", 
                                        #                   selected=c("IL4", "IL13", "IL4+IL13"), 
                                        #                   c("IL4", "IL13", "IL4+IL13")),
                                        # tags$hr(),
                                        h5(helpText("Mouse Primary")),
                                        tableOutput("mouse_means"),
                                        h5(helpText("Human Primary")),
                                        tableOutput("human_means"),
                                        h5(helpText("Human WAT")),
                                        tableOutput("human_WAT_means")
                           ),
                           
                           mainPanel(width = 9,
                                     fluidRow(
                                       column(4, align="left",
                                              plotOutput("mousePlot", height="500px") %>% withSpinner(color="#5b768e")),
                                       column(4, align="left",
                                              plotOutput("humanPlot", height="500px") %>% withSpinner(color="#5b768e")),
                                       column(4, align="left",
                                              plotOutput("FACSPlot", height="500px") %>% withSpinner(color="#5b768e"))
                                     )
                           )
                         )
                ),
                
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         tags$hr(),
                         h3("Statistics - mouse primary macrophages"),
                         dataTableOutput("mouse_stats_table"),
                         tags$hr(),
                         h3("Statistics - human primary macrophages"),
                         dataTableOutput("human_stats_table"),
                         tags$hr(),
                         h3("Statistics - human adipose tissue macrophages"),
                         dataTableOutput("human_WAT_stats_table")
                ),
                
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         tags$hr(),
                         h3("Datasets"),
                         dataTableOutput("datasets")
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected='F13A1', 
                       options=NULL)
  
  ##################################################################################################################
  mousePlot_function <- function(){
    genename <- "F13a1"
    genename <- firstup(tolower(input$inputGeneSymbol))
    data <- mouse_data[genename,]
    validate(need(isFALSE(rowSums(is.na(data)) == ncol(data)), "\nImpossible to find data for this gene name."))
    
    
    #get samples
    genedata <- t(mouse_data[genename,])
    samples <- gsub("_.*", "", colnames(mouse_data))
    GEO <- gsub(".*_", "", colnames(mouse_data))
    GEO <- gsub("\\..*", "", GEO)
    treatment <- gsub("M1_", "", colnames(mouse_data))
    treatment <- gsub("M2_", "", treatment)
    treatment <- gsub("_.*", "", treatment)
    treatment <- gsub("\\.", "+", treatment)
    treatment <- gsub("M0", "Control", treatment)
    treatment <- factor(treatment, levels=c("Control", "LPS", "IFNg", "LPS+IFNg", "IL4", "IL13", "IL4+IL13"))
    testdata <- data.frame(GEO, samples, treatment, genedata)
    colnames(testdata) <- c('GEO', 'sample', 'Treatment', 'gene')
    testdata$sample <- factor(testdata$sample, levels=c("M0", "M1", "M2"))
    
    binwidth <- max(testdata$gene, na.rm=T) - min(testdata$gene, na.rm=T)
    
    #filter according to selected categories
    #testdata <- testdata[testdata$Treatment %in% c("Control", input$treatment_M1, input$treatment_M2),]
    
    #normalize to M0
    #M0mean <- median(testdata[testdata$sample %in% 'M0','gene'], na.rm=T)
    #testdata$gene <- testdata$gene - M0mean
    
    #get stats from limma
    stat.limma <- rbind(mouse_M1vsM0[genename,],
                        mouse_M2vsM0[genename,],
                        mouse_M1vsM2[genename,])

    #get stats from ggpubr
    stat.test <- t_test(testdata, gene ~ sample)
    
    #replace stats by the ones from limma
    stat.test$p <- stat.limma$P.Value
    stat.test$p.adj <- stat.limma$adj.P.Val
    stat.test$p.adj.signif <- "ns"
    stat.test$p.adj.signif[stat.test$p.adj<0.05] <- "*"
    stat.test$p.adj.signif[stat.test$p.adj<0.01] <- "**"
    stat.test$p.adj.signif[stat.test$p.adj<0.001] <- "***"

    # Add position for adjusted p-values
    stat.test <- stat.test %>% add_xy_position(x = "sample")
    stat.test$y.position <- stat.test$y.position * 1.05
    
    #boxplot
    plot_mouse <- ggboxplot(testdata, x = "sample", y = "gene") + 
      theme_bw() + theme + guides(fill=guide_legend(nrow=4)) +
      labs(title="Mouse\nprimary macrophages",
           x=element_blank(),
           y=paste(genename, "log2(relative mRNA)")) +
      geom_boxplot() +
      geom_dotplot(aes(x=sample, y=gene, fill=Treatment),
                   binwidth = binwidth/60, alpha=0.8,
                   binaxis = "y", stackdir = "center", binpositions="all") +
      scale_shape_manual(values=c(15,16,17,15,16,17,15,16,17,15,16,17,15,16,17,15,16,17)) +
      scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#CC79A7", "#009E73", "#0072B2", "#D55E00", "#CC79A7")) +
      stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = T, size=7, bracket.size = 0.5) 
    plot_mouse
    return(plot_mouse)
  }
  
  output$mousePlot <- renderPlot({
    mousePlot_function()
  })
  
  output$mouse_means <- renderTable(align="c", spacing="xs", rownames=T, colnames=T,{
    genename <- firstup(tolower(input$inputGeneSymbol))
    data.frame(mean=as.character(round(as.numeric(mouse_means[genename,1:3]), 2)),
               Sd=as.character(round(as.numeric(mouse_means[genename,4:6]), 2)),
               n=as.character(mouse_means[genename,7:9]),
               row.names = c("M0", "M1", "M2"))
  })
  
  
  ##################################################################################################################
  humanPlot_function <- function(){
    #genename <- "Arg1"
    genename <- toupper(input$inputGeneSymbol)
    data <- human_data[genename,]
    validate(need(isFALSE(rowSums(is.na(data)) == ncol(data)), "\nImpossible to find data for this gene name."))
    
    #get samples
    genedata <- t(human_data[genename,])
    samples <- gsub("_.*", "", colnames(human_data))
    GEO <- gsub(".*_", "", colnames(human_data))
    GEO <- gsub("\\..*", "", GEO)
    treatment <- gsub("M1_", "", colnames(human_data))
    treatment <- gsub("M2_", "", treatment)
    treatment <- gsub("_.*", "", treatment)
    treatment <- gsub("\\.", "+", treatment)
    treatment <- gsub("M0", "Control", treatment)
    treatment <- factor(treatment, levels=c("Control", "LPS", "IFNg", "LPS+IFNg", "IL4", "IL13", "IL4+IL13"))
    testdata <- data.frame(GEO, samples, treatment, genedata)
    colnames(testdata) <- c('GEO', 'sample', 'Treatment', 'gene')
    testdata$sample <- factor(testdata$sample, levels=c("M0", "M1", "M2"))
    
    binwidth <- max(testdata$gene, na.rm=T) - min(testdata$gene, na.rm=T)
    
    #filter according to selected categories
    #testdata <- testdata[testdata$Treatment %in% c("Control", input$treatment_M1, input$treatment_M2),]
    
    #normalize to M0
    #M0mean <- median(testdata[testdata$sample %in% 'M0','gene'], na.rm=T)
    #testdata$gene <- testdata$gene - M0mean
    
    #get stats from limma
    stat.limma <- rbind(human_M1vsM0[genename,],
                        human_M2vsM0[genename,],
                        human_M1vsM2[genename,])
    
    #get stats from ggpubr
    stat.test <- t_test(testdata, gene ~ sample)
    
    #replace stats by the ones from limma
    stat.test$p <- stat.limma$P.Value
    stat.test$p.adj <- stat.limma$adj.P.Val
    stat.test$p.adj.signif <- "ns"
    stat.test$p.adj.signif[stat.test$p.adj<0.05] <- "*"
    stat.test$p.adj.signif[stat.test$p.adj<0.01] <- "**"
    stat.test$p.adj.signif[stat.test$p.adj<0.001] <- "***"
    
    # Add position for adjusted p-values
    stat.test <- stat.test %>% add_xy_position(x = "sample")
    stat.test$y.position <- stat.test$y.position * 1.05
    
    plot_human <- ggboxplot(testdata, x = "sample", y = "gene") + 
      theme_bw() + theme + guides(fill=guide_legend(nrow=4)) +
      labs(title="Human\nprimary macrophages",
           x=element_blank(),
           y=paste(genename, "log2(relative mRNA)")) +
      geom_boxplot() +
      geom_dotplot(aes(x=sample, y=gene, fill=Treatment),
                   binwidth = binwidth/60, alpha=0.8,
                   binaxis = "y", stackdir = "center", binpositions="all") +
      scale_shape_manual(values=c(15,16,17,15,16,17,15,16,17,15,16,17,15,16,17,15,16,17)) +
      scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#CC79A7", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
      stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = T, size=7, bracket.size = 0.5) 
    plot_human
    return(plot_human)
  }
  
  output$humanPlot <- renderPlot({
    humanPlot_function()
  })
  
  output$human_means <- renderTable(align="c", spacing="xs", rownames=T, colnames=T,{
    genename <- toupper(input$inputGeneSymbol)
    data.frame(mean=as.character(round(as.numeric(human_means[genename,1:3]), 2)),
               Sd=as.character(round(as.numeric(human_means[genename,4:6]), 2)),
               n=as.character(human_means[genename,7:9]),
               row.names = c("M0", "M1", "M2"))
  })
  
  
  ##################################################################################################################
  FACSPlot_function <- function(){
    genename <- "CD209"
    genename <- toupper(input$inputGeneSymbol)
    data <- human_WAT_data[genename,]
    validate(need(isFALSE(rowSums(is.na(data)) == ncol(data)), "\nImpossible to find data for this gene name."))
    
    #get samples
    genedata <- t(human_WAT_data[genename,])
    samples <- gsub("\\..*", "", colnames(human_WAT_data))
    treatment <- gsub(" .*", "", colnames(human_WAT_data))
    treatment <- gsub("CD14\\+", "All", treatment)
    testdata <- data.frame(samples, treatment, genedata)
    colnames(testdata) <- c('sample', 'Treatment', 'gene')
    binwidth <- max(testdata$gene, na.rm=T) - min(testdata$gene, na.rm=T)
    
    #add markers
    testdata$Markers <- "CD45/14/206+"
    testdata$Markers[testdata$Treatment == "M1"] <- "CD45/14/206+/CD11c+"
    testdata$Markers[testdata$Treatment == "M2"] <- "CD45/14/206+/CD11c-"
    testdata$Markers <- factor(testdata$Markers, levels=c("CD45/14/206+", "CD45/14/206+/CD11c+", "CD45/14/206+/CD11c-"))

    #normalize to M0
    #M0mean <- median(testdata[testdata$sample %in% 'M0','gene'], na.rm=T)
    #testdata$gene <- testdata$gene - M0mean
    
    #get stats from limma
    stat.limma <- rbind(human_WAT_M1vsM0[genename,],
                        human_WAT_M2vsM0[genename,],
                        human_WAT_M1vsM2[genename,])
    
    #get stats from ggpubr
    stat.test <- t_test(testdata, gene ~ sample)
    
    #replace stats by the ones from limma
    stat.test$p <- stat.limma$P.Value
    stat.test$p.adj <- stat.limma$adj.P.Val
    stat.test$p.adj.signif <- "ns"
    stat.test$p.adj.signif[stat.test$p.adj<0.05] <- "*"
    stat.test$p.adj.signif[stat.test$p.adj<0.01] <- "**"
    stat.test$p.adj.signif[stat.test$p.adj<0.001] <- "***"
    
    # Add position for adjusted p-values
    stat.test <- stat.test %>% add_xy_position(x = "sample")
    stat.test$y.position <- stat.test$y.position * 1.05
    
    #plot
    plot_FACS <- ggboxplot(testdata, x = "Treatment", y = "gene", fill="Markers") + 
      theme_bw() + theme + guides(fill=guide_legend(nrow=4)) +
      labs(title="Human adipose\ntissue macrophages",
           x=element_blank(),
           y=paste(genename, "log2(relative mRNA)")) +
      geom_boxplot(aes(x=Treatment, y=gene, fill=Markers)) +
      geom_dotplot(aes(x=Treatment, y=gene), fill="black",
                   binwidth = binwidth/60, alpha=0.9,
                   binaxis = "y", stackdir = "center", binpositions="all") +
      scale_fill_manual(values=c("#999999", "#CC79A7", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
      geom_point(size=0.5) +
      stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = T, size=7, bracket.size = 0.5) 
    plot_FACS
    return(plot_FACS)
    
  }
  
  output$FACSPlot <- renderPlot({
    FACSPlot_function()
  })
  
  output$human_WAT_means <- renderTable(align="c", spacing="xs", rownames=T, colnames=T,{
    genename <- toupper(input$inputGeneSymbol)
    data.frame(mean=as.character(round(as.numeric(human_WAT_means[genename,1:3]), 2)),
               Sd=as.character(round(as.numeric(human_WAT_means[genename,4:6]), 2)),
               n=as.character(human_WAT_means[genename,7:9]),
               row.names = c("All", "M1", "M2"))
  })
  
  
  ##################################################################################################################
  #Statistic tables
  output$human_stats_table <- renderDataTable(options=list(signif = 3),{
    datatable(human_stats_table, 
                  options = list(lengthMenu = c(5, 10, 50, 100), 
                                 pageLength = 5),
                  rownames= FALSE) %>% 
      formatSignif(columns = c(2:7), digits = 3)
  })
  
  output$mouse_stats_table <- renderDataTable({
    datatable(mouse_stats_table, options = list(lengthMenu = c(5, 10, 50, 100), 
                                                    pageLength = 5),
                  rownames= FALSE) %>% 
      formatSignif(columns = c(2:7), digits = 3)
  })
  
  output$human_WAT_stats_table <- renderDataTable({
    datatable(human_WAT_stats_table, options = list(lengthMenu = c(5, 10, 50, 100), 
                                                        pageLength = 5),
                  rownames= FALSE) %>% 
      formatSignif(columns = c(2:7), digits = 3)
  })
  
  
  ##################################################################################################################
  #datasets
  output$datasets <- renderDataTable({
    datatable(datasets, options = list(lengthMenu = c(10, 50, 100), 
                                           pageLength = 10),
                  rownames= FALSE)
  })
  
  
}

# Run the app ----
shinyApp(ui = ui, server = server)