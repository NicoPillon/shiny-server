#-----------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle response to hypoxia
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
AMPK_data <- readRDS("data/AMPK_data.Rds")


#get stats from limma
AMPK_stats_g3KO <- readRDS("data/AMPK_stats_g3KO.Rds")
AMPK_stats_g3TG <- readRDS("data/AMPK_stats_g3TG.Rds")
AMPK_stats_a1a2 <- readRDS("data/AMPK_stats_a1a2.Rds")
AMPK_stats_a1 <- readRDS("data/AMPK_stats_a1.Rds")
AMPK_stats_a2 <- readRDS("data/AMPK_stats_a2.Rds")
AMPK_stats_aicar <- readRDS("data/AMPK_stats_AICAR.Rds")


#--------------------------------------------------------------------------------------------------------
# List of genes
#--------------------------------------------------------------------------------------------------------
genelist <- unique(rownames(AMPK_data))


#function to gene name
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


#--------------------------------------------------------------------------------------------------------
# Sample
#--------------------------------------------------------------------------------------------------------
samples <- data.frame(sample=colnames(AMPK_data),
                      str_split_fixed(colnames(AMPK_data), "_|\\.", 3))[,1:3]
colnames(samples) <- c("sample", "GEO", "group")

samples$group <- gsub("AICAR.*", "AICAR", samples$group)

samples$model <- samples$group
samples$model[samples$GEO=="GSE4065"] <- "AMPKg3TG"
samples$model <- gsub("WT", "KO", samples$model)
samples$model <- gsub("CTRL|AICAR.*", "AICAR", samples$model)
samples$model <- factor(samples$model, levels=c("AICAR", "AMPKa1a2KO", "AMPKa1KO", "AMPKa2KO", "AMPKg3KO", "AMPKg3TG"))

samples$condition <- gsub("AMPK", "", samples$group)
samples$condition <- gsub("a1|a2|g3", "", samples$condition)
samples$condition <- gsub("AICAR.*", "AICAR", samples$condition)
samples$condition <- gsub("WT", "CTRL", samples$condition)
samples$condition <- factor(samples$condition, levels=c("CTRL", "KO", "TG", "AICAR"))



# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Transcriptomic response of skeletal muscle to AMPK modulation"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), "/ last update 2021-09-24")
                ),
                fluidRow(style="color:black;background-color:white;padding:1% 8% 0% 8%;",
                         "Skeletal muscle from mice with modulation of AMPK by knockout, mutation or AICAR."
                ),
                tags$hr(),
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                         plotOutput("ampkPlot", height="600px") %>% withSpinner(color="#5b768e")
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected=c("Bms1", "Mrps23", "Stk25", "Lyz2", "Gpx3", "Cox7b"), 
                       options=NULL)
  
  # AMPK plot
  output$ampkPlot <- renderPlot({
    validate(need(input$inputGeneSymbol, " "))
    genename <- c("Cd36", "Nos1", "Foxo1", "Ppargc1a","Srebf1",
                  "Bms1", "Mrps23", "Stk25", "Lyz2", "Gpx3", "Cox7b")
    genename <- input$inputGeneSymbol
    
    data <- data.frame()
    for (i in genename){
      asd <- AMPK_data[i,]
      asd <- data.frame(samples, as.numeric(asd), gene=i)
      colnames(asd) <- c(colnames(samples), "data", "gene")
      data <- rbind(data, asd)
    }
    
    #order data for stats calculations later
    data <- data[order(data$group),]
    
    #----------------------------------------------------------------------------------------------------------------
    # Statistics
    #----------------------------------------------------------------------------------------------------------------
    #get stats from limma
    stats_g3KO <- data.frame(gene=genename, AMPK_stats_g3KO[genename,], group1="AMPKg3KO", group2="AMPKg3WT")
    stats_g3TG <- data.frame(gene=genename, AMPK_stats_g3TG[genename,], group1="AMPKg3TG", group2="AMPKg3WT")
    stats_a1a2 <- data.frame(gene=genename, AMPK_stats_a1a2[genename,], group1="AMPKa1a2KO", group2="AMPKa1a2WT")
    stats_a1 <- data.frame(gene=genename, AMPK_stats_a1[genename,], group1="AMPKa1KO", group2="AMPKa1WT")
    stats_a2 <- data.frame(gene=genename, AMPK_stats_a2[genename,], group1="AMPKa2KO", group2="AMPKa2WT")
    stats_aicar <- data.frame(gene=genename, AMPK_stats_aicar[genename,], group1="AICAR", group2="CTRL")
    stat.limma <- rbind(stats_g3KO[,c(1,10,11,7,8)],
                        stats_g3TG[,c(1,10,11,7,8)],
                        stats_a1a2[,c(1,10,11,7,8)],
                        stats_a1[,c(1,10,11,7,8)],
                        stats_a2[,c(1,10,11,7,8)],
                        stats_aicar[,c(1,10,11,7,8)]
    )
    colnames(stat.limma)[1] <- "gene"
    
    #get stats from ggpubr
    stat.test <- compare_means(data=na.omit(data), data ~ group, 
                               group.by="gene")
    
    #merge with stats from limma
    stat.test <- full_join(stat.test, stat.limma)
    stat.test <- na.omit(stat.test)
    
    #replace stats by the ones from limma
    stat.test$p.adj.signif <- ""
    stat.test$p.adj.signif[stat.test$adj.P.Val<0.1] <- "*"
    stat.test$p.adj.signif[stat.test$adj.P.Val<0.05] <- "**"
    stat.test$p.adj.signif[stat.test$adj.P.Val<0.001] <- "***"
    
    # Add position for adjusted p-values
    stat.test <- stat.test %>% add_x_position(x = "model")
    stat.test$model <- stat.test$group1
    
    #manually get y position for each gene
    max.y <- data.frame(aggregate(data$data, by = list(paste(data$gene, data$model)), max, na.rm=T))
    max.y <- data.frame(t(data.frame(str_split(max.y$Group.1, " "))), max.y$x)
    colnames(max.y) <- c("gene", "model", "y.position")
    max.y$y.position <- max.y$y.position + 0.2
    stat.test <- full_join(stat.test, max.y)
    
    #barplot
    plotdata <- data[!grepl("CTRL", data$condition),]
    ggbarplot(plotdata, x="model", y="data", fill="model", add="mean_se",
              facet.by = "gene", scales="free") +
      geom_hline(yintercept = 0, linetype=1) +
      #geom_point(size=0.5) +
      labs(x=element_blank(),
           y="mRNA, log2(relative to CTRL)") +
      theme(plot.title  = element_text(face="bold", color="black", size=13, angle=0),
            axis.text.x = element_text(color="black", size=10, angle=90, hjust=1, vjust=0.5),
            axis.text.y = element_text(color="black", size=9, angle=0, vjust=0.3),
            axis.title  = element_text(face="bold", color="black", size=11, angle=0),
            strip.text.x = element_text(face="bold", color="black", size=13, angle=0),
            legend.text = element_text(color="black", size=12, angle=0),
            legend.title = element_blank(),
            legend.position="none",
            legend.key.size = unit(10, "pt")) +
      add_pvalue(stat.test, bracket.size = NA,
                 xmin = "group1", xmax = "group1",
                 label = "p.adj.signif",
                 y.position = "y.position", label.size=6) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))
    
  })
  
}


# Run the app ----
shinyApp(ui = ui, server = server)