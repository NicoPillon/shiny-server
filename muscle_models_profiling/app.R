#-----------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle cells
#
#----------------------------------------------------------------------
# Load data and libraries
library(shinycssloaders)
library(ggplot2)
theme <- theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
               legend.position = "none",
               legend.key.size = unit(30, "pt"))
library(stringr)
library(plyr)
library(DT)
library(dplyr)
library(DescTools)
library(ggpubr)
library(rstatix)
library(ggprism)
library(matrixStats)
library(pheatmap)

# Load data ----
data_all <- readRDS('data/Muscle_Models_Profiling_data.Rds')
data_all <- cbind(data_all[grepl('HumanCell', colnames(data_all))],
                  data_all[grepl('MouseC2C12', colnames(data_all))],
                  data_all[grepl('RatL6', colnames(data_all))],
                  data_all[grepl('HumanTissue', colnames(data_all))],
                  data_all[grepl('MouseTissue', colnames(data_all))],
                  data_all[grepl('RatTissue', colnames(data_all))])

#list of all samples (columns)
samples_list <- gsub("_G.*", "", colnames(data_all))
samples_list <- gsub(".*_", "", samples_list)

#find samples
samples_names <- data.frame(full.name=colnames(data_all))

samples_names <- gsub("_G.*", "", colnames(data_all))
samples_names <- gsub(".*_", "", samples_names)
samples_names <- data.frame(sample=unique(samples_names))

samples_names$model <- "Cell"
samples_names$model[samples_names$sample %in% c("HumanTissue", "RatTissue", "MouseTissue")] <- "Tissue"
samples_names$model.color <- c(rep("black", 3), rep("grey50", 3))

samples_names$species <- "Human"
samples_names$species[samples_names$sample %in% c("RatTissue", "RatL6")] <- "Rat"
samples_names$species[samples_names$sample %in% c("MouseTissue", "MouseC2C12")] <- "Mouse"
samples_names$species.color <- rep(c("#E69F00", "#56B4E9", "#CC79A7"), 2)

samples_names$fullname <- c("Human Primary Myotube", 
                            "Mouse C2C12 Myotube", 
                            "Rat L6 Myotube", 
                            "Human Muscle Tissue",
                            "Mouse Muscle Tissue", 
                            "Rat Muscle Tissue")
  
  
# List of genes
genelist <- unique(rownames(data_all))


#--------------------------------------------------------------------------------------------------------
# Shiny app
#--------------------------------------------------------------------------------------------------------
# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Gene expression in mouse, rat and human skeletal muscle tissue and cells"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), "/ last update 2022-03-08")
                ),
                
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;;text-align:center",
                         em(h5("If you use this tool, please cite: ", 
                            a("Abdelmoez et al. Am J Physiol Cell Physiol. 2020 Mar 1;318(3):C615-C626.",
                              href="https://doi.org/10.1152/ajpcell.00540.2019", target="_blank", style="color:#5B768E"))),
                         tags$hr()
                ),

                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                         actionButton("updatePlot", "Refresh plot", icon("refresh"))
                ),
                
                fluidRow(style="color:black;background-color:white;padding:2% 8% 1% 8%;",
                         column(3, align="left",
                                plotOutput("geneHeatmap", height="500px") %>% withSpinner(color="#5b768e")),
                         column(8, align="left",
                                plotOutput("geneBoxplot", height="500px") %>% withSpinner(color="#5b768e"))

                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected=c("FN1", "SPP1", "MYL4", "MYH7", "ATP2A1", "MYL3"), 
                       options=NULL)
  
  #-----------------------------------------------------------------
  # Boxplots
  plotDataBox <- eventReactive(input$updatePlot, {
    validate(need(input$inputGeneSymbol, " "))
    genename <- c("NR4A3", "IL6")
    genename <- toupper(input$inputGeneSymbol)
    
    #plot
    data <- data.frame()
    for( i in genename) { 
      asd <- data.frame(x=samples_list, 
                        y=as.numeric(data_all[i,]),
                        Gene=i, 
                        stringsAsFactors=FALSE) #empty dataframe to collect data
      data  <- rbind.data.frame(data, asd)              #bind the new gene data at the bottom of the previous one
    }
    
    data <- merge(data, samples_names, by=1)
      
    data$fullname <- factor(data$fullname, levels=samples_names$fullname)
    
    plotBoxplot <- ggplot(data) +  geom_boxplot(aes(x=fullname, y=y, fill=species)) + 
      facet_wrap(~Gene) +
      theme_bw(14) + theme +
      labs(x="",
           y="mRNA expression (log2)",
           title=element_blank()) +
      scale_y_continuous(breaks = round(seq(-4, 8, by=2),1)) +
      geom_hline(aes(yintercept=0), linetype="dashed", show.legend=F, color="gray60") +
      scale_fill_manual(values=samples_names$species.color)

    return(plotBoxplot)
    plotBoxplot
  })
  
    output$geneBoxplot <- renderPlot({
    plotDataBox()
  })
  
  
  #-----------------------------------------------------------------
  # Heatmap
  plotDataHeatmap <- eventReactive(input$updatePlot, {
    validate(need(input$inputGeneSymbol, " "))
    genename <- c("NR4A3", "IL6")
    genename <- toupper(input$inputGeneSymbol)

    #make matrix with found gene names
    df <- data_all[genename,]
    df[df=="NaN"] <- NA
    
    #Calculate the median value for each sample and organize in a matrix
    df_mean <- data.frame(SYMBOL=rownames(df))
    for (i in 1:length(samples_names$sample)){
      df_mean <- data.frame(df_mean,
                                  rowMedians(as.matrix(df[grepl(samples_names$sample[i], colnames(df))]), na.rm=T))
    }
    df_mean <- data.frame(df_mean, row.names=1)
    colnames(df_mean) <- samples_names$fullname
    
    return(df_mean)
  })

  
  output$geneHeatmap <- renderPlot({
   
    df_mean <- plotDataHeatmap()
    
    my_sample_col <- samples_names[,c(6,2,4)]
    my_sample_col <- data.frame(my_sample_col, row.names = 1)

    my_colour = list(
      model = setNames(unique(samples_names$model.color), unique(samples_names$model)),
      species = setNames(unique(samples_names$species.color), unique(samples_names$species))
    )
    
    plotHeatmap <- pheatmap(df_mean, scale = "row",
                            annotation_col = my_sample_col,
                            annotation_colors = my_colour,
                            display_numbers = F,
                            cluster_rows = F,
                            cluster_cols = F,
                            fontsize=11,
                            color = hcl.colors(50, "BluYl", rev = T),
                            angle_col=90)
    
    return(plotHeatmap$gtable)
  })
  
}


# Run the app ----
shinyApp(ui = ui, server = server)