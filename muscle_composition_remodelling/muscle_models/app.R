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
                
                # Custom CSS to change checkbox tick color
                tags$style(HTML("
                  input[type='checkbox'] {
                    accent-color: #5B768E; /* Change the checkbox tick color here */
                  }
                ")),
                
                # Title ribbon
                fluidRow(style="color:white;background-color:#5B768E;padding:0% 1% 1% 1%;text-align:center",
                         h3("Gene expression in mouse, rat and human skeletal muscle tissue and cells"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), 
                            "/ last update 2024-03-07")
                ),
                
                # main page
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         sidebarLayout(
                           sidebarPanel(width = 3, 
                                        selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                                        actionButton("updatePlot", "Refresh plot", icon("refresh")),
                                        tags$hr(),
                                        em("If you use this tool, please cite: ", 
                                              a("Abdelmoez et al. Am J Physiol Cell Physiol. 2020 Mar 1;318(3):C615-C626.",
                                                href="https://doi.org/10.1152/ajpcell.00540.2019", target="_blank", style="color:#5B768E"))
                           ),
                           mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                                     fluidRow(style="color:black;background-color:white;",
                                              column(4, align="left",
                                                     plotOutput("geneHeatmap", height="700px") %>% withSpinner(color="#5B768E")),
                                              column(8, align="left",
                                                     plotOutput("geneBoxplot", height="700px") %>% withSpinner(color="#5B768E"))
                                     )
                           )
                           
                         )
                ),
                
                tags$hr(),
                
                # Table with datasets
                fluidRow(style="color:black;background-color:white;padding:0% 2% 1% 2%;",
                         h3("Datasets Included in the Analysis"),
                         # dataTableOutput("datasets")
                ),
                
                # Author section at the bottom
                fluidRow(style="color:white;background-color:#5B768E;padding:2% 1% 2% 1%;display: flex; align-items: top; ",
                         # column(2, align="right", 
                         #        tags$img(src = "https://ki.se/profile-image/nicpil", height = "120px", width = "120px")  # Insert image here
                         # ),
                         column(4, align="left", 
                                tags$b("About the author:"), tags$br(),
                                "Nicolas J. Pillon, PhD", tags$br(),
                                "Associate Professor, Karolinska Institutet", tags$br(),
                                icon("globe"), a("/inflammation-and-metabolism", href="https://ki.se/en/research/research-areas-centres-and-networks/research-groups/inflammation-and-metabolism-nicolas-pillons-research-group",
                                                 target="_blank", style="color:white"), tags$br(),
                                icon("linkedin"), a("/nicopillon", href="https://www.linkedin.com/in/nicopillon/",
                                                    target="_blank", style="color:white"), tags$br(),
                                tags$br(),
                                "Feel free to write to me with feedback or questions:", tags$br(),
                                icon("envelope"), a("nicolas.pillon@ki.se", href="mailto:nicolas.pillon@ki.se",
                                                    target="_blank", style="color:white"), tags$br(),
                                
                         ),
                         column(4, align="center",
                                #tags$b("© 2024 Nicolas Pillon"), tags$br(),
                                
                         ),
                         column(4, align="right",
                                tags$b("Disclaimer"), tags$br(),
                                em("The authors disclaim any responsibility for the use or interpretation of the data 
                                   presented in this application. Users are solely responsible for ensuring the appropriate 
                                   use of any data they choose to re-use."), tags$br(),
                                tags$br(),
                                tags$b("© 2024 Nicolas Pillon"),
                         ),
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