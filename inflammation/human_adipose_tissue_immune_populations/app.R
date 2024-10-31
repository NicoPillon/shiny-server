#-----------------------------------------------------------------------
#
# FACS of adipose tissue resident cells
#
#----------------------------------------------------------------------
# Load libraries
library(shinycssloaders)
library(stringr)
library(DT)
library(plyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(feather)
library(ggforce)
library(viridis)

#load data
datamatrix <- readRDS('data/matrix.Rds')
metadata <- readRDS("data/metadata.Rds")
genelist <- rownames(datamatrix)

#add group
metadata$group <- "Adipose cells"
metadata$group[metadata$CellType %in% c("Monocyte/Macrophage", "Leukocyte")] <- "Leukocytes"
metadata$group[metadata$CellType %in% c("Total T-cells", "CD4+ T-cells", "CD8+ T-cells")] <- "Lymphocytes"
metadata$group[metadata$CellType %in% c("M1 Macrophage", "M2 Macrophage", "CD14+ Myeloid")] <- "Macrophages"
metadata$group[metadata$CellType %in% c("SVF")] <- "SVF"

# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                
                # Google analytics
                tags$head(includeScript("google-analytics.html")),
                
                # Custom CSS to change checkbox tick color
                tags$style(HTML("
                  input[type='checkbox'] {
                    accent-color: #c93f1e; /* Change the checkbox tick color here */
                  }
                ")),
                
                # title ribbon
                fluidRow(style="color:white;background-color:#5B768E;padding:0% 1% 1% 1%;text-align:center; display:flex; justify-content:center; align-items:center;",
                         column(1, 
                                tags$a(href = "https://shiny.nicopillon.com", 
                                       icon("home", class = "fa-2x"), 
                                       style = "color:white; text-decoration:none;"
                                )
                         ),
                         column(10,
                                h3("Transcriptomic profile of cells from human adipose tissue"),
                                h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                           target="_blank", style="color:#D9DADB"), 
                                   "/ last update 2024-10-31")
                         ),
                         column(1,
                         )
                ),
                
                # main page
                fluidRow(style="color:black;background-color:white;padding:1% 8% 1% 8%;",
                         sidebarLayout(
                           sidebarPanel(width = 3,
                                        selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                                        actionButton("updatePlot", "Refresh plot", icon("refresh"))
                           ),
                           mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                                     plotOutput("heatmap", height="250px") %>% withSpinner(color="#5b768e", proxy.height=300),
                                     plotOutput("barPlot", height="600px") %>% withSpinner(color="#5b768e", proxy.height=300)

                           )
                         ),
                         
                ),
                
                tags$hr(),
                
                # Table with datasets
                fluidRow(style="color:black;background-color:white;padding:1% 8% 0% 8%;",
                         "Cells from digested human adipose tissue: stroma-vascular fraction and mature adipocytes. From",
                         a("GSE80654", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80654", 
                           target="_blank", style="color:#5B768E"), 
                         "(HTA-2.0) and", 
                         a("GSE100795", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100795", 
                           target="_blank", style="color:#5B768E"), 
                         "(RNAseq)"
                ),
                
                tags$br(),
                
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
                                tags$b("Disclaimer:"), tags$br(),
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
                       selected=c("CD3G", "CD3D", #CD3 
                                  "CD4", 
                                  "CD8A", #CD8 
                                  "CD14", 
                                  "CD34", 
                                  "PECAM1", #CD31 
                                  "PTPRC", #CD45 
                                  "MRC1", #CD206
                                  "LEP",
                                  "ITGAX"),
                       options=NULL)
  
  plotData <- eventReactive(input$updatePlot, {
    validate(need(input$inputGeneSymbol, " "))
    genename <- c("ADIPOQ", "MSR1", "VWF")
    genename <- input$inputGeneSymbol
    
    data <- data.frame()
    for (i in genename){
      asd <- datamatrix[i,]
      asd <- data.frame(metadata, as.numeric(asd), Gene=i)
      colnames(asd) <- c(colnames(metadata), "data", "Gene")
      data <- rbind(data, asd)
    }
    
    return(data)

  })

  
  output$heatmap <- renderPlot({
    data <- plotData()
    
    # summarize data
    mean_expression <- data %>%
      group_by(Gene, CellType) %>%
      summarize(mean_expression = median(data, na.rm = TRUE))
    
    # Apply min-max scaling by Gene
    mean_expression <- mean_expression %>%
      group_by(Gene) %>%
      mutate(scaled_expression = (mean_expression - min(mean_expression)) / 
               (max(mean_expression) - min(mean_expression))) %>%
      ungroup()
    
    # Determine the order of genes based on median expression for "Adipocyte"
    gene_order <- mean_expression %>%
      filter(CellType == "Adipocyte") %>%
      arrange(desc(mean_expression)) %>%
      pull(Gene)
    
    # Reorder Gene factor levels based on this order
    mean_expression$Gene <- factor(mean_expression$Gene, levels = gene_order)
    
    # Reverse the levels of CellType to invert y-axis
    unique(mean_expression$CellType)
    mean_expression$CellType <- factor(mean_expression$CellType, 
                                       levels = c("Total T-cells", "CD8+ T-cells", "CD4+ T-cells",
                                                  "Leukocyte",
                                                  "Monocyte/Macrophage", 
                                                  "CD14+ Myeloid",
                                                  "M2 Macrophage", "M1 Macrophage",
                                                  "SVF", "Progenitor", "Adipocyte"))
    
    # Create the heatmap
    
    ggplot(mean_expression, aes(x = Gene, y = CellType, fill = scaled_expression )) +
      geom_tile(color = "white") +
      scale_fill_viridis_c(option = "viridis") + # Use viridis color scale
      theme_minimal(16) +
      labs(x = element_blank(), y = element_blank(), fill = "Scaled Expression") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for readability
  })


  output$barPlot <- renderPlot({
    data <- plotData()
    
    data$CellType <- factor(data$CellType, 
                            levels = rev(c("Total T-cells", "CD8+ T-cells", "CD4+ T-cells",
                                       "Leukocyte",
                                       "Monocyte/Macrophage", 
                                       "CD14+ Myeloid",
                                       "M2 Macrophage", "M1 Macrophage",
                                       "SVF", "Progenitor", "Adipocyte")))
    
    gg <- ggplot(data, aes(x=CellType, y=data, fill=group)) +  
      geom_boxplot(size = 0.25, outlier.size = 0.2)  + theme_bw(16) +
      theme(legend.position = "bottom") +
      facet_wrap(~Gene, scales="free_y") +
      labs(x=element_blank(),
           y="mRNA expression, log2",
           fill = element_blank()) +
      scale_fill_viridis_d(option = "magma") + # Use viridis color scale
      theme(axis.text.x = element_text(color="black", angle=90, hjust=1, vjust=0.5))
    gg
    
  })
  
}




# Run the app ----
shinyApp(ui = ui, server = server)
