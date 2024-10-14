#-----------------------------------------------------------------------
#
# Mouse Obesity
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

#load data
human_datamatrix <- readRDS('data/human_datamatrix.Rds')
human_metadata <- readRDS("data/human_metadata.Rds")
human_genelist <- rownames(human_datamatrix)
human_datasets <- human_metadata[,c("GEO", "platform")]
human_datasets$species <- "human"
human_datasets <- human_datasets[!duplicated(human_datasets$GEO),]

# function to format p values
p_value_formatter <- function(p) {
  sapply(p, function(x) {
    if (x < 0.001) {
      return("italic(p) < 0.001")
    } else {
      return(sprintf("italic(p) == %.3f", x))
    }
  })
}

# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Gene expression in skeletal muscle from humans with obesity and type 2 diabetes"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), 
                            "/ last update 2024-10-14")
                ),
                
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;;text-align:center",
                         em(h5("If you know other datasets that could be incorporated in this tool, please",
                               a("let me know!", href="mailto:nicolas.pillon@ki.se", 
                                 target="_blank", style="color:#5B768E"))),
                         tags$hr()
                ),
                
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         sidebarLayout(
                           sidebarPanel(width = 2, 
                                        selectizeInput("inputHumanGeneSymbol", "Gene Symbols:", choices=NULL, multiple=F, width=1000),
                                        checkboxGroupInput("diagnosis_diabetes", 
                                                           label = "Diabetes diagnosis", 
                                                           selected = c("Healthy", "Prediabetes", "T2D"),
                                                           choices = c("Healthy", "Prediabetes", "T2D")),
                                        em(h5("Prediabetes is defined as either impaired glucose tolerance measured during an OGTT or increased insulin
                                        resistance measured with euglycemic hyperinsulinemic clamps."))
                         ),
                         mainPanel(width = 9,
                                   fluidRow(style="color:black;background-color:white;",
                                            plotOutput("genePlotObesity", height="800px") %>% withSpinner(color="#5b768e")
                                            ),
                                   
                                   fluidRow(style="color:black;background-color:white;",
                                            h3("Datasets Included in the Analysis"),
                                            dataTableOutput("datasets")
                                            )
                )
                )
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputHumanGeneSymbol', 
                       choices=human_genelist, 
                       server=TRUE, 
                       selected=c("LEP"), 
                       options=NULL)
  
  output$genePlotObesity <- renderPlot({
    validate(need(input$inputHumanGeneSymbol, " "))
    genename <- c("AKR1C3")
    genename <- input$inputHumanGeneSymbol

    # Obesity
    plotdata <- data.frame(human_metadata, 
                           genedata = as.numeric(human_datamatrix[genename,]))
    plotdata$bmi_category <- factor(plotdata$bmi_category, 
                                    levels=c("Lean", "Overweight", "Obesity"))

    #filter according to selected categories
    plotdata <- plotdata[plotdata$diagnosis %in% input$diagnosis_diabetes 
                         ,]
    
    # label with n size
    plotdata$sex <- gsub("^male", paste0("Male, n = ", nrow(plotdata[plotdata$sex == "male",])), plotdata$sex)
    plotdata$sex <- gsub("^female", paste0("Female, n = ", nrow(plotdata[plotdata$sex == "female",])), plotdata$sex)
    
    cowplot::plot_grid(
      ggplot(plotdata, aes(x=bmi, y=genedata)) +  
        geom_point(aes(color = diagnosis, shape = diagnosis), size = 3, alpha = 0.25)  + 
        geom_smooth(method = "lm", color = "black", se = FALSE) +
        theme_bw(16) +  
        theme(legend.position = "none") +
        facet_wrap(.~sex, ncol = 1) +
        labs(x="BMI, kg/m2",
             y="mRNA expression, log2") +
        scale_shape_manual(values=rep(c(15,16,17), 20)) +
        scale_color_manual(values = c("darkgreen", "orange", "darkred")) +
        scale_y_continuous(expand = c(0,1)) +
        stat_cor(size = 4, 
                 vjust = -1, 
                 label.x = 20),
      
      ggplot(plotdata, aes(x=bmi_category, y=genedata)) +  
        geom_jitter(aes(color = diagnosis, shape = diagnosis), size = 2, alpha = 0.25, width = 0.1) +
        geom_boxplot(outlier.size = 0.1, fill = "gray80", alpha = 0.5)  + 
        theme_bw(16) + 
        theme(legend.position = "right") +
        facet_wrap(.~sex, ncol = 1) +
        labs(x=element_blank(),
             y="mRNA expression, log2") +
        scale_shape_manual(values=rep(c(15,16,17), 20)) +
        scale_color_manual(values = c("darkgreen", "orange", "darkred")) +
        scale_y_continuous(expand = c(0,1)) +
        stat_compare_means(aes(label = after_stat(p_value_formatter(..p..))), 
                           ref.group = "Lean",
                           parse = TRUE,
                           size = 4, 
                           vjust = -1),
      
      ncol = 2
    )
    
  })
  
  ##################################################################################################################
  #Dataset tables
  output$datasets <- renderDataTable(options=list(signif = 3),{
    DT::datatable(human_datasets,
                  options = list(lengthMenu = c(10, 50, 100),
                                 pageLength = 10),
                  rownames= FALSE)
  })
  
}


# Run the app ----
shinyApp(ui = ui, server = server)
