#-----------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle and immune cells in response to exercse
#
#----------------------------------------------------------------------
# Load data and libraries
library(shinycssloaders)
library(tidyverse)
library(ggpubr)

#--------------------------------------------------------------------------------------------------------
# Load data
#--------------------------------------------------------------------------------------------------------
human_matrix <- cbind(
  readRDS("data/human_matrix_part1.Rds"),
  readRDS("data/human_matrix_part2.Rds"),
  readRDS("data/human_matrix_part3.Rds")
)
human_metadata <- readRDS("data/human_metadata.Rds")
mouse_matrix <- readRDS("data/mouse_matrix.Rds")
mouse_metadata <- readRDS("data/mouse_metadata.Rds")

scRNAseq <- readRDS("data/scRNAseq.Rds")

# clean metadata
human_metadata$Time <- gsub("H00", "POSTEX0", human_metadata$Time)
human_metadata$Time <- gsub("H01", "POSTEX1", human_metadata$Time)
human_metadata$Time <- gsub("H03", "POSTEX3", human_metadata$Time)

mouse_metadata$Time <- gsub("CTRL", "Control", mouse_metadata$Time)
mouse_metadata$Time <- gsub("H00", "Exercise", mouse_metadata$Time)

#--------------------------------------------------------------------------------------------------------
# List of genes
#--------------------------------------------------------------------------------------------------------
genelist <- c(rownames(human_matrix),
              rownames(mouse_matrix)) %>%
  str_to_title() %>%
  unique()



# Define UI ----
ui <- fluidPage(theme = "bootstrap.css",
                fluidRow(style="color:white;background-color:#5b768e;padding:0% 1% 1% 1%;text-align:center",
                         h3("Transcriptomic response of immune and muscle cells to metabolic stress"),
                         h5("By", a("Nicolas J. Pillon", href="https://staff.ki.se/people/nicolas-pillon", 
                                    target="_blank", style="color:#D9DADB"), "/ last update 2024-05-14")
                ),
                tags$hr(),
                fluidRow(style="color:black;background-color:white;padding:0% 8% 1% 8%;",
                         selectizeInput("inputGeneSymbol", "Gene Symbol:", choices=NULL, multiple=F, width=600)
                ),
                fluidRow(style="color:black;background-color:white;padding:2% 8% 1% 8%;",
                         plotOutput("PalPlot", height = 700) %>% withSpinner(color="#5b768e")
                )
)


# Define server logic ----
server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'inputGeneSymbol', 
                       choices=genelist, 
                       server=TRUE, 
                       selected=c("Ptgs2"), 
                       options=NULL)
  
  # AMPK plot
  plotData <- eventReactive(input$inputGeneSymbol, {
    validate(need(input$inputGeneSymbol, " "))
    gene_of_interest <- "Alox12"
    gene_of_interest <- input$inputGeneSymbol
    
   # subset gene of interest
    human_plotdata <- data.frame(
      human_metadata,
      data = as.numeric(human_matrix[toupper(gene_of_interest),])
    )
    
    mouse_plotdata <- data.frame(
      mouse_metadata,
      data = as.numeric(mouse_matrix[gene_of_interest,])
    )
    
    ###############################################################################
    ## Skeletal muscle - exercise
    ###############################################################################
    muscle_exercise <- human_plotdata[grepl("exercise", human_plotdata$Protocol),]
    muscle_exercise <- muscle_exercise[grepl("Muscle", muscle_exercise$Tissue),]
    
    muscle_exercise$Time <- factor(muscle_exercise$Time,
                                   levels = c("PRE", "POSTEX0", "POSTEX1", "POSTEX3"))
    muscle_exercise$dataset <- "MetaMEx"
    
    normdata_muscle <- data.frame()
    for( i in unique(muscle_exercise$GEO)){
      asd <- muscle_exercise[muscle_exercise$GEO == i,]
      pre <- mean(asd[asd$Time == "PRE", "data"])
      asd[,"data"] <- asd[,"data"] - pre
      normdata_muscle <- rbind(normdata_muscle, asd)
    }
    
    # remove 24h and rename and order
    normdata_muscle$Time <- gsub("POSTEX", "", normdata_muscle$Time)
    normdata_muscle$Time <- factor(normdata_muscle$Time, levels = c("PRE", "0", "1", "3"))

    # Plot
    fig_exercise_muscle <- ggline(normdata_muscle, x = "Time", y = "data", 
                                  color = "dataset", shape = "dataset", linetype = "dataset",
                                  add = "mean_se",
                                  size = 0.25) +
      theme_bw() +
      geom_hline(yintercept = 0, linetype = 1, linewidth = 0.2, color = "gray40") +
      labs(x = "Time after exercise (h)",
           y = paste(toupper(gene_of_interest), "mRNA\nlog2(relative to basal)"),
           subtitle = "Skeletal muscle tissue") +
      geom_hline(yintercept = 0, linewidth = 0.2, color = "gray50") +
      scale_color_manual(values = c("gray10", "gray20", "gray40")) +
      scale_linetype_manual(values = c(1,2,1)) +
      stat_compare_means(size = 4,
                         label.y = 1.75)

    
    ###############################################################################
    ## PBMC - exercise
    ###############################################################################
    dat_PBMC <- human_plotdata[grepl("exercise", human_plotdata$Protocol),]
    dat_PBMC <- dat_PBMC[grepl("PBMC", dat_PBMC$Tissue),]

    normdata_PBMC <- data.frame()
    for( i in unique(dat_PBMC$GEO)){
      asd <- dat_PBMC[dat_PBMC$GEO == i,]
      pre <- mean(asd[asd$Time == "PRE", "data"])
      asd[,"data"] <- asd[,"data"] - pre
      normdata_PBMC <- rbind(normdata_PBMC, asd)
    }
    
    # remove 24h and rename and order
    normdata_PBMC <- normdata_PBMC[!normdata_PBMC$Time %in% "POSTEX24",]
    normdata_PBMC$Time <- gsub("POSTEX", "", normdata_PBMC$Time)
    normdata_PBMC$Time <- factor(normdata_PBMC$Time, levels = c("PRE", "0", "1", "6"))

    # plot
    fig_exercise_pbmc <- ggline(normdata_PBMC, x = "Time", y = "data", add = "mean_se",
                                color = "GEO", shape = "GEO", linetype = "GEO",
                                size = 0.25) +
      theme_bw() +
      geom_hline(yintercept = 0, linetype = 1, linewidth = 0.2, color = "gray40") +
      labs(x = "Time after exercise (h)",
           y = paste(toupper(gene_of_interest), "mRNA\nlog2(relative to basal)"),
           subtitle = "PBMC") +
      scale_color_manual(values = c("gray10", "gray20", "gray40")) +
      scale_linetype_manual(values = c(1,2,1)) +
      stat_compare_means(size = 4,
                         label.y = 1.75)

    
    ###############################################################################
    ## Whole blood - exercise
    ###############################################################################
    dat_blood <- human_plotdata[grepl("exercise", human_plotdata$Protocol),]
    dat_blood <- dat_blood[grepl("Blood", dat_blood$Tissue),]
    
    normdata_blood <- data.frame()
    for( i in unique(dat_blood$GEO)){
      asd <- dat_blood[dat_blood$GEO == i,]
      pre <- mean(asd[asd$Time == "PRE", "data"])
      asd[,"data"] <- asd[,"data"] - pre
      normdata_blood <- rbind(normdata_blood, asd)
    }
    
    # remove 24h and rename
    normdata_blood <- normdata_blood[!normdata_blood$Time %in% "POSTEX24",]
    normdata_blood$Time <- gsub("POSTEX", "", normdata_blood$Time)

    # plots
    fig_exercise_blood <- ggline(normdata_blood, x = "Time", y = "data", add = "mean_se",
                                 color = "GEO", shape = "GEO", linetype = "GEO",
                                 size = 0.25) +
      theme_bw() +
      geom_hline(yintercept = 0, linetype = 1, linewidth = 0.2, color = "gray40") +
      labs(x = "Time after exercise (h)",
           y = paste(toupper(gene_of_interest), "mRNA\nlog2(relative to basal)"),
           subtitle = "Whole blood") +
      scale_color_manual(values = c("gray10", "gray20", "gray40")) +
      scale_linetype_manual(values = c(1,2,1)) +
      stat_compare_means(size = 4,
                         label.y = 1.75)

    ###############################################################################
    ## Adipose tissue - exercise
    ###############################################################################
    dat_WAT <- human_plotdata[grepl("exercise", human_plotdata$Protocol),]
    dat_WAT <- dat_WAT[grepl("Adipose", dat_WAT$Tissue),]

    normdata_WAT <- data.frame()
    for( i in unique(dat_WAT$GEO)){
      asd <- dat_WAT[dat_WAT$GEO == i,]
      pre <- mean(asd[asd$Time == "PRE", "data"], na.rm = T)
      asd[,"data"] <- asd[,"data"] - pre
      normdata_WAT <- rbind(normdata_WAT, asd)
    }
    
    # remove 24h and rename and order
    normdata_WAT$Time <- factor(normdata_WAT$Time,
                                levels = c("PRE", "POSTEX0", "POSTEX1", "POSTEX3"))
    normdata_WAT$Time <- gsub("POSTEX", "", normdata_WAT$Time)
    normdata_WAT$Time <- factor(normdata_WAT$Time, levels = c("PRE", "0", "1", "3"))
    
    fig_exercise_WAT <- ggline(normdata_WAT, x = "Time", y = "data", 
                               add = "mean_se", color = "GEO", shape = "GEO", linetype = "GEO",
                               size = 0.25) +
      theme_bw() +
      theme(legend.position = "bottom") +
      geom_hline(yintercept = 0, linetype = 1, linewidth = 0.2, color = "gray40") +
      labs(x = "Time after exercise (h)",
           y = paste(toupper(gene_of_interest), "mRNA\nlog2(relative to basal)"),
           subtitle = "Adipose tissue") +
      geom_hline(yintercept = 0, linewidth = 0.2, color = "gray50") +
      scale_color_manual(values = c("gray10", "gray20", "gray40")) +
      scale_linetype_manual(values = c(1,2,1)) +
      stat_compare_means(size = 4,
                         label.y = 1.75)

    ###############################################################################
    # Muscle sincle cell RNAseq
    ###############################################################################
    scRNAseq_selected <- scRNAseq[scRNAseq$Gene_name %in% toupper(gene_of_interest),]
    scRNAseq_selected$signif <- ifelse(scRNAseq_selected$pvals_adj < 0.05, TRUE, FALSE)
    fig_scRNAseq <- ggbarplot(scRNAseq_selected, x="cell.type", y="scores", fill="signif",
                              position = position_dodge(0.7),
                              size = 0.2) +
      theme_bw() + theme(legend.position = "right") +
      labs(x = element_blank(),
           y = paste(toupper(gene_of_interest), "mRNA\nscore(exercise induction)"),
           subtitle = "Skeletal muscle single cells",
           fill = "FDR < 0.05") +
      #geom_hline(yintercept = 0, linetype = 1, linewidth = 0.2, color = "gray40") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                         limits = c(0,NA)) +
      scale_fill_manual(values = c("gray50", "#2B8076"))
    fig_scRNAseq
    
    ###############################################################################
    # Mouse Tissues
    ###############################################################################
    fig_mouse <- ggboxplot(mouse_plotdata, x = "Time", y = "data", fill = "Time",
                           facet.by = "Tissue", ncol = 5,
                           size = 0.2,
                           outlier.size = 0.2, outlier.shape = NA) +
      geom_jitter(aes(shape = GEO), width = 0.07, size = 0.2) +
      labs(x = element_blank(),
           y = paste(toupper(gene_of_interest), "mRNA\nlog(a.u.)"),
           title = "Mouse tissues",
           fill = element_blank()) +
      theme_bw() +
      stat_compare_means(label = "p.format",
                         ref.group = "Control", 
                         size = 4,
                         vjust = -1.5) +
      scale_y_continuous(expand = c(0, 0.7)) +
      scale_fill_manual(values = c("gray90", "#2B8076")) +
      scale_shape_manual(values = c(22,23,24,25,16))
    
    ###############################################################################
    # Full figure
    ###############################################################################
    fig_human <- cowplot::plot_grid(
      fig_exercise_blood + 
        theme_bw(14) + theme(legend.position = "bottom",
                            legend.title = element_blank(),
                            legend.key.height = unit(10, "pt"),
                            plot.margin = margin(5.5, 5.5, 15, 5.5, "pt")) +
        guides(color=guide_legend(ncol = 1)) +
        expand_limits(y=c(-1.2, 1.5)),
      fig_exercise_pbmc + 
        theme_bw(14) + theme(legend.position = "bottom",
                            legend.title = element_blank(),
                            legend.key.height = unit(10, "pt"),
                            plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt")) +
        guides(color=guide_legend(ncol = 1)) +
        expand_limits(y=c(-1.2, 1.5)),
      fig_exercise_muscle + 
        theme_bw(14) + theme(legend.position = "bottom",
                            legend.title = element_blank(),
                            legend.key.height = unit(10, "pt"),
                            plot.margin = margin(5.5, 5.5, 26, 5.5, "pt")) +
        guides(color = guide_legend(ncol = 1)) +
        expand_limits(y=c(-1.2, 1.5)),
      fig_scRNAseq + 
        theme_bw(14) + theme(legend.position = "right",
                            axis.text.x = element_text(angle = 45, hjust=1),
                            plot.margin = margin(5.5, 5.5, 10, 5.5, "pt")),
      ncol = 4,
      rel_widths = c(1,1,1,1.25)
    )
    
    cowplot::plot_grid(
      fig_human,
      cowplot::plot_grid(
        fig_exercise_WAT +
          theme_bw(14) + theme(legend.position = "bottom",
                              legend.title = element_blank(),
                              legend.key.height = unit(10, "pt"),
                              plot.margin = margin(5.5, 5.5, 26, 5.5, "pt")) +
          guides(color=guide_legend(ncol = 1)) +
          expand_limits(y=c(-1.2, 1.5)),
        fig_mouse + 
          theme_bw(14) + theme(legend.position = "right",
                              legend.title = element_blank(),
                              legend.key.height = unit(10, "pt"),
                              plot.margin = margin(5.5, 5.5, 26, 5.5, "pt")),
        ncol = 2, rel_widths = c(1,3)
      ),
      ncol = 1
    )

    
    })
  
  
  
  output$PalPlot <- renderPlot({
    plotData()
  })

}


# Run the app ----
shinyApp(ui = ui, server = server)