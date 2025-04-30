#-----------------------------------------------------------------------
#
# Transcriptomic profile of human and mouse macrophages
#
#----------------------------------------------------------------------
# Load data and libraries
# Load data and libraries
library(feather)
library(shinycssloaders)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(cowplot)
library(DT)
library(matrixStats)
library(pheatmap)
library(arrow)

# special libraries for this app:
library(rstatix)
library(data.table)
library(zip)


references <- readRDS("data/datasets.Rds")


#Mouse data
mouse_data <- readRDS("data/mouse_data.Rds")
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


#gene list
genelist <- unique(c(rownames(mouse_data), rownames(human_data)))

#function to gene name
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Define UI ----
ui <- fluidPage(
  # CSS for style
  tags$head(includeCSS("../../www/style.css")),
  
  # title ribbon
  fluidRow(
    style = "color:black; padding:0% 2% 0% 2%; text-align:left;",
    h3("Transcriptomic signatures of macrophage polarization states"),
    h5("Last update 2022-04-07"),
    tags$hr()
  ),
  
  # Main layout with plot and controls
  fluidRow(style = "color:black;background-color:white;padding:0% 5% 1% 5%;",
           sidebarLayout(
             sidebarPanel(width = 3, 
                          selectizeInput("inputGeneSymbol", "Gene Symbol:", choices=NULL, multiple=F),   
                          tags$hr(),
                          tags$b("Statistics:"),
                          em("Expression differences were analyzed using limma with empirical Bayes moderation. Plots display Bonferroni-adjusted p-values. Significance is indicated as: * FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001."),
                          tags$hr(),
                          downloadButton("downloadGeneData", "Download Gene Data")
             ),
             
             mainPanel(width = 9,
                       fluidRow(
                         column(6, align="left",
                                plotOutput("mousePlot", height="500px") %>% withSpinner(color="#5b768e")),
                         column(6, align="left",
                                plotOutput("humanPlot", height="500px") %>% withSpinner(color="#5b768e"))
                       )
             )
           )
  ),
  
  # Description of methods
  fluidRow(
    style = "color:black; background-color:white; padding:0% 2% 0% 2%;",
    tags$hr(),
    h3("Methods"),
    p("Transcriptomic datasets of human primary blood–derived macrophages were identified through the GEO repository. These included macrophages polarized in vitro to M1 (stimulated with lipopolysaccharide, interferon-γ, or both) or M2 (stimulated with IL-4, IL-13, or both) phenotypes. Eight RNA-seq and thirteen microarray datasets were selected, totaling 57 M0 (unstimulated), 92 M1, and 87 M2 samples."),
    p("Gene expression data were downloaded, merged into a unified matrix, and filtered to exclude low-informative features such as microRNAs, long non-coding RNAs, ribosomal proteins, and hypothetical transcripts. Samples were annotated based on phenotype (M0, M1, M2), treatment condition, and GEO dataset identifier. Genes not detected in at least 20 samples per condition were excluded. Batch effects across studies were corrected using the `removeBatchEffect` function from the limma package."),
    p("Differential gene expression between groups (M1 vs. M0, M2 vs. M0, and M1 vs. M2) was assessed using empirical Bayes statistics implemented in the limma package. P-values were adjusted for multiple testing using the Bonferroni method."),
    p("Genes were considered specific to the M1 or M2 polarization states if they met the following criteria: FDR < 10⁻³ and log₂(fold change) > 0 compared to M0, and FDR < 10⁻³ and log₂(fold change) > 2 when comparing M1 and M2. The resulting gene signature lists were used to test for enrichment in external RNA sequencing datasets using Fisher’s exact test."),
    downloadButton("downloadDatasetHuman", "Download the human signature data and statistics"),
    downloadButton("downloadDatasetMouse", "Download the mouse signature data and statistics")
  ),
  
  # Table with references
  fluidRow(
    style = "color:black;background-color:white;padding:0% 2% 0% 2%;",
    tags$hr(),
    h3("Datasets Included in the Analysis"),
    dataTableOutput("references"),
    tags$p(
      tags$b("Are we missing a relevant study? Please "),
      a("let us know!", href = "mailto:nicolas.pillon@ki.se", target = "_blank")
    )
  ),
  
  # statistics
  fluidRow(style="color:black;background-color:white;padding:0% 2% 1% 2%;",
           h3("Statistics"),
           h4("Mouse primary macrophages"),
           dataTableOutput("mouse_stats_table"),
           h4("Human primary macrophages"),
           dataTableOutput("human_stats_table")
  ),
  
  # Citation
  fluidRow(
    style = "color:black;background-color:white;padding:0% 2% 2% 2%;",
    tags$hr(),
    h3("Citation"),
    p(
      "Pillon NJ, Smith JAB, Alm PS, Chibalin AV, Alhusen J, Arner E, Carninci P, Fritz T, Otten J, Olsson T, van Doorslaer de Ten Ryen S, Deldicque L, Caidahl K, Wallberg-Henriksson H, Krook A, Zierath JR. ",
      a("Distinctive exercise-induced inflammatory response and exerkine induction in skeletal muscle of people with type 2 diabetes.",
        href = "https://www.science.org/doi/10.1126/sciadv.abo3192", 
        target = "_blank", style = "color:#5B768E"),
      "Sci Adv. 2022 Sep 9;8(36):eabo3192."
    )
  ),
  
  # Code to send height to resizing iframe
  tags$head(
    tags$script(HTML("
    Shiny.addCustomMessageHandler('resizeFrame', function(message) {
      const height = document.documentElement.scrollHeight;
      parent.postMessage({ frameHeight: height }, '*');
    });

    window.addEventListener('resize', function() {
      const height = document.documentElement.scrollHeight;
      parent.postMessage({ frameHeight: height }, '*');
    });
  "))
  )
)


# Define server logic ----
server <- function(input, output, session) {
  
  # Code to send the height of the app to adjust iframe  
  session$onFlushed(function() {
    session$sendCustomMessage("resizeFrame", list())
  }, once = FALSE)
  
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
      theme_bw(16) + 
      theme(legend.position = "bottom",
            legend.title = element_blank()) +
      guides(fill=guide_legend(nrow=4)) +
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
      theme_bw(16) + 
      theme(legend.position = "bottom",
            legend.title = element_blank()) +
      guides(fill=guide_legend(nrow=4)) +
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
  
  
  
  #-----------------------------------------------------------------
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
  
  
  #-----------------------------------------------------------------
  # Dataset tables
  output$references <- renderDataTable({
    datatable(references, options = list(lengthMenu = c(10, 50, 100), 
                                         pageLength = 10),
              rownames= FALSE)
  })
  
  #-----------------------------------------------------------------
  # Download button - signature human
  # Download button - signature human
  output$downloadDatasetHuman <- downloadHandler(
    filename = function() {
      paste0("MuscleOmics.org_MacrophagePolarization_HumanSignature.zip")
    },
    content = function(file) {
      withProgress(message = "Preparing download...", value = 0, {
        tmpdir <- tempdir()
        paths <- c()
        
        # Step 1: Write expression matrix
        incProgress(0.1, detail = "Writing expression matrix...")
        human_data_file <- file.path(tmpdir, "data.csv")
        fwrite(cbind(SYMBOL = rownames(human_data), human_data), human_data_file)
        paths <- c(paths, human_data_file)
        
        # Step 2: Write limma results
        incProgress(0.3, detail = "Writing differential expression results...")
        file_M1vsM0 <- file.path(tmpdir, "limma_M1_vs_M0.csv")
        file_M2vsM0 <- file.path(tmpdir, "limma_M2_vs_M0.csv")
        file_M1vsM2 <- file.path(tmpdir, "limma_M1_vs_M2.csv")
        
        fwrite(human_M1vsM0, file_M1vsM0)
        incProgress(0.45)
        fwrite(human_M2vsM0, file_M2vsM0)
        incProgress(0.6)
        fwrite(human_M1vsM2, file_M1vsM2)
        incProgress(0.75)
        
        paths <- c(paths, file_M1vsM0, file_M2vsM0, file_M1vsM2)
        
        # Step 3: Zip
        incProgress(0.85, detail = "Compressing files...")
        zip::zipr(zipfile = file, files = paths, root = tmpdir)
        
        # Final step
        incProgress(1, detail = "Done.")
      })
    }
  )
  
  
  # Download button - signature mouse
  output$downloadDatasetMouse <- downloadHandler(
    filename = function() {
      paste0("MuscleOmics.org_MacrophagePolarization_MouseSignature.zip")
    },
    content = function(file) {
      withProgress(message = "Preparing download...", value = 0, {
        tmpdir <- tempdir()
        paths <- c()
        
        # Step 1: Write expression matrix
        incProgress(0.1, detail = "Writing expression matrix...")
        mouse_data_file <- file.path(tmpdir, "data.csv")
        fwrite(cbind(SYMBOL = rownames(mouse_data), mouse_data), mouse_data_file)
        paths <- c(paths, mouse_data_file)
        
        # Step 2: Write limma results
        incProgress(0.3, detail = "Writing differential expression results...")
        file_M1vsM0 <- file.path(tmpdir, "limma_M1_vs_M0.csv")
        file_M2vsM0 <- file.path(tmpdir, "limma_M2_vs_M0.csv")
        file_M1vsM2 <- file.path(tmpdir, "limma_M1_vs_M2.csv")
        
        fwrite(mouse_M1vsM0, file_M1vsM0)
        incProgress(0.45)
        fwrite(mouse_M2vsM0, file_M2vsM0)
        incProgress(0.6)
        fwrite(mouse_M1vsM2, file_M1vsM2)
        incProgress(0.75)
        
        paths <- c(paths, file_M1vsM0, file_M2vsM0, file_M1vsM2)
        
        # Step 3: Zip
        incProgress(0.85, detail = "Compressing files...")
        zip::zipr(zipfile = file, files = paths, root = tmpdir)
        
        # Final step
        incProgress(1, detail = "Done.")
      })
    }
  )
  
  
  
  
  # Download button - selected gene human and mouse
  output$downloadGeneData <- downloadHandler(
    filename = function() {
      paste0("MuscleOmics.org_MacrophagePolarization_", input$inputGeneSymbol, ".csv")
    },
    content = function(file) {
      genename <- input$inputGeneSymbol
      genename_mouse <- firstup(tolower(genename))
      genename_human <- toupper(genename)
      
      df_list <- list()
      
      # Human expression data
      if (genename_human %in% rownames(human_data)) {
        df_human <- data.frame(
          SampleID = colnames(human_data),
          Expression = as.numeric(human_data[genename_human, , drop = FALSE]),
          Species = "human",
          SYMBOL = genename_human
        )
        df_list[["human"]] <- df_human
      }
      
      # Mouse expression data
      if (genename_mouse %in% rownames(mouse_data)) {
        df_mouse <- data.frame(
          SampleID = colnames(mouse_data),
          Expression = as.numeric(mouse_data[genename_mouse, , drop = FALSE]),
          Species = "mouse",
          SYMBOL = genename_mouse
        )
        df_list[["mouse"]] <- df_mouse
      }
      
      # Combine and annotate
      df <- bind_rows(df_list)
      df <- df %>%
        mutate(
          Polarization = gsub("_.*", "", SampleID),
          Treatment = gsub("^[^_]+_", "", SampleID),
          Treatment = gsub("_GSE.*", "", Treatment),
          GEO = gsub(".*_GSE", "GSE", SampleID)
        ) %>%
        select(SampleID, GEO, Polarization, Treatment, Species, SYMBOL, Expression)
      
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  
  
}

# Run the app ----
shinyApp(ui = ui, server = server)