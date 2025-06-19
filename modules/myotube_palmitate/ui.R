#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle cells - UI
#
#--------------------------------------------------------------------------------------------------------

ui <- fluidPage(
  title = "MyotubePalmitate",         # Title of the browser tab
  style = "padding:0%",              # Remove default padding around the page
  
  #---------------------------------------------------------
  # Load custom CSS for consistent style across the app
  tags$head(includeCSS("../../www/style.css")),
  
  #---------------------------------------------------------
  # Navigation bar layout with multiple tabs
  navbarPage(
    # Custom title: logo image + text
    title = HTML('
      <div style="display: flex; align-items: center;margin: -10px;">
        <img src="../../www/img/snippet/myotube_palmitate.png" style="height: 40px; margin-right: 10px;">
        <span style="font-size: 20px; font-weight: bold; margin-right: 10px;">Palmitate</span>
      </div>
    '),
    
    #=====================================================================
    #============================ TAB 1 =================================
    #=====================================================================
    tabPanel("Analysis",  # Main tab for exploring data and results
             
             # Introductory text with link to the methods section
             p(HTML('Skeletal muscle cells response to palmitate exposure using integrated transcriptomic data. \n To understand how the data was generated, read the 
             <a href="#" onclick="$(\'.navbar-nav a:contains(\\\'Methods\\\')\').click()">methods</a>.')),
             
             tags$br(),
             
             # Layout: Sidebar on the left + Main content on the right
             sidebarLayout(
               
               #-----------------------------
               # Sidebar panel with input filters
               sidebarPanel(width = 3,
                            h4(tags$b("Select your criteria")),
                            
                            # Gene selector (multi-select)
                            selectizeInput("inputGeneSymbol", 
                                           "Genes of Interest", 
                                           selected = c("PDK4", "PXMP4", "LSM2", "ANGPTL4", "CPT1A", "ACAA2"),
                                           choices = gene_list$TARGET, 
                                           multiple = TRUE, 
                                           width = 600),
                            
                            # Sliders for palmitate concentration and exposure duration
                            sliderTextInput("concentration", "Palmitate Concentration (µmol/L)",
                                            choices = c("100", "200", "400", "500"),
                                            selected = c("100", "500"),
                                            grid = TRUE, force_edges = TRUE),
                            sliderTextInput("duration", "Palmitate Exposure (hours)",
                                            choices = c("12", "16", "18", "20", "24", "30", "36", "42", "48", "54", "96"),
                                            selected = c("12", "96"),
                                            grid = TRUE, force_edges = TRUE),
                            
                            # # Sliders for palmitate concentration and exposure duration
                            # sliderInput("concentration", "Palmitate Concentration (µmol/L)",
                            #             min = 100, max = 500, value = c(100, 500), step = 100, sep = ""),
                            # sliderInput("duration", "Palmitate Exposure (hours)",
                            #             min = 12, max = 96, value = c(12, 96), step = 12, sep = ""),
                            
                            # Checkbox filters for cell type and species
                             checkboxGroupInput("cell_type", 
                                               label = "Cell Type", 
                                               selected = c("human primary", "human LHCN-M2", "mouse C2C12", "rat primary"),
                                               choices = c("human primary", "human LHCN-M2", "mouse C2C12", "rat primary")),
                            
                            # Button to reset all filters
                            actionButton("resetInputs", "Reset Filters", icon = icon("undo")),

                            tags$hr(),
                            
                            # Download section
                            h4(tags$b("Download")),
                            downloadButton("downloadPlot", "Plot (.png)"), tags$br(), 
                            downloadButton("downloadData", "Data (.csv)"), tags$br(), 
                            downloadButton("downloadStats", "Statistics (.csv)")
               ),
               
               #-----------------------------
               # Main panel with output: plot + tables
               mainPanel(width = 9, style = "padding:0% 4% 0% 1%;",
                         
                         # Plot of gene expression (boxplot)
                         plotOutput("geneBoxplot", height = "550px") %>% withSpinner(color = "#5B768E", type = 8),
                         
                         #tags$hr(),
                         tags$em(textOutput("filterSummaryText")),
                         
                         tags$hr(),
                         
                         # Table 1: Differential expression results
                         DT::dataTableOutput("statistics1"),
                         tags$br(),
                         
                         # Table 2: Summary statistics (mean, SD, n)
                         DT::dataTableOutput("statistics2"),
                         tags$br()
               )
             )
    ),
    
    #=====================================================================
    #============================ TAB 2 =================================
    #=====================================================================
    tabPanel("Methods",  # Tab describing methodology
             
             style = "padding:0% 5% 1% 5%;",
             
             # Introduction
             h3("Why this app?"),
             p("This app offers an accessible resource to explore how saturated fats affect skeletal muscle cells. 
         Palmitate, the most abundant saturated fatty acid in Western diets, is closely linked to inflammation 
         and metabolic diseases. By examining its impact on the muscle cell transcriptome, this tool opens new 
         avenues for research, therapeutic discovery, and prevention strategies targeting metabolic disorders."),
             
             # Data processing and integration methodology
             h3("Methods"),
             p("Publicly available datasets were downloaded from the Gene Expression Omnibus (GEO), including both 
        microarray and RNA-seq platforms from human, mouse, and rat skeletal muscle cell models. Each dataset 
        was processed using platform-appropriate normalization methods — Robust Multi-array Average (RMA) for 
        microarrays and Variance Stabilizing Transformation (VST) for RNA-seq. Lowly expressed genes were filtered 
        out, and only samples with paired control and palmitate-treated conditions were retained."),
             p("Orthologous gene annotations were obtained using BiomaRt to enable integration across species, using human 
        gene identifiers as the reference. Expression matrices were centered, quantile-normalized, and batch-corrected 
        using study ID as a covariate. Missing values were imputed using k-nearest neighbor imputation. Quality control 
        included the inspection of known palmitate-responsive genes (e.g., ANGPTL4, PDK4) and the computation of pairwise 
        correlations across sample pairs. Sample pairs with poor correlation (lowest 5%) were flagged as outliers and removed."),
             
             # Statistical tests
             h3("Statistics"),
             p("Statistical analysis is based on the Wilcoxon signed-rank test comparing palmitate-treated samples to controls. 
        The results table presents both unadjusted and Bonferroni-adjusted p-values to correct for multiple testing across 
        all transcripts in the database, ensuring a highly conservative approach. Significance is indicated as follows: 
        * for FDR < 0.05, ** for FDR < 0.01, and *** for FDR < 0.001. 'ns' denotes non-significant results."),
             
             # Citation
             h3("Citation"),
             HTML('
        <div style="margin: 0rem;">
          <div>
            Pillon NJ, Sardón Puig L, Altıntaş A, Kamble PG, Casaní-Galdón S, Gabriel BM, Barrès R, Conesa A, Chibalin AV, Näslund E, Krook A, Zierath JR.
            <strong>Palmitate impairs circadian transcriptomics in muscle cells through histone modification of enhancers.</strong>
            <em>Life Sci Alliance.</em> 2022 Oct 27;6(1):e202201598.
            <a href="https://doi.org/10.26508/lsa.202201598" target="_blank">https://doi.org/10.26508/lsa.202201598</a>
          </div>
          <div style="margin-top: 1rem;">
            <script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>
            <span class="__dimensions_badge_embed__" data-id="pub.1152270218" data-legend="always"></span>
          </div>
        </div>
      '),
             
             # References table for included datasets
             h3("Datasets"),
             tags$p(
               tags$em("Have we missed a relevant study? Please ",
                       a("let us know!", href = "mailto:nicolas.pillon@ki.se", target = "_blank")
               )
             ),
             dataTableOutput("references")
    ),
    
    #=====================================================================
    #======================= IFRAME HEIGHT SYNC ==========================
    #=====================================================================
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
)
