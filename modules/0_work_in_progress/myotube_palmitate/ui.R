#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle cells
#
#--------------------------------------------------------------------------------------------------------

ui <- fluidPage(title="MyotubePalmitate",
                style="padding:0%",
                
                # CSS for style
                #tags$head(includeCSS("../../www/style.css")),
                
                # main page
                navbarPage(
                  title = HTML('
  <div style="display: flex; align-items: center;margin: -10px;">
    <img src="../../../www/img/snippet/myotube_palmitate.png" style="height: 40px; margin-right: 10px;">
    <span style="font-size: 20px; font-weight: bold; margin-right: 10px;">Palmitate</span>
  </div>
'),
                  # Panel for plots
                  tabPanel("Explore Data",

                           p(HTML('Explore how skeletal muscle cells respond to palmitate exposure using integrated transcriptomic data. 
       <a href="#" onclick="$(\'.navbar-nav a:contains(\\\'Read Me\\\')\').click()">Click here</a> to learn more about the methods and statistical analyses.')),
                           
                           tags$br(),
                           sidebarLayout(
                             sidebarPanel(width = 3,
                                          selectizeInput("inputGeneSymbol", 
                                                         "Select your genes of interest:", 
                                                         choices=NULL, multiple=T, width=600),
                                          sliderInput("concentration", "Palmitate Concentation (µmol/L)",
                                                      min = 100, max = 500, value = c(100,500), step = 100, sep = ""),
                                          sliderInput("duration", "Palmitate Exposure (hours)",
                                                      min = 12, max = 96, value = c(12,96), step = 12, sep = ""),
                                          checkboxGroupInput("cell_type", 
                                                             label = "Cell type", 
                                                             selected = c("C2C12", "LHCN-M2", "primary"),
                                                             choices = c("C2C12", "LHCN-M2", "primary")),
                                          checkboxGroupInput("species", 
                                                             label = "Species", 
                                                             selected = c("human", "mouse", "rat"),
                                                             choices = c("human", "mouse", "rat")),
                                          tags$b("Download your results:"), tags$br(),
                                          downloadButton("downloadPlot", "Plot (.png)"),
                                          downloadButton("downloadData", "Data (.csv)"),
                                          downloadButton("downloadStats", "Statistics (.csv)"),
                                          tags$hr(),
                                          actionButton("resetInputs", "Reset Filters", icon = icon("undo"))
                             ),
                             mainPanel(width = 9, style="padding:0% 4% 0% 1%;",
                                       plotOutput("geneBoxplot", height="550px") %>% withSpinner(color="#5B768E", type = 8),
                                       tags$hr(),
                                       DT::dataTableOutput("statistics1"),
                                       tags$br(),
                                       DT::dataTableOutput("statistics2"),
                                       tags$br()
                             )
                           )
                  ),
                  
                  
                  # description of methods
                  tabPanel("Read Me",
                           style="padding:0% 5% 1% 5%;",
                           
                           h3("Why this app?"),
                           p("This app offers a powerful and accessible resource to explore how saturated fats affect skeletal muscle cells. Palmitate, the most abundant saturated fatty acid in Western diets, is closely linked to inflammation and metabolic diseases. By examining its impact on the muscle cell transcriptome, this tool opens new avenues for research, therapeutic discovery, and prevention strategies targeting metabolic disorders."),
                           
                           h3("Methods"),
                           p("Publicly available datasets were downloaded from the Gene Expression Omnibus (GEO), including both microarray and RNA-seq platforms from human, mouse, and rat skeletal muscle cell models. Each dataset was processed using platform-appropriate normalization methods — Robust Multi-array Average (RMA) for microarrays and Variance Stabilizing Transformation (VST) for RNA-seq. Lowly expressed genes were filtered out, and only samples with paired control and palmitate-treated conditions were retained."),
                           p("Orthologous gene annotations were obtained using BiomaRt to enable integration across species, using human gene identifiers as reference. Expression matrices were centered, quantile-normalized, and batch-corrected using study ID as a covariate. Missing values were imputed using k-nearest neighbor imputation. Quality control included the inspection of known palmitate-responsive genes (e.g., ANGPTL4, PDK4) and the computation of pairwise correlations across sample pairs. Sample pairs with poor correlation (lowest 5%) were flagged as outliers and removed."),
                           
                           h3("Statistics"),
                           p("Statistical analysis is based on Wilcoxon signed-rank tests comparing palmitate-treated samples to controls. The results table presents both unadjusted p-values and Bonferroni-adjusted p-values to correct for multiple testing across all transcripts in the database, ensuring a highly conservative approach. Significance is indicated as follows: * for FDR < 0.05, ** for FDR < 0.01, and *** for FDR < 0.001. 'ns' denotes non-significant results."),
                           
                           h3("Datasets"),
                           tags$p(
                             tags$em("Have we missed a relevant study? Please ",
                                     a("let us know!", href = "mailto:nicolas.pillon@ki.se", target = "_blank")
                             )
                           ),    
                           dataTableOutput("references"),
                           
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
')
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
)


