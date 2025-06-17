#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle cells
#
#--------------------------------------------------------------------------------------------------------

ui <- fluidPage(title="MyotubePalmitate",
                style="padding:0%",
                
                # CSS for style
                tags$head(includeCSS("../../www/style.css")),
                
                # main page
                navbarPage("MyotubePalmitate",
                           
                           # description of methods
                           tabPanel("Description",

                                    h3("Objective"),
                                    p("This integrated dataset provides a robust resource to study palmitate-induced transcriptomic changes in skeletal muscle cells."),
                                    
                                    h3("Methods"),
                                    p("Publicly available datasets from the Gene Expression Omnibus (GEO) were downloaded, including both microarray and RNA-seq platforms from human, mouse, and rat studies. Each dataset was processed using platform-appropriate normalization procedures, such as RMA for microarrays and VST transformation for RNA-seq count data. Lowly expressed genes were filtered based on manually defined expression thresholds, and only samples with paired control and palmitate treatments were retained."),
                                    p("Cross-species orthologs were annotated using biomaRt, allowing all datasets to be merged using human gene identifiers. Expression matrices were median-centered, quantile-normalized, and batch-corrected using study ID as the batch factor. Missing values were imputed using K-nearest neighbor imputation. Quality control was performed by inspecting known palmitate-responsive genes such as ANGPTL4 and PDK4, and by computing pairwise correlations of palmitate-induced gene expression changes across sample pairs. Outlier pairs, defined as those with average correlation in the lowest 5%, were excluded from the final analysis."),
                                    p("Wilcoxon ranked signed test comparing palmitate to control. The result table shows the p values as well as the adjusted p values corrected for multipple testing using Bonferroni correction."),
                                    
                                    h3("Datasets"),
                                    tags$p(
                                      tags$b("Are we missing a relevant study? Please "),
                                      a("let us know!", href = "mailto:nicolas.pillon@ki.se", target = "_blank")
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
                           
                           # Panel for plots
                           tabPanel(
                             title = "Data",
                             sidebarLayout(
                               sidebarPanel(width = 3,
                                            selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                                            sliderInput("concentration", tags$b("Concentation (µmol/L)"),
                                                        min = 100, max = 500, value = c(100,500), step = 100, sep = ""),
                                            sliderInput("duration", tags$b("Duration (hours)"),
                                                        min = 12, max = 96, value = c(12,96), step = 12, sep = ""),
                                            checkboxGroupInput("cell_type", 
                                                               label = "Cell type", 
                                                               selected = c("C2C12", "LHCN-M2", "primary"),
                                                               choices = c("C2C12", "LHCN-M2", "primary")),
                                            checkboxGroupInput("species", 
                                                               label = "Species", 
                                                               selected = c("human", "mouse", "rat"),
                                                               choices = c("human", "mouse", "rat")),
                                            actionButton("updatePlot", "Refresh plot", icon("refresh")),
                                            tags$br(),tags$br(),
                                            downloadButton("downloadGeneData", "Download Data")
                               ),
                               mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                                         plotOutput("geneBoxplot", height="500px") %>% withSpinner(color="#5B768E", type = 8),
                                         DT::dataTableOutput("statistics1", width = "90%") %>% withSpinner(color="#5B768E", type = 8),
                                         DT::dataTableOutput("statistics2", width = "90%") %>% withSpinner(color="#5B768E", type = 8))

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
)

                           
                