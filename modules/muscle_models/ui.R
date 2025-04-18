#--------------------------------------------------------------------------------------------------------
#
# Transcriptomic profile of skeletal muscle cells
#
#--------------------------------------------------------------------------------------------------------

ui <- fluidPage(
    
  # CSS for style
  tags$head(includeCSS("../../www/style.css")),
  
  # title ribbon
  fluidRow(
    style = "color:black; padding:0% 2% 0% 2%; text-align:left;",
    h3("Transcriptomic profiles of mouse C2C12, rat L6 and human skeletal myotubes"),
    h5("Last update 2025-04-11"),
    tags$hr()
  ),
  
  # main page
  fluidRow(style="color:black;background-color:white;padding:0% 5% 1% 5%;",
           sidebarLayout(
             sidebarPanel(width = 3, 
                          selectizeInput("inputGeneSymbol", "Gene Symbols:", choices=NULL, multiple=T, width=600),
                          actionButton("updatePlot", "Refresh plot", icon("refresh")),
                          tags$hr(),
                          downloadButton("downloadGeneData", "Download Data")
             ),
             mainPanel(width = 9, style="padding:0% 4% 1% 4%;",
                       fluidRow(style="color:black;background-color:white;",
                                column(5, align="left",
                                       plotOutput("geneHeatmap", height="600px") %>% withSpinner(color="#5B768E", type = 8)),
                                column(7, align="left",
                                       plotOutput("geneBoxplot", height="700px") %>% withSpinner(color="#5B768E", type = 8))
                       )
             )
             
           )
  ),
  
  # description of methods
  fluidRow(
    style = "color:black; background-color:white; padding:0% 2% 0% 2%;",
    tags$hr(),
    h3("Methods"),
    p("To compare gene expression across human, mouse, and rat skeletal muscle and myotube models, we performed a systematic integration and normalization of publicly available microarray datasets from multiple Affymetrix platforms. Raw CEL files were downloaded from the Gene Expression Omnibus (GEO) and organized by platform (e.g., GPL81, GPL570, GPL6244), with metadata indicating species and sample type (e.g., primary human muscle, L6 myotubes, C2C12 myotubes)."),
    
    p("For each platform, raw intensity values were background corrected and normalized using the Robust Multi-array Average (RMA) algorithm via the oligo R package. Probes were then mapped to Ensembl gene identifiers using platform-specific annotation packages (e.g., pd.hg.u133.plus.2) and supplemented with annotations retrieved from Ensembl via the biomaRt R package. A master annotation table was constructed for all platforms to standardize probe-to-gene mappings across arrays. Probe-level expression values were aggregated to the Ensembl gene level by computing the mean expression per gene."),
    
    p("To enable cross-species comparisons, we identified one-to-one orthologs between human, mouse, and rat using Ensembl. Gene-level ortholog annotations were retrieved using biomaRt, with priority given to Ensembl gene IDs. When ortholog mappings were unavailable at the Ensembl level, gene symbols were used to complete the mapping. Genes corresponding to pseudogenes, non-coding RNAs, ribosomal proteins, and poorly annotated loci (e.g., “LOC”, “RIK”, “MIR”) were excluded from the analysis. The final ortholog table included only protein-coding genes with reliable annotation in at least two species."),
    
    p("Expression data from all platforms were merged by gene symbol using the curated ortholog annotation table. To minimize technical variability, we excluded genes with more than 50% missing data across all samples. The resulting expression matrix was normalized across arrays using quantile normalization (normalizeBetweenArrays function from the limma package). No batch correction was applied across platforms to avoid artificially removing genuine biological differences between species, acknowledging the inherent trade-off between batch effect mitigation and preservation of species-specific signal.")
  ),
  
  # Table with datasets
  fluidRow(style="color:black;background-color:white;padding:0% 2% 0% 2%;",
           tags$hr(),
           h3("Datasets Included in the Analysis"),
           dataTableOutput("references"),
           tags$p(
             tags$b("Are we missing a relevant study? Please "),
             a("let us know!", href = "mailto:nicolas.pillon@ki.se", target = "_blank")
           )    
  ),
  
  # Citation
  fluidRow(style="color:black;background-color:white;padding:0% 2% 2% 2%;",
           tags$hr(),
           h3("Citation"),
           "Abdelmoez AM, Sardón Puig L, Smith JAB, Gabriel BM, Savikj M, Dollet L, Chibalin AV, Krook A, Zierath JR, Pillon NJ.",
           a("Comparative profiling of skeletal muscle models reveals heterogeneity of transcriptome and metabolism.",
             href="https://doi.org/10.1152/ajpcell.00540.2019", target="_blank", style="color:#5B768E"),
           "Am J Physiol Cell Physiol. 2020 Mar 1;318(3):C615-C626."
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
