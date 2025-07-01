########################################################################################
tabMethods <- tabPanel("Methods",
                       
                       style = "padding:0% 5% 1% 5%;",
                       
                       h3("Why this app?"),
                       p("This app was developed to provide an open-access, user-friendly interface to explore transcriptomic responses to exercise and inactivity across human and mouse studies. By harmonizing data from diverse platforms (microarrays and RNA-seq) and re-analyzing raw data when possible, this resource helps researchers identify gene expression patterns and compare responses across conditions, interventions, and populations. It enables dynamic filtering, visual exploration, and download of curated meta-analyses to support hypothesis generation and reproducibility in exercise biology."),
                       
                       #--------------------------------------------------------------------
                       h3("Methods"),
                       p("For data processing and statistical analysis, Affymetrix arrays were normalized using the Robust Multi-array Average (RMA) method via the ",
                         code("oligo"), " package, while other microarray platforms were normalized using quantile normalization through the ",
                         code("limma"), " package. RNA-seq data underwent variance-stabilizing transformation (VST) using ",
                         code("DESeq2"), ". Meta-analyses were performed using restricted maximum likelihood estimation with the ",
                         code("metafor"), " package. Statistical comparisons incorporated sample-size weighting, moderated t-tests using empirical Bayes statistics, and adjustments for multiple testing via Bonferroni correction. Studies were re-processed from raw files whenever possible."), 
                       
                       p("Older studies or custom arrays often include fewer probes, limiting the number of detectable genes. More recent RNA-seq datasets have greater depth and can detect non-coding RNAs not present in gene arrays. In some cases, limited metadata led to unpaired statistics and smaller subgroup sizes, reducing statistical power compared to original pooled analyses."),
                       
                       p("Duplicates, re-analyses, studies lacking controls or baseline, or confounded by non-exercise interventions (e.g., diet) were excluded. Some irrelevant subgroups (e.g., reloading, electrical stimulation) were also removed."),
                       
                       tags$ul(style = "margin-left: 2em;",
                               tags$li("GSE156247 was excluded because it contains the exact same data as GSE53598."),
                               tags$li("GSE163434 was excluded because it is a re-analysis of GSE157585."),
                               tags$li("GSE230002 and GSE137631 were not included because the intervention was a combination of exercise with weight loss (diet/surgery)."),
                               tags$li("GSE44818, GSE28998, and GSE24235 were excluded because they did not include pre-exercise or baseline biopsies."),
                               tags$li("GSE165630: This study compared trained athletes to untrained individuals, but the data clustered with acute exercise studies, suggesting that the biopsies were not taken long enough after the last exercise session."),
                               tags$li("GSE1718: This study compared individuals with or without changes in insulin sensitivity following a training intervention. The only data available was already processed in a way that made it impossible to separate conditions. In addition, the samples are a mix of male and female biopsies.")
                       ),
                       
                       p("A forest plot shows the results from multiple studies. Each study is listed with its log2 fold-change, false discovery rate (FDR), and sample size. The plot displays each result as a square (effect size) with horizontal lines representing the confidence intervals. The diamond at the bottom represents the combined meta-analytic estimate."),
                       
                       #--------------------------------------------------------------------
                       h3("Version History"),
                       tags$ul(style = "margin-left: 2em;",
                               tags$li(strong("MetaMEx 4.2508 – Aug 2025:"), " Machine learning for metadata imputation, switch to individual-level data"),
                               tags$li(strong("MetaMEx 3.2208 – Sept 2022:"), " RNA-seq update, bug fixes"),
                               tags$li(strong("MetaMEx 2.2101 – Jan 2021:"), " More studies added"),
                               tags$li(strong("MetaMEx 1.1805 – May 2018:"), " Initial release")
                       ),
                       
                       #--------------------------------------------------------------------
                       h3("Download"),

                       downloadButton("download_human_metadata", label = "Human - Metadata"),
                       downloadButton("download_human_genelist", label = "Human - Gene list"),
                       downloadButton("download_human_matrix", label = "Human - Gene matrix"),
                       downloadButton("download_mouse_metadata", label = "Mouse - Metadata"),
                       downloadButton("download_mouse_genelist", label = "Mouse - Gene list"),
                       downloadButton("download_mouse_matrix", label = "Mouse - Gene matrix"),
                       
                       #--------------------------------------------------------------------
                       h3("Citation"),
                       HTML('
    <div style="margin: 0rem;">
      <div>
        Pillon NJ, Gabriel BM, Dollet L, Smith JAB, Sardón Puig L, Botella J, Bishop DJ, Krook A, Zierath JR.
        <strong>Transcriptomic profiling of skeletal muscle adaptations to exercise and inactivity.</strong>
        <em>Nat Commun.</em> 2020 Jan 24;11(1):470.
        <a href="https://doi.org/10.1038/s41467-019-13869-w" target="_blank">https://doi.org/10.1038/s41467-019-13869-w</a>
      </div>
      <div style="margin-top: 0.5rem;">
        <script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>
        <span class="__dimensions_badge_embed__" data-id="pub.1124285483" data-legend="always"></span>
      </div>
    </div>
  '),
                       
                       #--------------------------------------------------------------------
                       h3("Datasets"),
                       p(em("Think we've missed a study? ",
                            a("Email us!", href = "mailto:nicolas.pillon@ki.se", target = "_blank"))),
                       
                       tabsetPanel(
                         tabPanel("Human Studies", 
                                  DT::dataTableOutput("reftable_human")
                         ),
                         tabPanel("Mouse Studies", 
                                  DT::dataTableOutput("reftable_mouse")
                         )
                       )
                       
)

