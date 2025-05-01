options(shiny.sanitize.errors=F) # sanitize errors

fluidPage(title="MetaMEx",
          useShinyjs(),
          style="padding:0%",
          
          # CSS for style
          tags$head(includeCSS("../../www/style.css")),
          
          #include google analytics & adjust progress bar position
          tags$head(includeHTML("google-analytics.html")),
          
          navbarPage(title="MetaMEx",  id="inTabset",
                     
                     ########################################################################################        
                     tabHome,
                     
                     ########################################################################################        
                     tabMetaanalysis_human,         
                     
                     ########################################################################################        
                     tabMetaanalysis_mouse,  
                     
                     ########################################################################################        
                     tabDatasets,
                     
                     ########################################################################################        
                     tabDownloads,
                     
                     ########################################################################################        
                     tabPanel("About", value="Tutorial",
                              tags$script(src = "myscript.js"),
                              fluidRow(style="background-color: #5b768e;margin:-18px -15px 0px -15px;
                                             padding:0.5% 1.5% 1% 1.3%",
                                       h5(style="color:white", "The online database of skeletal muscle 
                                                transcriptomic response to exercise and inactivity. On this website, 
                                                you will be able to explore how specific genes respond to acute exercise, 
                                                exercise training and inactivity in all available transcriptomic datasets.
                                                This section summarizes and explains the main analyses available on MetaMEx."
                                       ),
                              ),
                              fluidRow(style="padding:2% 5% 10% 5%",
                                       h3("FAQ"),
                                       includeMarkdown("annexes/tutorial.md"),
                                       tags$hr(),
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
')
                              )
                              
                              
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
