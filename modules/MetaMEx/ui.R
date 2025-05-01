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
                     tabPanel("Help", value="Tutorial",
                              tags$script(src = "myscript.js"),
                              fluidRow(style="background-color: #5b768e;margin:-18px -15px 0px -15px;
                                             padding:0.5% 1.5% 1% 1.3%",
                                       h2(style="color:white", "Welcome to MetaMEx"),
                                       h5(style="color:white", "The online database of skeletal muscle 
                                                transcriptomic response to exercise and inactivity. On this website, 
                                                you will be able to explore how specific genes respond to acute exercise, 
                                                exercise training and inactivity in all available transcriptomic datasets.
                                                This section summarizes and explains the main analyses available on MetaMEx."
                                       ),
                              ),
                              fluidRow(style="padding:2% 5% 10% 5%",
                                       includeMarkdown("annexes/tutorial.md"))
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
