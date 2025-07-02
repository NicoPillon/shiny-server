options(shiny.sanitize.errors=F) # sanitize errors

fluidPage(title="MetaMEx",
          useShinyjs(),
          style="padding:0%",
          
          #---------------------------------------------------------
          # CSS for style
          tags$head(includeCSS("../../www/style.css")),
          
          #---------------------------------------------------------
          # Code to send height to resizing iframe
          tags$head(
            tags$script(HTML("
    function sendHeight() {
      const height = document.documentElement.scrollHeight;
      parent.postMessage({ frameHeight: height }, '*');
    }

    Shiny.addCustomMessageHandler('resizeFrame', function(message) {
      sendHeight();
    });

    window.addEventListener('resize', sendHeight);

    document.addEventListener('shiny:connected', function() {
      sendHeight();
    });

    document.addEventListener('shiny:value', function() {
      setTimeout(sendHeight, 100); // slight delay to allow DOM to settle
    });

    // Optional: observe tab switches
    document.addEventListener('click', function() {
      setTimeout(sendHeight, 100);
    });
  "))
          ),
          
          #---------------------------------------------------------
          #include google analytics
          tags$head(includeHTML("google-analytics.html")),
          
          
          #---------------------------------------------------------
          # Navigation bar layout with multiple tabs
          navbarPage(id="inTabset",
                     
                     # Custom title: logo image + text
                     title = HTML('
      <div style="display: flex; align-items: center;margin: -10px;">
        <img src="../../www/img/snippet/MetaMEx.png" style="height: 40px; margin-right: 10px;">
        <span style="font-size: 20px; font-weight: bold; margin-right: 10px;">MetaMEx</span>
      </div>
    '),
                     
                     ## Tab 1 ######################################################################################        
                     tabHome,
                     
                     ## Tab 2 ######################################################################################
                     # navbarMenu("Analysis",
                     #            tabPanel("Human",
                     #                     tabMetaanalysis_human),
                     #            tabPanel("Mouse",
                     #                     tabMetaanalysis_mouse)
                     #            ),  
                     tabMetaanalysis_human,
                     tabMetaanalysis_mouse,
                     
                     ## Tab 4 ######################################################################################     
                     tabMethods,
          )

)
