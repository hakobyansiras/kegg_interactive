library(shiny)
library(shinyWidgets)
library(shinyjs)
library(visNetwork)

shinyUI(
  
  fluidPage(
    
    titlePanel("Network construction"),
    
    tags$script('
      $(document).mousedown(function(ev){
        if(ev.which == 3)
        {
          Shiny.onInputChange("right_click", Math.random());
        }
      });
    '),
    
    tags$script('
      $(document).mouseup(function(ev){
        if(ev.which == 1)
        {
          Shiny.onInputChange("left_click", Math.random());
        }
      });
    '),
    
    #### right click dialog prevention script ####
    tags$head(tags$script('document.addEventListener("contextmenu", function(e) {

                           e.preventDefault();
                           }, false);'
    )),
    
    #### pathway image ####
    fluidRow(
      
                   div(id = "visnet_editor_panel",
                       actionButton("network_submit", label = "Submit Network"),
                       visNetworkOutput("visnet_editor", height = "800px")
                   ),
                 useShinyjs(),
                 hidden(
                   div(id = "edge_attrs",
                       absolutePanel(
                         top = 500, right = 500,
                         width = 300,
                         draggable = TRUE,
                         fixed = TRUE,
                         cursor = "auto",
                         wellPanel(
                           style = "background-color: rgb(217, 217, 217, 0.8); z-index:2000;",
                           # actionButton("close_edge_", "X", style='position: absolute; right: 8px; top: 8px; text-align: center; width: 18px; height 15px; border-style: solid; border-radius: 10px; background-color: rgb(255, 0, 0, 0.9); font-size: 80%; padding: 0px; align :'),
                           selectizeInput("edge_type", label = "Interaction type", choices = c("Activation", "Inhibition", "Dissociation", ""), selected = ""),
                           selectizeInput("edge_mode", label = "Interaction state", choices = c("Direct", "Indirect", ""), selected = ""),
                           div(style = 'padding-left: 8px; padding-top: 8px;', htmlOutput('edge_adding_error')),
                           actionButton("submit_edge_attrs", label = "Save")
                         )
                       )
                   )
                 )
                 
    )
  )
)
