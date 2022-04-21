library(shiny)
library(shinyjs)
library(shinyBS)
library(visNetwork)
library(shinyjqui)

# construction_edges <- fread(file = "edge_table.tsv", sep = "\t", verbose = F)
# load("app_data_new.RData")
load("app_data_small.RData")
construction_edges <- empty_edge_table

shinyServer(function(input, output, session) {
  
  #### reactive values ####
  v <- reactiveValues(construction_edge_table = empty_edge_table,
                      edge_adding_error = ""
  )
  
  
  #### edge data collector ####
  observeEvent(input$visnet_editor_graphChange, {
    
    if(input$visnet_editor_graphChange$cmd == "addEdge") {
      
      updateSelectizeInput(session, inputId = "edge_type", selected = "")
      updateSelectizeInput(session, inputId = "edge_mode", selected = "")
      
      show('edge_attrs')
      
    } 
    if(input$visnet_editor_graphChange$cmd == "editEdge") {
      
      selected_edge_data <- v$construction_edge_table[from == input$visnet_editor_graphChange$from & to == input$visnet_editor_graphChange$to]
      
      updateSelectizeInput(session, inputId = "edge_type", selected = names(which(arrow_types == selected_edge_data$arrows.to.type)))
      
      updateSelectizeInput(session, inputId = "edge_mode", selected = names(which(line_types == selected_edge_data$dashes)))
      
      show('edge_attrs')
    }
    
    if(input$visnet_editor_graphChange$cmd == "deleteElements") {
      
      v$construction_edge_table <- v$construction_edge_table[id != input$visnet_editor_graphChange$edges[[1]]]
      
      # write.table(v$construction_edge_table, file = "edge_table.tsv", sep = "\t", quote = F, row.names = F)
      
    }
  })
  
  #### edge attra filling error text output ####
  output$edge_adding_error <- renderText({
    paste("<font color=\"#ff0000\"><b>", v$edge_adding_error, "</b></font>")
  })
  
  #### edge data submiting and saving ####
  observeEvent(input$submit_edge_attrs, {
    
    if(input$edge_type == "" | input$edge_mode == "") {
      v$edge_adding_error <- "Fill all fields"
      delay(3000, v$edge_adding_error <- "")
    } else {
      
      v$construction_edge_table <- rbind(v$construction_edge_table, list(input$visnet_editor_graphChange$id,
                                                                         input$visnet_editor_graphChange$from, 
                                                                         input$visnet_editor_graphChange$to, 
                                                                         line_colors[input$edge_type], line_types[input$edge_mode], TRUE, arrow_types[input$edge_type]))
      
      visNetworkProxy("visnet_editor") %>%
        visUpdateEdges(data.frame(id = input$visnet_editor_graphChange$id,
                                  from = input$visnet_editor_graphChange$from, 
                                  to = input$visnet_editor_graphChange$to, 
                                  color = line_colors[input$edge_type], dashes = line_types[input$edge_mode], 
                                  arrows.to.enabled = TRUE,
                                  arrows.to.type = arrow_types[input$edge_type]))
      
      
      # write.table(v$construction_edge_table, file = "edge_table.tsv", sep = "\t", quote = F, row.names = F)
      
      hide('edge_attrs')
    }
    
  })
  
  #### network submiting warning ####
  observeEvent(input$network_submit, {
    
    showModal(modalDialog(
      title = NULL,
      size = "m",
      footer = tagList(column(12, align="center", 
                              h5("After submiting the network you will not be able to get back to this task"),
                              actionButton("final_submit", label = "Submit"),
                              modalButton("Cancel")
      )
      ),
      easyClose = TRUE
    ))
    
  })
  
  #### network submiting ####
  observeEvent(input$final_submit, {
    removeModal()
    hide('visnet_editor_panel')
    hide('edge_attrs')
  })
  
  #### edge data editor ####
  # observe({
  #   if(!is.null(unlist(input$editor_clicked_node$edges))) {
  #     
  #     show('edge_attrs')
  #     
  #     visNetworkProxy("visnet_editor") %>%
  #       visGetEdges(input = "edge_list")
  #     
  #     edge <- c(input$edge_list[[unlist(input$editor_clicked_node$edges)]]$from, input$edge_list[[unlist(input$editor_clicked_node$edges)]]$to)
  #     
  #     
  #   } else {
  #     hide('edge_attrs')
  #   }
  # })
  
  #### network cunstructor rendering ####
  output$visnet_editor <-  renderVisNetwork({
    visNetwork(nodes = construction_nodes_unlabeled, 
               edges = construction_edges, width = "100%", height = "800px") %>%
      visEdges(smooth = FALSE) %>%
      visPhysics(enabled = F) %>%
      visInteraction(navigationButtons = TRUE, multiselect = F, selectConnectedEdges = F) %>%
      visLegend(addEdges = ledges, addNodes = lnodes, useGroups = FALSE, width = 0.1, stepY = 60, position = "right", zoom = FALSE) %>%
      visOptions(manipulation = list(enabled = TRUE, editNode = FALSE, addNode = FALSE, deleteNode = FALSE, initiallyActive = TRUE, 
                                     editEdge = FALSE
      )
      ) %>%
      visEvents(click = "function(nodes) {
                      console.info('click')
                      console.info(nodes)
                      Shiny.onInputChange('editor_clicked_node', {nodes : nodes.nodes, edges : nodes.edges});
                      ;}"
      )
  })
  
})
