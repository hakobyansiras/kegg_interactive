visnet_creator <- function(graphical_data, hide_processed_edges, node_colors = NULL) {
  graphical_data$edge_coords <- graphical_data$edge_coords[which(graphical_data$edge_coords$lty == "solid"),]
  graphical_data$node_coords <- graphical_data$node_coords[which(graphical_data$node_coords$exist),]
  
  node_shapes <- unname(sapply(graphical_data$node_coords$node_name, function(x) {
    if(grepl("cpd",x)) {
      "dot"
    } else {
      "box"
    }
  }))
  
  if(!is.null(node_colors)) {
    color <- unname(sapply(graphical_data$node_coords$node_id, function(x) {
      if(x %in% node_colors$node_id) {
        node_colors[which(node_colors$node_id == x),"col"]
      } else {
        "#BFFFBF"
      }
    }))
    
    font_color <- unname(sapply(graphical_data$node_coords$node_id, function(x) {
      if(x %in% node_colors$node_id) {
        node_colors[which(node_colors$node_id == x),"text_col"]
      } else {
        "#000000"
      }
    }))
    
  } else {
    color <- unname(sapply(graphical_data$node_coords$sink, function(x) {
      if(x) {
        "#0099cc"
      } else {
        "#BFFFBF"
      }
    }))
    
    font_color <- rep("#000000", nrow(graphical_data$node_coords))
    
  }
  
  size <- unname(sapply(graphical_data$node_coords$node_name, function(x) {
    if(grepl("cpd",x)) {
      10
    } else {
      25
    }
  }))
  
  nodes <- data.frame(id = graphical_data$node_coords$node_id,
                      label = graphical_data$node_coords$gr_name,
                      shape = node_shapes, color = color, 
                      title = graphical_data$node_coords$hover_name,
                      font.size = rep(22, nrow(graphical_data$node_coords)), size = size,
                      font.color = font_color,
                      x = (graphical_data$node_coords$x_start + graphical_data$node_coords$x_end)/2,
                      y = (graphical_data$node_coords$y_start + graphical_data$node_coords$y_end)/2
  )
  
  arrows_type <- c("arrow","bar")
  names(arrows_type) <- c("simple", "T")
  
  if(hide_processed_edges) {
    graphical_data$edge_coords$col[which(graphical_data$edge_coords$label == "Done")] <- "white"
    
    graphical_data$edge_coords$label[which(graphical_data$edge_coords$label == "Done")] <- ""
    
  }
  
  edges <- data.frame(from = graphical_data$edge_coords$from, to = graphical_data$edge_coords$to,
                      color = graphical_data$edge_coords$col,
                      label = graphical_data$edge_coords$label,
                      arrows.to.enabled = rep(TRUE, length(graphical_data$edge_coords$label)),
                      arrows.to.type = arrows_type[graphical_data$edge_coords$arr.type]
  )
  
  
  return(list(nodes = nodes, edges = edges))
  
}

#### visnet node click ####

observe({
  if(!is.null(unlist(input$clicked_node$nodes))) {
    # v$selected_node_data <- as.character(v$subnetwork[[2]][unlist(input$clicked_node$nodes),2])
    
    updateNumericInput(session, "fc_num", label = paste0(unlist(unname(graph::nodeData(v$psf_and_colors$psf_graph[[1]][[v$pathway_name]]$graph, unlist(input$clicked_node$nodes),attr = "label"))), " Exp value"),
                       value = round(unlist(unname(graph::nodeData(v$psf_and_colors$psf_graph[[1]][[v$pathway_name]]$graph, unlist(input$clicked_node$nodes),attr = "expression"))), digits = 3)
    )
    
    
  } else {
    updateNumericInput(session, "fc_num", label = "", value = 0)
  }
})

## application needed static variables

graphical_data_generator <- function(pathway, include_changes = FALSE) {
  
  entrez_id <- unname(sapply(pathway$graph@nodes, function(y) {
    ifelse(is.null(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))),
           as.character(unlist(graph::nodeData(pathway$graph, y, attr = "label"))), 
           paste(as.character(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))), collapse = ","))
  }))
  
  hover_name <- unname(sapply(pathway$graph@nodes, function(y) {
    ifelse(is.null(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))),
           paste0(as.character(unlist(graph::nodeData(pathway$graph, y, attr = "label"))), " (",
                  unlist(strsplit(kegg_compounds_to_full_name[as.character(unlist(graph::nodeData(pathway$graph, y, attr = "label"))),], split = ";"))[1],
                  ")"
           ),
           paste0(
             paste0(as.character(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))), 
                    paste0("(", entrez_to_symbol[as.character(unlist(graph::nodeData(pathway$graph, y, attr = "genes"))),], ")")),
             collapse = ", "
           )
    )
  }))
  
  kegg_coords_new <- data.frame(
    pathway = rep(pathway$attrs$name, length(graph::nodeData(pathway$graph, attr = "kegg.gr.x"))),
    x_start = as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.x"))) - as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.width")))*0.5, 
    y_start = as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.y"))) - as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.height")))*0.5,
    x_end  = as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.x"))) + as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.width")))*0.5,
    y_end =  as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.y"))) + as.integer(unlist(graph::nodeData(pathway$graph, attr = "kegg.gr.height")))*0.5,
    node_name = as.character(unlist(graph::nodeData(pathway$graph, attr = "kegg.name"))), 
    gr_name = as.character(unlist(graph::nodeData(pathway$graph, attr = "label"))),
    node_id = as.character(unlist(graph::nodeData(pathway$graph, attr = "kegg.id"))),
    entrez_id = entrez_id, hover_name = hover_name,
    sink = (unlist(graph::nodeData(pathway$graph, attr = "kegg.id")) %in% pathway$sink.nodes), 
    exist = (unlist(graph::nodeData(pathway$graph, attr = "existence")) == "exist"),
    changed = (unlist(graph::nodeData(pathway$graph, attr = "data_source")) != "kegg"),
    node_class = as.character(unlist(graph::nodeData(pathway$graph, attr = "type"))), stringsAsFactors = F 
  )
  
  kegg_arrows_type <- c("simple", "simple", "simple", "simple", "T", "simple", "simple", "simple", "T", "simple", "", "simple", "simple", "T", "simple", "simple")
  names(kegg_arrows_type) <- c("activation", "binding/association", "compound", "dephosphorylation", "dissociation", "expression", "glycosylation", "indirect effect", "inhibition", "missing interaction", "n/a", "phosphorylation", "reaction", "repression", "state change", "ubiquitination" )
  
  line_col <- c("red", "red", "red", "red", "blue","red", "red", "red", "blue", "red", "red", "red", "red", "blue", "red", "red")
  names(line_col) <- c("activation", "binding/association", "compound", "dephosphorylation", "dissociation", "expression", "glycosylation", "indirect effect", "inhibition", "missing interaction", "n/a", "phosphorylation", "reaction", "repression", "state change", "ubiquitination" )
  
  # from <- as.character(graph::edgeMatrix(pathway$graph)[1,])
  # to <- as.character(graph::edgeMatrix(pathway$graph)[2,])
  # 
  # if(sum(is.na(c(match(from, kegg_coords_new$node_id), match(to, kegg_coords_new$node_id)))) > 0) {
  #   splitted_interactions <- strsplit(graph::edgeNames(pathway$graph), split = "~")
  #   from <- as.character(sapply(splitted_interactions, "[[", 1))
  #   to <- as.character(sapply(splitted_interactions, "[[", 2))
  # }
  
  if(length(names(pathway$graph@edgeData@data)) > 0) {
    splitted_interactions <- strsplit(names(pathway$graph@edgeData@data), split = "|", fixed = T)
    from <- as.character(sapply(splitted_interactions, "[[", 1))
    to <- as.character(sapply(splitted_interactions, "[[", 2))
  } else {
    return(list(node_coords = kegg_coords_new, edge_coords = NULL))
  }
  
  from_index <- match(from, kegg_coords_new$node_id)
  to_index <- match(to, kegg_coords_new$node_id)
  
  col <- line_col[unname(unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "subtype1")))]
  col[which(unname(unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "data_source"))) != "kegg" & col == "red")] <- "purple"
  col[which(unname(unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "data_source"))) != "kegg" & col == "blue")] <- "yellow"
  
  lty_type <- c("solid", "solid", "dotted")
  names(lty_type) <- c("exist", "added", "removed")
  
  edge_label_type <- c("", "In progress", "Done")
  names(edge_label_type) <- c("none", "in_progress", "done")
  
  arrow_data <- data.frame(
    from = from, to = to,
    x0 = (kegg_coords_new$x_start[from_index] + kegg_coords_new$x_end[from_index])/2,
    x1 = (kegg_coords_new$x_start[to_index] + kegg_coords_new$x_end[to_index])/2,
    y0 = (kegg_coords_new$y_start[from_index] + kegg_coords_new$y_end[from_index])/2,
    y1 = (kegg_coords_new$y_start[to_index] + kegg_coords_new$y_end[to_index])/2,
    arr.type = kegg_arrows_type[unname(unlist(graph::edgeData(pathway$graph, from = from, to = to,attr = "subtype1")))],
    label = edge_label_type[unname(unlist(graph::edgeData(pathway$graph, from = from, to = to, attr = "expansion_curation_stage")))],
    col = col, lty = lty_type[unname(unlist(graph::edgeData(pathway$graph, from = from, to = to,attr = "existence")))],
    stringsAsFactors = F
  )
  
  if(!include_changes) {
    arrow_data <- arrow_data[which(arrow_data$lty != "dotted"),]
  }
  
  return(list(node_coords = kegg_coords_new, edge_coords = arrow_data))
}


graphical_data = graphical_data_generator(edited_pathways_new[["Chemokine_signaling_pathway"]])

graphical_data$node_coords$gr_name <- chemokine_pathway_with_graphic_names$kegg_name

graphical_data$node_coords <- graphical_data$node_coords[-c(5,26,56),]

# graphical_data$node_coords <- as.data.table(graphical_data$node_coords)

# graphical_data$node_coords$gr_name <- as.list(graphical_data$node_coords$gr_name)

graphical_data$node_coords["11", "gr_name"] <- expression(paste("PLC", beta))

graphical_data$node_coords["17", "gr_name"] <- expression(paste("G", alpha, "i"))

graphical_data$node_coords["18", "gr_name"] <- expression(paste("G", beta, "y"))

graphical_data$node_coords["24", "gr_name"] <- expression(paste(beta, "-arrestin"))

graphical_data$node_coords["36", "gr_name"] <- expression(paste("PKC", zeta))

graphical_data$node_coords["52", "gr_name"] <- expression(paste("PKC", zeta))


setNames(object = c(expression(paste("PLC", beta)), expression(paste("G", alpha, "i")), expression(paste("G", beta, "y")), expression(paste(beta, "-arrestin")), expression(paste("PKC", zeta)), expression(paste("PKC", zeta))),
         nm = c("11", "17", "18", "24", "36", "52")
         )

node_coords_extended <- cbind(graphical_data$node_coords, t(
  sapply(graphical_data$node_coords$sink, function(x) {
    
    if(x) {
      c("blue", "dashed")
    } else {
      c("red", "solid")
    }
    
  })
))

colnames(node_coords_extended)[15:16] <- c("border_color", "lty_type")

node_coords_extended$border_color <- as.character(node_coords_extended$border_color)

node_coords_extended$lty_type <- as.character(node_coords_extended$lty_type)

node_coords_extended$x_center <- node_coords_extended$x_start + (node_coords_extended$x_end - node_coords_extended$x_start)/2

node_coords_extended["10","gr_name"] <- "PI3K IA"

node_coords_extended["69","gr_name"] <- "PI3K IB"

node_coords_extended["18","gr_name"] <- "Gβ"

greek_maping_nodes <- node_coords_extended[c("11", "17", "18", "24", "36", "52"),]

greek_maping_nodes$gr_name <- c(expression(paste("PLC", beta)), expression(paste("G", alpha, "i")), expression(paste("G", beta)), expression(paste(beta, "-arrestin")), expression(paste("PKC", zeta)), expression(paste("PKC", zeta)))

node_coords_extended <- node_coords_extended[which(!(node_coords_extended$node_id %in% c("11", "17", "18", "24", "36", "52"))),]

graphical_data$node_coords <- node_coords_extended


network_edge_proxy <- dataTableProxy('fc_table_out', deferUntilFlush = FALSE)

observeEvent(input$fc_table_out_cell_edit, {
  
  v$v$fc_table$FC[input$fc_table_out_cell_edit$row] <- input$fc_table_out_cell_edit$value
  
})



default_node_colors <- data.frame(node_id = graphical_data$node_coords$node_id[graphical_data$node_coords$node_class == "gene"], 
                                  col = rep("#BFFFBF", 57),
                                  text_col = rep("#000000", 57),
                                  stringsAsFactors = F
)


c("β", "α", "ζ")


plot.new()
rect(c(0,0.2), 1, expression(paste("G", alpha, "i")))

textbox(c(0,0.2), 1, c("many words","more words","why not?",
                       "keep going",rep("and going",10)))
