function(group_graphics, kegg_pathway_graphics, pathway_name, pathway_image, highlight.genes = NULL, edge_mapping = FALSE, advanced_edge_ampping = FALSE, 
         show_changes = FALSE, edge_in_mode = FALSE, highlight_color = "green", opacity = 0)
  
  
function(group_graphics, node_graphics, highlight.genes = NULL, highlight_color = "green", opacity = 0, color.genes = NULL, 
         color_bar_psf_mode = F, col_legend_title, color_bar_lims = NULL, draw_color_bar = F)
  
  
kegg_node_mapper <- function(group_graphics, kegg_pathway_graphics, pathway_name, pathway_image, highlight.genes = NULL, 
                             edge_mapping = FALSE, advanced_edge_ampping = FALSE, show_changes = FALSE, 
                             edge_in_mode = FALSE, highlight_color = "green", opacity = 0, color.genes = NULL, 
                             color_bar_psf_mode = F, col_legend_title, color_bar_lims = NULL, draw_color_bar = F
                             ) {
  
  img <- image_draw(pathway_image)
  
  if(is.null(color.genes)) {
    for(i in 1:nrow(kegg_pathway_graphics$node_coords)) {
      
      if(kegg_pathway_graphics$node_coords[i,"changed"]) {
        lty_type <- "dotted"
      }
      
      if(kegg_pathway_graphics$node_coords[i,"sink"]) {
        border_color <- "blue"
        lty_type <- "dashed"
      } else {
        border_color <- "red"
        lty_type <- "solid"
      }
      
      if(show_changes) {
        if(kegg_pathway_graphics$node_coords[i,"changed"]) {
          border_color <- "purple"
        }
        if(!kegg_pathway_graphics$node_coords[i,"exist"]) {
          border_color <- "red"
          lty_type <- "dashed"
        }
      } else {
        if(!kegg_pathway_graphics$node_coords[i,"exist"]) {
          border_color <- NA
        }
      }
      
      if(kegg_pathway_graphics$node_coords[i,"node_class"] == "bio_event") {
        
        plotrix::textbox(x = c(kegg_pathway_graphics$node_coords[i,"x_start"], kegg_pathway_graphics$node_coords[i,"x_end"]),
                         y = c(kegg_pathway_graphics$node_coords[i,"y_start"], kegg_pathway_graphics$node_coords[i,"y_end"]), 
                         textlist = kegg_pathway_graphics$node_coords[i,"gr_name"], 
                         fill = "#32CD32", lty = lty_type, lwd = 2, leading = 0.3,
                         justify = 'c', border = border_color, adj = c(0,0.25))
        
      } else {
        rect( kegg_pathway_graphics$node_coords[i,"x_start"], 
              kegg_pathway_graphics$node_coords[i,"y_start"], 
              kegg_pathway_graphics$node_coords[i,"x_end"], 
              kegg_pathway_graphics$node_coords[i,"y_end"], 
              border = border_color, lty = lty_type, lwd=2)
      }
      
      
      if(kegg_pathway_graphics$node_coords[i,"node_id"] %in% highlight.genes[,"node_id"]) {
        rect( kegg_pathway_graphics$node_coords[i,"x_start"], 
              kegg_pathway_graphics$node_coords[i,"y_start"], 
              kegg_pathway_graphics$node_coords[i,"x_end"], 
              kegg_pathway_graphics$node_coords[i,"y_end"],   
              border = highlight_color, lty = "solid", lwd=2, col = adjustcolor( "#a3297a", alpha.f = opacity))
      }
      
      # text( kegg.pathway$pathway.info[[i]]$graphics$x, height - kegg.pathway$pathway.info[[i]]$graphics$y, kegg.pathway$pathway.info[[i]]$graphics$label, 
      #       cex = ifelse( nchar(kegg.pathway$pathway.info[[i]]$graphics$label)<=6,.3/2,.24/2), col = "black", family="sans" )
    }
  } else {
    
    ### node exp coloring
    
    node_graphics <- cbind(kegg_pathway_graphics$node_coords, t(
      sapply(kegg_pathway_graphics$node_coords$sink, function(x) {
        
        if(x) {
          c("blue", "dashed")
        } else {
          c("red", "solid")
        }
        
      })
    ))
    
    colnames(node_graphics)[15:16] <- c("border_color", "lty_type")
    
    node_graphics$border_color <- as.character(node_graphics$border_color)
    
    node_graphics$lty_type <- as.character(node_graphics$lty_type)
    
    node_graphics$x_center <- node_graphics$x_start + (node_graphics$x_end - node_graphics$x_start)/2
    
      
    rownames(color.genes) <- color.genes$node_id
      
    rownames(node_graphics) <- node_graphics$node_id
      
    rect( node_graphics$x_start, 
          node_graphics$y_start, 
          node_graphics$x_end, 
          node_graphics$y_end, 
          border = node_graphics$border_color, lty = node_graphics$lty_type, lwd=2)
    
      
    if(any(node_graphics$node_id %in% color.genes$node_id)) {
      coloring_set <- node_graphics[color.genes$node_id,]
      rect( coloring_set$x_start,
            coloring_set$y_start,
            coloring_set$x_end,
            coloring_set$y_end,
            border = coloring_set$border_color, lty = coloring_set$lty_type, lwd=2,
            col = adjustcolor( color.genes$col, alpha.f = 1)
      )
      
      text(x = coloring_set$x_center,
           y = coloring_set$y_start,
           labels = coloring_set$gr_name,
           col = color.genes$text_col, adj = c(0,0.5) + c(0.5, 1))
      
    }
    
  }
  
  ### scale color bar
  
  if(draw_color_bar) {
    if(color_bar_psf_mode) {
      color_legend_maker(x = 1000, y = 50, leg = 200, cols = c(pal1(10), pal2(10)), title = col_legend_title, lims = color_bar_lims, digits=3, prompt=FALSE,
                         lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
    } else {
      color_legend_maker(x = 1000, y = 50, leg = 200, cols = exp_pal(20), title = col_legend_title, lims = exp_color_all$exp_lims, digits=3, prompt=FALSE,
                         lwd=4, outline=TRUE, subtitle = "", fsize = 1.3)
    }
  }
  
  
  ## color grop nodes
  if(length(group_graphics) > 0 ) {
    lapply(group_graphics, function(z) {
      rect( z$graphics$x-z$graphics$width*0.5, 
            z$graphics$y+z$graphics$height*0.5, 
            z$graphics$x+z$graphics$width*0.5, 
            z$graphics$y-z$graphics$height*0.5, 
            border = "yellow", lty = "dashed", lwd=2)
    })
  }
  
  if(edge_mapping) {
    if(advanced_edge_ampping) {
      
      from <- kegg_pathway_graphics$node_coords[which(kegg_pathway_graphics$node_coords$node_id == highlight.genes[1,"node_id"]),]
      to <- kegg_pathway_graphics$node_coords[which(kegg_pathway_graphics$node_coords$node_id == highlight.genes[2,"node_id"]),]
      
      shape::Arrows(
        x0 = (from$x_start + from$x_end)/2,
        x1 = (to$x_start + to$x_end)/2,
        y0 = (from$y_start + from$y_end)/2,
        y1 = (to$y_start + to$y_end)/2,
        col = "red", 
        lwd=2, arr.length = 0.2, arr.type = "simple"
      )
      
    } else {
      
      if(!is.null(kegg_pathway_graphics$edge_coords)) {
        
        if(!show_changes) {
          kegg_pathway_graphics$edge_coords <- kegg_pathway_graphics$edge_coords[which(kegg_pathway_graphics$edge_coords$lty == "solid"),]
          # kegg_pathway_graphics$edge_coords$lty <- "solid"
          kegg_pathway_graphics$edge_coords$col[which(kegg_pathway_graphics$edge_coords$col == "purple")] <- "red"
          kegg_pathway_graphics$edge_coords$col[which(kegg_pathway_graphics$edge_coords$col == "yellow")] <- "blue"
        }
        
        if(nrow(highlight.genes) == 0) {
          shape::Arrows(
            x0 = kegg_pathway_graphics$edge_coords$x0, 
            x1 = kegg_pathway_graphics$edge_coords$x1, 
            y0 = kegg_pathway_graphics$edge_coords$y0, 
            y1 = kegg_pathway_graphics$edge_coords$y1,
            col = kegg_pathway_graphics$edge_coords$col,
            lty = kegg_pathway_graphics$edge_coords$lty,
            lwd=2, arr.length = 0.2, arr.type = "simple"
          )
        } else {
          if(edge_in_mode) {
            index <- unique(which(kegg_pathway_graphics$edge_coords$from %in% highlight.genes[,"node_id"] & kegg_pathway_graphics$edge_coords$to %in% highlight.genes[,"node_id"]))
          } else {
            index <- unique(c(which(kegg_pathway_graphics$edge_coords$from %in% highlight.genes[,"node_id"]), 
                              which(kegg_pathway_graphics$edge_coords$to %in% highlight.genes[,"node_id"])))
          }
          
          shape::Arrows(
            x0 = kegg_pathway_graphics$edge_coords$x0[index], 
            x1 = kegg_pathway_graphics$edge_coords$x1[index], 
            y0 = kegg_pathway_graphics$edge_coords$y0[index], 
            y1 = kegg_pathway_graphics$edge_coords$y1[index],
            col = kegg_pathway_graphics$edge_coords$col[index], 
            lty = kegg_pathway_graphics$edge_coords$lty[index],
            lwd=2, arr.length = 0.2, arr.type = "simple"
          )
        }
      }
    }
  }
  
  
  
  dev.off()
  
  return(img)
  
}



