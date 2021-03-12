map.gene.data_new <- function (g, entrez.fc) 
{
  gene.data <- graph::nodeData(g)
  gene.data <- lapply(gene.data, function(x, y) {
    genes.in.node = which(rownames(y) %in% x$genes)
    expression.values <- y[genes.in.node, ]
    if (length(expression.values) > 0) {
      if (x$type == "gene") {
        expression <- mean(expression.values)
      }
      else {
        if (x$type == "group")
          expression <- min(expression.values)
      }
    }
    else expression <- 1
    x$expression <- expression
    return(x)
  }, entrez.fc)
  for (i in 1:length(gene.data)) {
    if (gene.data[[i]]$type == "gene") 
      graph::nodeData(g, names(gene.data)[i], "expression") = gene.data[[i]]$expression
  }
  return(g)
}


map.gene.data_new(g, chronic_exps$fc_new)

map.gene.data_new(edited_pathways_new$Chemokine_signaling_pathway$graph, chronic_exps$fc_new)
