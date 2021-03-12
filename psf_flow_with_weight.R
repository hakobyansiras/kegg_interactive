psf_flow_with_weight <- function (g, node.ordering, sink.nodes, split = TRUE, sum = FALSE) 
{
  node.order <- node.ordering$node.order
  node.rank <- node.ordering$node.rank
  nods = names(node.order)
  symb.exprs = vector("list")
  eval.exprs = vector("list")
  E = data.frame(as.numeric(graph::nodeData(g, nods, attr = "expression")), 
                 row.names = nods)
  I = matrix(data = NA, nrow = length(nods), ncol = length(nods))
  W = matrix(data = NA, nrow = length(nods), ncol = length(nods))
  rownames(I) = nods
  colnames(I) = nods
  rownames(W) = nods
  colnames(W) = nods
  for (node in nods) {
    l = length(graph::edgeL(g)[[node]]$edges)
    for (e in 1:l) {
      from = node
      to = graph::nodes(g)[graph::edgeL(g)[node][[1]]$edges][e]
      if (!is.na(to) && length(to) > 0) {
        weight = graph::edgeData(g, from = from, to = to, 
                                 attr = "weight")[[1]]
        W[from, to] = weight
        impact = graph::edgeData(g, from = from, to = to, 
                                 attr = "impact")[[1]]
        I[from, to] = impact
      }
    }
  }
  recalc.sinks = F
  if (is.null(sink.nodes)) {
    recalc.sinks = T
    cat("\nsink nodes were not supplied. Those will be recalculated ...\n")
  }
  for (i in 1:length(nods)) {
    parent.nodes <- graph::inEdges(nods[i], g)[[1]]
    node = nods[i]
    if (length(parent.nodes) > 0) {
      node.exp = E[i, 1]
      in.signal <- unlist(graph::nodeData(g, parent.nodes, 
                                          attr = "signal"))
      pi = which(nods %in% parent.nodes)
      node.signal <- graph::nodeData(g, nods[i], attr = "signal")[[1]]
      impact = I[parent.nodes, nods[i]]
      weight = W[parent.nodes, nods[i]]
      if (sum) {
        node.signal = (node.exp) + sum(in.signal * weight * 
                                         impact)
      }
      else {
        if (split) {
          a = 2000
          proportion = 1
          if (sum(in.signal) != 0) {
            proportion = in.signal/sum(in.signal)
            node.signal <- sum((proportion * weight * 
                                  node.exp) *  a*(2/(1+exp((-2*in.signal^impact)/a)) - 1))
          }
        }
        else {
          proportion = 1
          node.signal <- node.exp * weight * prod(in.signal^impact)
        }
      }
      graph::nodeData(g, nods[i], attr = "signal") <- node.signal
      in.sign.impacts = vector("list")
      signal.base.denoms = vector("list")
      if (length(parent.nodes) > 1) {
        for (parent in parent.nodes) {
          if (!is.null(eval.exprs[[parent]])) {
            in.sign.impacts[[parent]] = sprintf("(%s)^(1+I['%s','%s'])", 
                                                paste("eval.exprs[['", parent, "']]", sep = ""), 
                                                parent, node)
            signal.base.denoms[[parent]] = sprintf("(%s)", 
                                                   paste("eval.exprs[['", parent, "']]", sep = ""))
          }
          else {
            in.sign.impacts[[parent]] = sprintf("(E['%s',1])^(1+I['%s','%s'])", 
                                                parent, parent, node)
            signal.base.denoms[[parent]] = sprintf("(E['%s',1])", 
                                                   parent)
          }
        }
        signal.base.denom = paste(signal.base.denoms, 
                                  collapse = "+")
        signal.base = sprintf("E[%d,1]/(%s)", i, signal.base.denom)
        in.signal.impact = paste(in.sign.impacts, collapse = "+")
        eval.exprs[[node]] = sprintf("(%s)*(%s)", signal.base, 
                                     in.signal.impact)
      }
      else {
        parent = parent.nodes[1]
        if (!is.null(eval.exprs[[parent]])) {
          in.signal.impact = sprintf("(%s)^(I['%s','%s'])", 
                                     paste("eval.exprs[['", parent, "']]", sep = ""), 
                                     parent, node)
        }
        else {
          in.signal.impact = sprintf("(E['%s',1])^(I['%s','%s'])", 
                                     parent, parent, node)
        }
        signal.base = sprintf("E[%d,1]", i)
        eval.exprs[[node]] = sprintf("(%s)*(%s)", signal.base, 
                                     in.signal.impact)
      }
    }
    else {
      node.exp = E[i, 1]
      graph::nodeData(g, nods[i], attr = "signal") <- node.exp
      eval.exprs[[node]] = sprintf("E[%d,1]", i)
    }
    if (recalc.sinks) {
      if (length(parent.nodes) > 0) {
        child.nodes <- unlist(graph::edges(g, names(node.order)[i]))
        if (length(child.nodes) > 0) {
          child.node.ranks <- node.rank[child.nodes]
          rank.comp = child.node.ranks > node.rank[i]
          if (all(rank.comp)) 
            sink.nodes <- c(sink.nodes, names(node.order)[i])
        }
        else {
          sink.nodes <- c(sink.nodes, names(node.order)[i])
        }
      }
    }
  }
  signal.at.sink = NULL
  if (!is.null(sink.nodes)) {
    signal.at.sink <- unlist(graph::nodeData(g, sink.nodes, 
                                             attr = "signal"))
    all_signals <- unlist(graph::nodeData(g, nods, 
                                          attr = "signal"))
  }
  return(list(graph = g, order = node.ordering, sink.nodes = sink.nodes, 
              signal.at.sink = signal.at.sink, all_signals = all_signals, eval.exprs = eval.exprs, 
              I = I, E = E))
}
