parse.KGML_new <- function (kgml) 
{
  if (!file.exists(kgml)) 
    stop("Provided KGML file: ", kgml, ", does not exist")
  doc <- XML::xmlParseDoc(kgml)
  r <- XML::xmlRoot(doc)
  root.children = XML::xmlChildren(r)
  entry.nodes = root.children[which(names(root.children) == 
                                      "entry")]
  entry.attr.names = c("kegg.id", "kegg.name", "kegg.type", 
                       "kegg.link")
  names(entry.attr.names) = c("id", "name", "type", "link")
  graphics.attr.names = c("kegg.gr.name", "kegg.gr.fgcolor", 
                          "kegg.gr.bgcolor", "kegg.gr.type", "kegg.gr.x", "kegg.gr.y", 
                          "kegg.gr.width", "kegg.gr.height")
  names(graphics.attr.names) = c("name", "fgcolor", "bgcolor", 
                                 "type", "x", "y", "width", "height")
  nodelist.mat = matrix("NA", nrow = length(entry.nodes), 
                        ncol = length(c(entry.attr.names, graphics.attr.names)) + 
                          1)
  colnames(nodelist.mat) = c(entry.attr.names, graphics.attr.names, 
                             "components")
  i = 1
  for (entry in entry.nodes) {
    entry.attrs = XML::xmlAttrs(entry)
    for (name in names(entry.attr.names)) {
      if (name %in% names(entry.attrs)) 
        nodelist.mat[i, entry.attr.names[[name]]] = as.character(entry.attrs[[name]])
    }
    entry.children = XML::xmlChildren(entry)
    graphics = entry.children$graphics
    graphics.attrs = XML::xmlAttrs(graphics)
    for (name in names(graphics.attr.names)) {
      if (name %in% names(graphics.attrs)) 
        nodelist.mat[i, graphics.attr.names[[name]]] = as.character(graphics.attrs[[name]])
    }
    if (entry.attrs[["type"]] == "group") {
      components = entry.children[which(names(entry.children) ==
                                          "component")]
      if (length(components) > 0) {
        componentids = ""
        for (component in components) {
          id = as.character(XML::xmlGetAttr(component,
                                            "id"))
          componentids = paste(componentids, id, sep = ";")
        }
      }
      nodelist.mat[i, "components"] = componentids
    }
    if (graphics.attrs[["type"]] == "line") {
      coords = XML::xmlGetAttr(graphics, "coords")
      coords = as.numeric(unlist(strsplit(coords, split = ",")))
      xmean = base::mean(c(coords[1], coords[3]))
      ymean = base::mean(c(coords[2], coords[4]))
      nodelist.mat[i, "kegg.gr.x"] = xmean
      nodelist.mat[i, "kegg.gr.y"] = ymean
    }
    i = i + 1
  }
  rownames(nodelist.mat) = nodelist.mat[, "kegg.id"]
  g = graph::graphNEL(nodes = c(rownames(nodelist.mat)), edgemode = "directed")
  graph::nodeDataDefaults(g, attr = "genes") <- list()
  graph::nodeDataDefaults(g, attr = "expression") <- 1
  graph::nodeDataDefaults(g, attr = "signal") <- 1
  graph::nodeDataDefaults(g, attr = "type") <- "gene"
  graph::nodeDataDefaults(g, attr = "label") <- "NA"
  graph::nodeDataDefaults(g, attr = "components") <- "NA"
  for (attr in entry.attr.names) {
    graph::nodeDataDefaults(g, attr = attr) <- NA
  }
  for (attr in graphics.attr.names) {
    graph::nodeDataDefaults(g, attr = attr) <- NA
  }
  for (attr in entry.attr.names) {
    graph::nodeData(g, graph::nodes(g), attr = attr) <- nodelist.mat[, 
                                                                     attr]
  }
  for (attr in graphics.attr.names) {
    graph::nodeData(g, graph::nodes(g), attr = attr) <- nodelist.mat[, 
                                                                     attr]
  }
  graph::nodeData(g, graph::nodes(g), "components") = nodelist.mat[, 
                                                                   "components"]
  graph::nodeData(g, graph::nodes(g), "type") = nodelist.mat[, 
                                                             "kegg.type"]
  for (node in graph::nodes(g)) {
    type = nodelist.mat[node, "kegg.type"]
    if (type == "gene") {
      kegg.name = nodelist.mat[node, "kegg.name"]
      if (!is.na(kegg.name)) {
        entrezIds = unlist(base::strsplit(kegg.name, 
                                          "hsa:", fixed = T))
        entrezIds = unlist(base::strsplit(entrezIds, 
                                          " ", fixed = T))
        na.ind = which(entrezIds == "")
        if (length(na.ind) > 0) 
          entrezIds = entrezIds[-na.ind]
      }
      if (length(entrezIds) > 0) 
        graph::nodeData(g, node, "genes") = list(entrezIds)
    }
    kegg.gr.name = nodelist.mat[node, "kegg.gr.name"]
    gr.names = unlist(strsplit(kegg.gr.name, ", "))
    na.ind = which(gr.names == "")
    if (length(na.ind) > 0) 
      gr.names = gr.names[-na.ind]
    label = gr.names[1]
    graph::nodeData(g, node, "label") = label
  }
  edge.attrs = list(impact = "impact", type = "type", subtype1 = "subtype1", 
                    subtypeValue1 = "subtypeValue1", subtypeValue2 = "subtypeValue2", 
                    subtype2 = "subtype2")
  graph::edgeDataDefaults(g, attr = edge.attrs$impact) <- 1
  graph::edgeDataDefaults(g, attr = edge.attrs$type) <- NA
  graph::edgeDataDefaults(g, attr = edge.attrs$subtype1) <- NA
  graph::edgeDataDefaults(g, attr = edge.attrs$subtypeValue1) <- NA
  graph::edgeDataDefaults(g, attr = edge.attrs$subtype2) <- NA
  graph::edgeDataDefaults(g, attr = edge.attrs$subtypeValue2) <- NA
  relations = root.children[which(names(root.children) == 
                                    "relation")]
  for (relation in relations) {
    entry1 = as.character(XML::xmlGetAttr(relation, "entry1"))
    entry2 = as.character(XML::xmlGetAttr(relation, "entry2"))
    type = as.character(XML::xmlGetAttr(relation, "type"))
    g = psf::add.kegg.edge(entry1, entry2, type, subtype = "n/a", 
                           edge.attrs, g)
    relation.children = XML::xmlChildren(relation)
    subtypes = relation.children[which(names(relation.children) == 
                                         "subtype")]
    if (length(subtypes) >= 1) {
      subtype1atts = XML::xmlAttrs(subtypes[[1]])
      if (entry2 %in% graph::edges(g)[[entry1]]) {
        graph::edgeData(g, entry1, entry2, edge.attrs$subtype1) = as.character(subtype1atts[["name"]])
        graph::edgeData(g, entry1, entry2, edge.attrs$subtypeValue1) = as.character(subtype1atts[["value"]])
        if (length(subtypes) == 2) {
          subtype2atts = XML::xmlAttrs(subtypes[[2]])
          graph::edgeData(g, entry1, entry2, edge.attrs$subtype2) = as.character(subtype2atts[["name"]])
          graph::edgeData(g, entry1, entry2, edge.attrs$subtypeValue2) = as.character(subtype2atts[["value"]])
        }
      }
      else {
        cat("Edge ", entry1, ":", entry2, " not found\n")
      }
    }
    else {
      cat("No subtypes found for edge ", entry1, "|", 
          entry2, "\n")
    }
  }
  reactions = root.children[which(names(root.children) == 
                                    "reaction")]
  for (reaction in reactions) {
    reaction.id = XML::xmlGetAttr(reaction, "id")
    reaction.name = XML::xmlGetAttr(reaction, "name")
    reaction.entries = XML::xmlChildren(reaction)
    subst.ind = which(names(reaction.entries) == "substrate")
    prod.ind = which(names(reaction.entries) == "product")
    if (length(subst.ind) > 1) 
      cat("\nreaction ", reaction.name, " had ", length(subst.ind), 
          " substrates\n")
    for (s.id in subst.ind) {
      substrate = reaction.entries[[s.id]]
      substrate.id = XML::xmlGetAttr(substrate, "id")
      entry1 = substrate.id
      entry2 = reaction.id
      g = psf::add.kegg.edge(entry1, entry2, "ECrel", 
                             "reaction", edge.attrs, g)
    }
    if (length(prod.ind) > 1) 
      cat("\nreaction ", reaction.name, " had ", length(prod.ind), 
          " products\n")
    for (p.id in prod.ind) {
      product = reaction.entries[[p.id]]
      product.id = XML::xmlGetAttr(product, "id")
      entry1 = reaction.id
      entry2 = product.id
      g = psf::add.kegg.edge(entry1, entry2, "ECrel", 
                             "reaction", edge.attrs, g)
    }
  }
  # for (node in graph::nodes(g)) {
  #   if (graph::nodeData(g, node, "type") == "group") {
  #     components = as.character(graph::nodeData(g, node, 
  #                                               "components"))
  #     components = unlist(strsplit(components, split = ";"))
  #     na.ind = which(components == "")
  #     if (length(na.ind) > 0) 
  #       components = components[-na.ind]
  #     g = process.groupNode(node, components, g, edge.attrs)
  #     g = graph::removeNode(node, g)
  #   }
  # }
  g = process.compounds(g, edge.attrs)
  g = correctEdgeDirections(g, edge.attrs)
  g = set.edge.impacts(g)
  return(g)
}


g <- parse.KGML_new("pathway_kgmls/kgml.122")
