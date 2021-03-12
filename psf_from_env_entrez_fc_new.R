psf.from.env.entrez.fc_new <- function (entrez.fc, kegg.collection, split = TRUE, calculate.significance = T, 
          bst.steps = 200, sum = FALSE) 
{
  psf.results.collection = list()
  psf.results.collections = list()
  for (c in 1:ncol(entrez.fc)) {
    for (i in 1:length(kegg.collection)) {
      cat("\nKEGG collection number: ", i, " exp matrix column: ", 
          c)
      pathway = names(kegg.collection)[i]
      if (!length(kegg.collection[[i]]) == 0) {
        entrez.column = as.matrix(entrez.fc[, c])
        g = map.gene.data(kegg.collection[[i]]$graph, 
                          entrez.column)
        psf.results.collection[[pathway]] = psf_flow_with_weight(g, 
                                                     kegg.collection[[i]]$order, kegg.collection[[i]]$sink.nodes, 
                                                     split, sum = sum)
        psf.results.collection[[pathway]] = c(psf.results.collection[[pathway]], 
                                              attrs = list(kegg.collection[[i]]$attrs), 
                                              order = list(kegg.collection[[i]]$order))
      }
    }
    if (calculate.significance) {
      cat("\nPerforming bootstrap calculations with", 
          bst.steps, " steps\n")
      psf.results = bootstrap.significance(psf.results.collection, 
                                           entrez.fc, bst.steps)
      psf.results.processed = process.psf.results(psf.results)
    }
    else {
      psf.results.processed = psf.results.collection
    }
    psf.results.collections[[c]] = psf.results.processed
  }
  return(psf.results.collections)
}

