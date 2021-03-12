kegg_data_downloader <- function(pathway_name) {
  
  pathway_id <- pathway_codes_new[pathway_name]
  
  kgml <- paste0("http://rest.kegg.jp/get/", pathway_id, "/kgml")
  
  image <- paste0("http://rest.kegg.jp/get/", pathway_id, "/image")
  
  kgml_path <- tempfile()
  
  download.file(kgml, destfile = kgml_path, method = "auto")
  
  image_path <- tempfile()
  
  download.file(image, destfile = image_path, method = "auto")
  
  img <- magick::image_read(image_path)
  
  graph <- psf::generate.kegg.collection.from.kgml(kgml_path)[[1]]
  
  graph::edgeDataDefaults(graph$graph, attr = "existence") <- "exist"
  graph::nodeDataDefaults(graph$graph, attr = "existence") <- "exist"
  graph::edgeDataDefaults(graph$graph, attr = "change_info") <- "no_change"
  graph::nodeDataDefaults(graph$graph, attr = "change_info") <- "no_change"
  graph::edgeDataDefaults(graph$graph, attr = "data_source") <- "kegg"
  graph::nodeDataDefaults(graph$graph, attr = "data_source") <- "kegg"  
  
  return(list(pathway = graph, image = img))
  
}

chemokine_test <- kegg_data_downloader(pathway_name = "Chemokine_signaling_pathway")
