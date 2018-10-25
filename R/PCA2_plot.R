#' Make Scatter Plot of PCA
#'
#' This functions creates a scatter plot of most significant Principal Components
#'
#' @export
#' @param file_1 main genomic file
#' @param file_2 metadata file
PCA2_plot <- function(file_1,file_2,loadings=FALSE,...) {
  library(plotly)

  if(!grepl(".csv$", file_1)){
    stop("Uploaded genomic file must be a .csv file!")
  }
  if(!grepl(".csv$", file_2)){
    stop("Uploaded meta_file must be a .csv file!")
  }

  genomic_data = read.csv(file_1)
  metadata = read.csv(file_2)

  #make gene names as row names
  genomic_data = as.data.frame(sapply(genomic_data,as.numeric))
  genomic_data = na.omit(genomic_data)
  rownames(genomic_data) = make.names(genomic_data[,2],unique = TRUE)
  genomic_data = genomic_data[,3:ncol(genomic_data)]

  t_genomic_data = as.data.frame(t(genomic_data))
  t_genomic_data = na.omit(t_genomic_data)

  genomic_data.pca = prcomp(genomic_data,scale. = TRUE,center = TRUE)


  if(!loadings){
    pca_graph_1 = plot_ly(as.data.frame(genomic_data.pca$x), x = ~PC1, y = ~PC2,text = paste("Gene : ", rownames(genomic_data))) %>%
    layout(
      title = "Rotate Data into PC axis"
      )

    htmlwidgets::saveWidget(pca_graph_1, "pca_graph_1.html", selfcontained = FALSE)
  }else{
    pca_graph_2 = plot_ly(as.data.frame(genomic_data.pca[2]), x = ~rotation.PC1, y = ~rotation.PC2,mode="lines+marker",color = as.character(metadata$Time),text = paste("Sample : ",metadata$sIdx,"\nTime : ",metadata$Time," ",metadata$Unit)) %>%
    layout(
      title = "Loadings graph with Principal Components",
      scene = list(
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2")

      ))
    htmlwidgets::saveWidget(pca_graph_2, "pca_graph_2.html", selfcontained = FALSE)
  }



}







