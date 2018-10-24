#' Make Scatter Plot of PCA
#'
#' This functions creates a scatter plot of most significant Principal Components
#'
#' @export
#' @param file_1 main genomic file
#' @param file_2 metadata file
PCA2_plot <- function(file_1,file_2,...) {
  library(factoextra)
  library(ggfortify)

  if(!grepl(".csv$", file_1)){
    stop("Uploaded genomic file must be a .csv file!")
  }
  if(!grepl(".csv$", file_2)){
    stop("Uploaded meta_file must be a .csv file!")
  }

  genomic_data = read.csv(file_1)
  metadata = read.csv(file_2)

  #make gene names as row names
  rownames(genomic_data) = make.names(genomic_data[,1],unique = TRUE)
  genomic_data = genomic_data[,2:ncol(genomic_data)]
  genomic_data = as.data.frame(sapply(genomic_data,as.numeric))
  genomic_data = na.omit(genomic_data)
  t_genomic_data = as.data.frame(t(genomic_data))
  t_genomic_data = na.omit(t_genomic_data)

  t_genomic_data.pca = prcomp(t_genomic_data,scale. = TRUE,center = TRUE)

  #fviz_pca_var(t_genomic_data.pca,
  #             col.var = "contrib", # Color by contributions to the PC
  #             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  #             repel = TRUE )


  autoplot(t_genomic_data.pca,data=metadata,colour="Time")

}
