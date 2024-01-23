library(HDF5Array)
library(SummarizedExperiment)
library(igraph)
library(rhdf5)
library(Matrix)


#' Read the taxonomy graph from a h5
#' @param h5_file path to the h5 file containing the taxonomy. The taxonomy should be stored into path "/dendrogram"
#' @return an igraph
read_tax_from_h5 <- function(h5_file) {
  tax <- data.frame(h5read(h5_file,"dendrogram"),stringsAsFactors = FALSE,check.names = FALSE)
  e <- tax[tax$parent_accession!="",c("cell_set_accession","parent_accession")]
  v <- cbind(vid=tax$cell_set_accession,tax)
  g <- graph_from_data_frame(e,vertices = v)
  g
}


#' Compute ancestor matrix
#' @param tax an igraph
#' @return a square logical matrix. M[i,j] is true if there is a path from node i to node j in the graph
tax_ancestors <- function(tax) {
  is.finite(distances(tax,mode="out"))
}

# Wrap HDF5 dataset into a SummarizedExpermient
read_from_h5 <- function(h5_file) {
  meta <- as.data.frame(h5read(h5_file,"col_meta"))
  rownames(meta) <- meta$sample_name
  meta$sample_name <- NULL

  x <- HDF5Array(h5_file,"x",type="integer")
  rownames(x) <- as.vector(h5read(h5_file,"col_names"))
  colnames(x) <- as.vector(h5read(h5_file,"row_names"))
  
  cells <- intersect(rownames(meta),colnames(x))
  x <- x[,cells]
  
  e <- SummarizedExperiment(
    x,
    colData = meta[cells,],
    metadata = meta
  )
}






