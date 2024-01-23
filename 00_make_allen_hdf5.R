
source("lib_hdf5.R")


#' Read json file storing the allen taxonomy and convert into a data.frame
#' @param path to the json file containing the taxonomy in allen20 format
#' @return a data.frame
read_dend_json <- function(json_file) {
  x <- rjson::fromJSON(file=json_file,simplify = FALSE)

  # recursive function converting into dcf format
  to_dcf <- function(n,path="root") {
    # Check if the node is either a leaf or an internal node
    stopifnot(xor(identical(names(n),"leaf_attributes"),identical(names(n),c("node_attributes","children"))))
    
    # retrieve attributes
    att <- n[[grep("_attributes",names(n))]][[1]]
    att$path <- path # add path
    att <- paste(names(att),att,sep=":",collapse="\n")
    l <- mapply(to_dcf,n$children,file.path(path,seq_along(n$children)))
    c(att,unlist(l))
  }
  
  dcf <- to_dcf(x)
  dcf <- read.dcf(textConnection(paste(dcf,collapse="\n\n")))
  dcf <- as.data.frame(dcf)
  dcf$parent_accession <- dcf$cell_set_accession[match(dirname(dcf$path),dcf$path)]
  dcf$path <- NULL
  
  # write and reparse the data.frame to convert numeric values
  local({
    con <- textConnection(NULL,"w")
    write.table(dcf,file=con,sep="\t",quote=FALSE,na = "")
    read.table(text=textConnectionValue(con),sep="\t",quote="",comment="",check.names = FALSE,stringsAsFactors = FALSE)
  })
}


#' Generate a file in HDF5 format with an allen20 dataset
#' @param dir path to the directory containing the dataset files (matrix.csv.gz, metadata.csv.gz, tsne.csv.gz, dend.json.gz)
#' @param h5file name of the hdf5 file to create
make_allen_h5 <- function(dir,h5file) {
  # Create HDF5 assay
  h5_from_csv(gzfile(file.path(dir,"matrix.csv.gz")),h5file)
  
  # Add cells meta-data
  meta <- read.csv(file.path(dir,"metadata.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names=NULL)
  meta <- meta[sapply(meta,function(x) length(unique(x)))>1]
  h5write(meta,h5file,"col_meta",DataFrameAsCompound=FALSE)

  # Add tsne coordinates
  tsne <- read.csv(file.path(dir,"tsne.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names=NULL)
  tsne <- tsne[match(h5read(h5file,"col_meta/sample_name"),tsne$sample_name),]
  h5write(tsne$Lim1,h5file,"col_meta/tsne_1")
  h5write(tsne$Lim2,h5file,"col_meta/tsne_2")
  
  # Add dendrogram
  # dend <- read_dend_json(file.path(dir,"dend.json"))
  # rhdf5::h5write(dend,h5file,"dendrogram",DataFrameAsCompound=FALSE)
}
make_allen_h5("mouse-whole-cortex-and-hippocampus-10x_2022/","mouse-whole-cortex-and-hippocampus-10x_2022.h5")
make_allen_h5("mouse-whole-cortex-and-hippocampus-smart-seq_2022/","mouse-whole-cortex-and-hippocampus-smart-seq_2022.h5")
#make_allen_h5("data/mouse-whole-cortex-and-hippocampus-smart-seq/","out/mouse-whole-cortex-and-hippocampus-smart-seq.h5")
#make_allen_h5("data/human-multiple-cortical-areas-smart-seq/","out/human-multiple-cortical-areas-smart-seq.h5")
#make_allen_h5("data/human-m1-10x/","out/human-m1-10x.h5")




