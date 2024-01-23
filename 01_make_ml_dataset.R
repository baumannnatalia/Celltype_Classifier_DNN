source("lib_allen.R")

write_h5_dataset <- function(e,h5f,chunkdim=c(nrow(e),1L),level=1L) {
  h5write(colnames(e),h5f,"/row_names")
  h5write(rownames(e),h5f,"/col_names")
  h5createGroup(h5f,"/col_meta")
  h5write(e$donor_sex_label,h5f,"/col_meta/donor_sex_label")
  h5write(e$class_label,h5f,"/col_meta/class_label")
  h5write(e$subclass_label,h5f,"/col_meta/subclass_label")
  h5write(e$cluster_label,h5f,"/col_meta/cluster_label")
  h5write(e$region_label,h5f,"/col_meta/region_label")
  h5write(e$cell_type_accession_label,h5f,"/col_meta/cell_type_accession_label")
  #h5write(metadata(e)$tax,h5f,"dendrogram",DataFrameAsCompound=FALSE)
  writeHDF5Array(assay(e),h5f,"/x",chunkdim=chunkdim,level=level)
}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Make a randomized train/valid/test sets of 10000 cells
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
make_shuffled_allen_datasets <- function() {
  # Load the data object
  e <- read_from_h5("mouse-whole-cortex-and-hippocampus-10x_2022.h5")
  #Filter cells
  e$subclass_label[grep("L4 ",e$cluster_label)] <- "L4 IT CTX"
  e <- e[,e$region_label=="SSp"]
  # Randomize cells
  set.seed(12345)
  e <- e[,sample(ncol(e))]
  write_h5_dataset(e,"allen10x_all_lev1_dim1_NB.h5",level=1L,chunkdim = c(nrow(e),1L))
  
  # Generate datasets
  # write_h5_dataset(e[,1:10000],"out/datasets/allen10x_train_lev1_dim1.h5",level=1L,chunkdim = c(nrow(e),1L)) # Generate a h5 of 69.1Mo, which is 4x faster to load
  # write_h5_dataset(e[,10000+1:5000],"out/datasets/allen10x_valid_lev1_dim1.h5",level=1L,chunkdim = c(nrow(e),1L))
  # write_h5_dataset(e[,15000+1:5000],"out/datasets/allen10x_test_lev1_dim1.h5",level=1L,chunkdim = c(nrow(e),1L))
}

make_shuffled_allen_datasets()
