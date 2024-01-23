
library(rhdf5)


#' Read rows of a count matrix from a connection to a CSV
#'
#' The first column should be a character with the rowname.
#' The CSV is expected to contain integers, and the function must be adapted for numeric values.
#'
#' @param con a connection to the CSV lines to reads
#' @param nrow desired number of row to read
#' @param ncol number of column (must be specified)
#' @return an integer matrix of size (nrow x ncol), possibly truncated if end of stream is reached. 
#'         The returned matrix has 0-rows if the function is further called after after end-of-stream.
scan_csv_matrix <- function(con,nrow,ncol) {
  m <- matrix(0L,nrow,ncol)
  rownames(m) <- seq_len(nrow)
  i <- 0L
  while(i<nrow) {
    row_name <- scan(con,character(),n = 1L,sep = ",",quiet = TRUE) # read the rowname
    if (length(row_name)==0) break # stop if rowname cannnot be read because of end-of-stream
    i <- i + 1L;
    rownames(m)[i] <- row_name
    
    n <- scan(con,integer(),n = ncol,sep = ",",nlines = 1L,quiet = TRUE)  # read the integer values of the row
    m[i,] <- n
  }
  m[seq_len(i),,drop=FALSE]
}


#' Helper function creating a compressed h5 dataset with unlimited rows
#' 
#' @param h5f an H5File object
#' @param name a length-1 character: name of the dataset to create
#' @param ncol number of column in the dataset
#' @param transpose a logical: TRUE if the dataset should be transposed into the H5
#' @param h5_type a character: type of elements stored into the dataset
h5_create_unlimited_row_dataset <- function(h5f,name,ncol,transpose=TRUE,h5_type="H5T_NATIVE_INT") {
  if (transpose) {
    fs <- H5Screate_simple(c(ncol,0L),maxdims = c(ncol,H5Sunlimited()))
  } else {
    fs <- H5Screate_simple(c(0L,ncol),maxdims = c(H5Sunlimited(),ncol))
  }
  prop <- H5Pcreate("H5P_DATASET_CREATE")
  H5Pset_chunk(prop,c(128L,128L))
  H5Pset_deflate(prop,6) # Set compression level
  h5x <- H5Dcreate(h5f,name,h5_type,fs,dcpl = prop) # Set value type to integers
  H5Pclose(prop)
  H5Sclose(fs)
  return(h5x)
}


#' Write rows of the matrix 'm' at the end of the given dataset
#' 
#' @param h5x the dataset object to update
#' @param row_index
#' @param ncol number of column in the dataset
#' @param transpose a logical: TRUE if the dataset should be transposed into the H5
#' @param h5_type a character: type of elements stored into the dataset
h5_rbind <- function(h5x,m,transpose=TRUE) {
  fs <- H5Dget_space(h5x)
  fsdim <- H5Sget_simple_extent_dims(fs)$size
  H5Sclose(fs)
  if (transpose) {
    stopifnot(ncol(m)==fsdim[1L])
    H5Dset_extent(h5x,c(fsdim[1L],fsdim[2L]+nrow(m)))
    fs <- H5Dget_space(h5x)
    m <- t(m)
    ms <- H5Screate_simple(dim(m))
    H5Sselect_hyperslab(fs,start=c(1L,fsdim[2L]+1L),count = dim(m))
  } else {
    stopifnot(ncol(m)==fsdim[2L])
    H5Dset_extent(h5x,c(fsdim[1L]+nrow(m),fsdim[2L]))
    fs <- H5Dget_space(h5x)
    ms <- H5Screate_simple(dim(m))
    H5Sselect_hyperslab(fs,start=c(fsdim[1L]+1L,1L),count = dim(m))
  }
  H5Dwrite(h5x,m,h5type=NULL,ms,fs)
  H5Sclose(ms)
  H5Sclose(fs)
}

#' Incrementally transform a CSV file containing a count matrix to H5 file
#' 
#' The script read the given csv file line by line and create H5 file with the matrix 
#' stored in the root object /x
#' @param con a connection to the CSV file
#' @pram h5file name of the h5file to create
h5_from_csv <- function(con,h5file) {
  open(con,"r")
  
  # read header line of the CSV (remove first element)
  col_names <- scan(con,sep = ",",nlines=1L,what = character())[-1]
  num_col <- length(col_names)
  
  # create h5 file
  h5f <- H5Fcreate(h5file)
  h5x <- h5_create_unlimited_row_dataset(h5f,"x",num_col)
  
  # iterate over rows of the CSV
  row_names <- character()
  while(TRUE) {
    m <- scan_csv_matrix(con,56L,num_col)
    if (nrow(m)==0L) break # stop if no data were read because of the end-of-stream
    h5_rbind(h5x,m)
    row_names <- c(row_names,rownames(m))
    print(length(row_names)) # display progression
  }
  H5Dclose(h5x)
  H5Fclose(h5f)
  close(con)
  h5write(col_names,h5file,"col_names")
  h5write(row_names,h5file,"row_names")
  
}

