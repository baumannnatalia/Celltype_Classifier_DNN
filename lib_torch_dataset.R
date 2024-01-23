
library(torch)
library(SummarizedExperiment)



#' Random sampling of a given number of instance per class (with replacement)
#' @param y a class label factor
#' @param n number of instance per class to 
#' @return an randomized integer vector
#' NA are ignored
balanced_bootstrap <- function(y,n) {
  i <- split(seq_along(y),y)
  bsi <- lapply(i,sample,size=n,replace=TRUE)
  sample(unlist(bsi,use.names = FALSE))
}



se_dataset <- dataset(
  "SummarizedExperimentDataset",
  initialize = function(o,y_true=NULL,assay_index=1L) {
    # check parameters
    stopifnot(is(o,"SummarizedExperiment"))
    if (is.character(y_true) && (length(y_true)==1L)) {
      y_true <- o[[y_true]]
    } else if (is.null(y_true)) {
      y_true <- gl(1,ncol(o),labels = "?")
    }
    stopifnot(is.factor(y_true) && (length(y_true)==ncol(o)))

    self$o <- o
    self$assay_index <- assay_index
    self$y_true <- y_true
  },
  
  .getbatch = function(j) {
    x <- as.array(t(assay(self$o,self$assay_index)[,j,drop=FALSE]))
    x[is.na(x)] <- 0
    x <- torch_tensor(x,torch_float())
    
    # RPM
    x <- nnf_normalize(x,p = 1,dim=2L,eps = 100)$mul(10000)
    
    y <- torch_tensor(self$y_true[j],torch_long())
    list(x=x,y=y)
  },
  .length = function() {ncol(self$o)}
)



