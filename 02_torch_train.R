library(torch)
rm(list=ls())
source("lib_cell_classifier.R")
source("lib_torch_dataset.R")
source("lib_torch.R")
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Create datasets
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
make_datasets <- function() {
  x <- t(HDF5Array::HDF5Array("allen10x_all_lev1_dim1_NB.h5","/x"))
  rownames(x) <- as.character(rhdf5::h5read("allen10x_all_lev1_dim1_NB.h5","/row_names"))
  colnames(x) <- as.character(rhdf5::h5read("allen10x_all_lev1_dim1_NB.h5","/col_names"))
  
  y <- rhdf5::h5read("allen10x_all_lev1_dim1_NB.h5","/col_meta/class_label")
  y <- paste0(y,"|",rhdf5::h5read("allen10x_all_lev1_dim1_NB.h5","/col_meta/subclass_label"))
  y[y %in% c("","|")] <- NA
  y <- factor(y)

  # remove classes with low number of cell
  y <- factor(y,levels(y)[tabulate(y) >= 100])
  i <- which(!is.na(y))
  allen <- SummarizedExperiment(assays = t(x),
    colData = data.frame(class=y, row.names = rownames(x))
  )
    # Compute target labels to predict
  allen$y_true <- allen$class
  allen$y_true[allen$y_true %in% c("","|")] <- NA

  # Split cells into 3 datasets (trn/val/tst)
  set.seed(12345L)
  allen$dataset <- sample(factor(c("trn","val","tst"),c("trn","val","tst")),size = ncol(allen),replace = TRUE,prob = c(0.5,0.25,0.25))
  
  trn_ds <- allen[,allen$dataset %in% "trn"]
  trn_ds <- trn_ds[,!is.na(trn_ds$y_true)]
  
  
  val_ds <- allen[,allen$dataset %in% "val"]
  val_ds <- val_ds[,!is.na(val_ds$y_true)]

  trn_ds <- se_dataset(trn_ds,"y_true")
  val_ds <- se_dataset(val_ds,"y_true")
  
  saveRDS(trn_ds,file="trn_ds.rds")
  saveRDS(val_ds,file="val_ds.rds")
}
make_datasets()
trn_ds <- readRDS("trn_ds.rds")
val_ds <- readRDS("val_ds.rds")

DEVICE <- if (cuda_is_available()) "cuda" else "cpu"
# trn_ds$x <- torch_tensor(trn_ds$x,torch_float(),device=DEVICE,pin_memory=TRUE)
# trn_ds$y <- torch_tensor(trn_ds$y,torch_long(),device=DEVICE,pin_memory=TRUE)
# val_ds$x <- torch_tensor(val_ds$x,torch_float(),device=DEVICE,pin_memory=TRUE)
# val_ds$y <- torch_tensor(val_ds$y,torch_long(),device=DEVICE,pin_memory=TRUE)

trn_dl <- dataloader(trn_ds,batch_size=64L,num_workers=0L,worker_packages = c("HDF5Array","SummarizedExperiment"))
val_dl <- dataloader(val_ds,batch_size=64L,num_workers=0L,worker_packages = c("HDF5Array","SummarizedExperiment"))

model <- advanced_cell_classifier(rownames(trn_ds$o),levels(trn_ds$y_true))
fit1_log <- fit(model,trn_dl,val_dl,num_epoch = 50L,device=DEVICE)
torch_save(model,"out/Natalia/final_model_unpruned.torch")

invisible(model$compute_prune_mask())
fit2_log <- fit(model,trn_dl,val_dl,num_epoch = 50L,device=DEVICE)
torch_save(model,"out/Natalia/final_model_prune.torch")

