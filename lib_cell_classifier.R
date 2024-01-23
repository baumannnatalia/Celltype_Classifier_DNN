



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Define a basic model
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
basic_cell_classifier <- nn_module(
  "basic_cell_classifier",
  
  initialize = function(feature_names,class_names) {
    self$feature_names <- feature_names
    self$class_names <- class_names
    self$layers <- nn_sequential(
      nn_dropout(0.5), nn_linear(in_features = length(self$feature_names), out_features = 1024L,bias = FALSE), nn_hardtanh(),
      nn_dropout(0.3), nn_linear(in_features = 1024L, out_features = 512L), nn_relu(),
      nn_dropout(0.3), nn_linear(in_features = 512L, out_features = 256L), nn_relu(),
      nn_dropout(0.3), nn_linear(in_features = 256L, out_features = length(self$class_names))
    )
  },
  
  forward = function(rpm) {
      rpm$log1p() %>% self$layers()
  }
)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Define an advanced model with layer1 scaling and pruning
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
advanced_cell_classifier <- nn_module(
  "advanced_cell_classifier",
  
  initialize = function(feature_names,class_names) {
    self$feature_names <- feature_names
    self$class_names <- class_names

    self$layer1_dropout <- nn_dropout(0.5)
    self$layer1_linear <- nn_linear(in_features = length(self$feature_names), out_features = 1024L,bias = FALSE)
    self$prune_mask <- nn_buffer(torch_ones_like(self$layer1_linear$weight,torch_bool()))
    
    self$layer2 <- nn_sequential(
      nn_hardtanh(),
      nn_dropout(0.3), nn_linear(in_features = 1024L, out_features = 512L), nn_relu(),
      nn_dropout(0.3), nn_linear(in_features = 512L, out_features = 256L), nn_relu(),
      nn_dropout(0.3), nn_linear(in_features = 256L, out_features = length(self$class_names))
    )
  },
  
  compute_prune_mask = function(n=50L) {
    lo <- self$layer1_linear$weight$topk(n,largest = FALSE)[[1]][,n,drop=FALSE]
    hi <- self$layer1_linear$weight$topk(n,largest = TRUE)[[1]][,n,drop=FALSE]
    self$prune_mask$set_((self$layer1_linear$weight<=lo) | (self$layer1_linear$weight>=hi))
    NULL
  },
  
  forward = function(rpm) {
    # apply pruning mask and normalize weights to 1
    with_no_grad({
      self$layer1_linear$weight$mul_(self$prune_mask)
      self$layer1_linear$weight$div_(self$layer1_linear$weight$norm(dim=2L,keepdim=TRUE))
    })
    
    # preprocessing
    log_rpm <- rpm$log1p()

    # compute normalization factors for each node of layer1 and each instance
    sf <- torch_matmul(log_rpm,self$layer1_linear$weight$abs()$transpose(1,2))
    sf <- sf$clamp(0.01)
    
    log_rpm %>% 
      self$layer1_dropout() %>%
      self$layer1_linear() %>%
      torch_div(sf) %>%
      self$layer2()
  }
)



