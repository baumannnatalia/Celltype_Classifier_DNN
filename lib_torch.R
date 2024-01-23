library(torch) #options(timeout=300);torch::install_torch()
library(ggplot2)
library(R6)



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Helper class storing batch learning statistics
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
batch_logger <- R6Class("batch_logger", list(
  name = "",      # name of the logger (only used for progression display)
  num_batch = 0L, # number of batch that will be processed (only used for progression display)
  losses = numeric(0),
  y_pred = list(),
  y_true = list(),
  start_time = Sys.time(),
  initialize = function(name="",num_batch=0L) {self$name <- name;self$num_batch <- num_batch;self$losses<-numeric(0);self$y_pred<-self$y_true<-list();self$start_time = Sys.time()},
  loss = function() {mean(self$losses)},
  accuracy = function() {
    y_true <- unlist(self$y_true)
    y_pred <- unlist(self$y_pred)
    mean(y_pred == y_true)
  },
  contingency_matrix  = function() {
    y_true <- unlist(self$y_true)
    y_pred <- unlist(self$y_pred)
    table(y_true,y_pred)
  },
  add = function(loss,y_pred,y_true=NULL) {
    self$losses[[length(self$losses)+1L]] <- loss
    self$y_pred[[length(self$y_pred)+1L]] <- y_pred
    self$y_true[[length(self$y_true)+1L]] <- y_true
    #cat(sprintf("\r%s:%d/%d", self$name, length(self$losses), self$num_batch))
    cat(sprintf("\r%s:%d/%d batch_loss: %1.5f; avg_loss: %1.5f; accuracy: %.1f%%; %.3fs/batch", self$name, length(self$losses), self$num_batch, loss, self$loss(), 100*self$accuracy(),difftime(Sys.time(),self$start_time,units="secs")/length(self$losses)))
  },
  plot = function() {
    m <- self$contingency_matrix()
    m <- data.frame(y_true=factor(rownames(m)[row(m)],rownames(m)),y_pred=factor(colnames(m)[col(m)],colnames(m)),n=as.vector(m))
    ggplot(m,aes(x=y_pred,y=y_true,fill=n)) + 
      geom_tile(aes(fill=n)) +
      geom_text(aes(label=n)) +
      ggtitle(sprintf("accuracy: %.1f%%",100*self$accuracy())) +
      scale_fill_viridis_c() +
      coord_equal()
  }
))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Helper class storing epoch learning statistics
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
epoch_logger <- R6Class("batch_logger", list(
  epoch = 1L,
  metrics = data.frame(epoch=integer(0),metric=character(0),value=numeric(0),group=character(0)),
  initialize = function() {},
  end_of_epoch = function() {self$epoch <- self$epoch + 1L},
  add = function(metric,value,group="group") {
    e <- data.frame(epoch=self$epoch,metric,value,group)
    self$metrics <- rbind(self$metrics,e)
  },
  plot = function() {
    p <- ggplot(self$metrics) + facet_wrap(~metric,scales = "free_y") +
      geom_line(aes(x=epoch,y=value,color=group)) +
      ylab("value")
    print(p)
  }
))





# Fit a model on given training data
fit <- function(model,trn_dl,val_dl=NULL,num_epoch=30L,device='cpu',loss=nnf_multi_margin_loss) {
  model$to(device=device)
  
  log <- epoch_logger$new()  
  optimizer <- optim_adam(model$parameters,lr=0.0001,weight_decay = 0)
  for (epoch in seq(num_epoch)) {
    #-#-#-#-#-#-#-#-#-#-#
    # Epoch training
    #-#-#-#-#-#-#-#-#-#-#
    model$train() # activate train mode
    trn_log <- batch_logger$new(sprintf("Epoch %d:trn",log$epoch))
    coro::loop(for (b in trn_dl) {
      optimizer$zero_grad()
      y_pred <- model$forward(b$x$to(device=device))
      loss_val <- loss(y_pred, b$y$to(device=device))
      loss_val$backward()   # have gradients get calculated
      optimizer$step()  # have gradients get applied
      trn_log$add(loss_val$item(),as_array(y_pred$argmax(2L)$to(device="cpu")),as_array(b$y$to(device="cpu")))
    })
    log$add("loss",mean(trn_log$losses),"trn")
    log$add("acc",trn_log$accuracy(),"trn")
    cat("\n")
    gc()
    
    #-#-#-#-#-#-#-#-#-#-#
    # Epoch validation
    #-#-#-#-#-#-#-#-#-#-#
    model$eval() # activate eval mode
    val_log <- batch_logger$new(sprintf("Epoch %d:val",log$epoch))
    if (!is.null(val_dl)) {
      coro::loop(for (b in val_dl) {
        y_pred <- model$forward(b$x$to(device=device))
        loss_val <- loss(y_pred, b$y$to(device=device))
        val_log$add(loss_val$item(),as_array(y_pred$argmax(2L)$to(device="cpu")),as_array(b$y$to(device="cpu")))
      })
      log$add("loss",mean(val_log$losses),"val")
      log$add("acc",val_log$accuracy(),"val")
      cat("\n")
    }
    
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Epoch end
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    #torch_save(model,"out/last_model.torch")
    #saveRDS(log,file="out/last_log.rds")
    #plot(log)
    log$end_of_epoch()
  }
  list(log=log,trn_log=trn_log,val_log=val_log)
}




