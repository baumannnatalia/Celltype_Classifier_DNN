library(torch)
library(matrixStats)
source("lib_torch_dataset.R")
batch_predict <- function(m,x) {
  dl <- dataloader(x,batch_size=1024L,num_workers=0L,worker_packages = c("HDF5Array","SummarizedExperiment"))
  g <- coro::generator(function() for (o in dl) {print(".");yield(m(o$x))})
  P <- coro::collect(g())
  P <- torch_cat(P)
  P <- as_array(P)
  colnames(P) <- m$class_names
  P
}

predict_type_AB <- function(dataset){
  model <- torch_load("final_model_prune.torch")
  print(paste0(sum(model$feature_names%in%rownames(dataset))," Genes used for prediction out of ",length(model$feature_names)))
  x <- SummarizedExperiment(
    rbind(0,dataset@assays$RNA@counts)[1L+match(model$feature_names,rownames(dataset),nomatch = 0L),],
    colData = dataset@meta.data
  )
  assay(x) <- assay(x)/colSums(assay(x))[col(x)]*1e4
  dataset@reductions$Prediction_AB <- batch_predict(model,se_dataset(x))
  dataset$AB.cluster.predicted <- colnames(dataset@reductions$Prediction_AB)[max.col(dataset@reductions$Prediction_AB)]
  dataset$AB.prediction.confidence <- apply(dataset@reductions$Prediction_AB,1, function(l){
    l[order(l, decreasing = T)[1]]-l[order(l, decreasing = T)[2]]
  })
  
  CGE_derived <- c("GABAergic|Lamp5","GABAergic|Vip","GABAergic|Sncg")
  MGE_derived <- c("GABAergic|Pvalb","GABAergic|Sst","GABAergic|Sst Chodl")
  
  
  Lamp5 <- c("GABAergic|Lamp5","GABAergic|Vip","GABAergic|Sncg")
  MGE_derived <- c("GABAergic|Pvalb","GABAergic|Sst","GABAergic|Sst Chodl")
  Meis2_IN <- "GABAergic|Meis2"
  CR <- "Glutamatergic|CR"
  L6 <- c("Glutamatergic|L6 CT CTX","Glutamatergic|L6 IT CTX","Glutamatergic|L6b CTX","Glutamatergic|Car3")
  L5 <- c("Glutamatergic|L5 IT CTX","Glutamatergic|L5 PT CTX")
  L5_6_NP <- c("Glutamatergic|L5/6 NP CTX")
  L4_5 <- c("Glutamatergic|L4/5 IT CTX")
  L4 <- c("Glutamatergic|L4 IT CTX")
  L2_3 <- c("Glutamatergic|L2/3 IT CTX","Glutamatergic|L2/3 IT PPP")
  Astro <- c("Non-Neuronal|Astro")
  Endothelial <- ("Non-Neuronal|Endo")
  Oligo <- c("Non-Neuronal|Oligo")
  Others <- c("Non-Neuronal|Micro-PVM","Non-Neuronal|SMC-Peri","Non-Neuronal|VLMC")
  types <- list(CGE_derived=CGE_derived,MGE_derived=MGE_derived,
                L6=L6,L5_6_NP=L5_6_NP,L5=L5,L4_5=L4_5,L4=L4,L2_3=L2_3,
                Astro=Astro,Endothelial=Endothelial,Oligo=Oligo)
  
  dataset$AB.H3 <- NA
  dataset$AB.H3[dataset$AB.cluster.predicted%in%CGE_derived] <- "CGE-derived"
  dataset$AB.H3[dataset$AB.cluster.predicted%in%MGE_derived] <- "MGE-derived"
  dataset$AB.H3[dataset$AB.cluster.predicted%in%L6] <- "L6"
  dataset$AB.H3[dataset$AB.cluster.predicted%in%L5] <- "L5"
  dataset$AB.H3[dataset$AB.cluster.predicted%in%L4_5] <- "L4_5"
  dataset$AB.H3[dataset$AB.cluster.predicted%in%L5_6_NP] <- "L5_6_NP"
  dataset$AB.H3[dataset$AB.cluster.predicted%in%L2_3] <- "L2-3"
  dataset$AB.H3[dataset$AB.cluster.predicted%in%L4] <- "L4"
  dataset$AB.H3[dataset$AB.cluster.predicted%in%Astro] <- "Astro"
  dataset$AB.H3[dataset$AB.cluster.predicted%in%Endothelial] <- "Endothelial"
  dataset$AB.H3[dataset$AB.cluster.predicted%in%Oligo] <- "Oligo"
  dataset$AB.H3[dataset$AB.cluster.predicted%in%Others] <- "Other"
  dataset$AB.H3[dataset$AB.cluster.predicted%in%Meis2_IN] <- "Meis2_IN"
  
  dataset$AB.H1 <- NA
  dataset$AB.H1[dataset$AB.H3%in%c("Oligo","Endothelial","Astro","Other")] <- "Non-Neuronal"
  dataset$AB.H1[dataset$AB.H3%in%c("L6","L2-3","L4","L5", "L4_5","L5_6_NP")] <- "Glutamatergic"
  dataset$AB.H1[dataset$AB.H3%in%c("CGE-derived","MGE-derived","Meis2_IN")] <- "GABAergic"
  
  # dataset$AB.H1 <- "Neurons"
  # dataset$AB.H1[dataset$AB.H3%in%c("Oligo","Endothelial","Astro","Other")] <- "Non-Neuronal"
  
  
  reduced_table <- lapply(types,function(i){
    if(length(i)>1){
    rowMaxs(dataset@reductions$Prediction_AB[,i])
    }else{
      dataset@reductions$Prediction_AB[,i]
    }
  })
  reduced_table <- do.call(cbind, reduced_table)
  
  dataset$AB.H3.prediction.confidence <- apply(reduced_table,1, function(l){
    l[order(l, decreasing = T)[1]]-l[order(l, decreasing = T)[2]]
  })
  
  
  dataset
}
