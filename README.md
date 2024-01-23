# Celltype_Classifier_DNN
Scripts to train the feedforward neural network to predict cell type identity from Seurat object.

The original dataset used to train the model is available from the Allen Brain institute website : https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x

The scripts 00_make_allen_hdf5.R and 01_make_ml_dataset.R are used to create dataset ready to use it for training.

02_torch_train.R creates train and valid datasets and trains the model. 

Prediction_types_ABA_SSp.R contains the function to predict cell identity from the model (final_modelprune.toch) of a new seurat object. Be careful of having the "RNA" assay as default in your seurat object.

