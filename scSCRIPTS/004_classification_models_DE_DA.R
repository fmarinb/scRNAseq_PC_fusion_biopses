### PACKAGE LOADING ###

# Datasets handling and additional plots.
library(tidyverse)
library(reshape2)

# Pseudobulk data building and normalization
library(Seurat)
library(DESeq2)

# Machine learning methods
library(caret)
library(devtools) 
library(RCurl)
library(glmnet) # Elastic Net
library(C50) # C 5.0 tree


#IMPORTANT: IF THE USER HAS EXCEUTED PREVIOUSLY THE COMANDS IN scATOMIC_clustering_DE_ssGSEA_(Fig1).R
# the steps of LOADING THE PRE-PROCESSED SEURAT FILE and ADD NEW METADATA INFO can be ommited

#### LOADING THE PRE-PROCESSED SEURAT FILE ###
  options(parallelly.fork.enable = TRUE) # Set up Seurat pararell computing.
  plan("multicore", workers = parallel::detectCores())
  SCP_data <- readRDS(file = "SC_Prostate_processed.rds") # Load seurat object with 0.4 resolution.
  SCP_data <- JoinLayers(SCP_data, assay = "RNA") # Join the per-run splitted layers.
  original_ident<-SCP_data$orig.ident
  rownm<-rownames(SCP_data@meta.data)
  SCP_data <- CalculateCDR(SCP_data) # Calculate Cell Detection Rate (CDR) with our own function (code available in a separate file).
  
  
  ### ADD NEW METADATA INFO ###
  # In this section we will add new columns at the metadata slot of the Seurat object.
  # Some informatiom are clinical and come from the sc_risk_scores.csv.
  sc_amico<-read.csv('sc_risk_scores.csv', sep = ";") # Clinical data loading.
  
  # We'll also create a new column clasifying cells in three main groups based on scATOMIC's annptation.
  metadata <- SCP_data@meta.data
  metadata$CAF_can_norm <- ifelse(metadata$layer_3 == "Cancer Associated Fibroblasts", "CAF", metadata$pan_cancer_cluster) # Renaming value to CAF.
  metadata<-metadata%>%
    mutate(CAF_can_norm = ifelse(CAF_can_norm == "Normal", "No Cancer", CAF_can_norm)) # New column CAF_can_norm divides cells into three groups.
  
  # Add new clinical information (ISUP, AMICO risk score and new sample naming).
  metadata <- metadata %>%
    left_join(sc_amico %>% select(Sample_Name_Group, 
                                  Sample_Name_Group2, 
                                  Sample_group2, ISUP, 
                                  Amico_risk_score), by = "Sample_Name_Group")
  
  SCP_data@meta.data <- metadata # Adding updated metadata in Seurat's object.
  rownames(SCP_data@meta.data)<-rownm
  SCP_data <- subset(SCP_data, pan_cancer_cluster %in% na.omit(SCP_data$pan_cancer_cluster)) # Remove cells with NA scATOMIC classification
  
  
  
  

#### TRAINNING DATAFRAME BUILDING ###

  # The pseudobulk data of the 52 DE genes from the comparisons of the DA cell neihbourhoods.
  # in Cancer and subclusters 1 and 3 will be used.
  
  genes_milo<-c(markers2$GeneID, markers_cancerct1$GeneID, markers_cancerct3$GeneID) # Get the names of the 52 genes.
  SCP_pseudobulk_milo<-AggregateExpression(SCP_data,
                                           assays = 'RNA', 
                                           return.seurat=F,
                                           features=genes_milo, 
                                           group.by = c("Sample_Name_Group2", 'Sample_group2')) # Create pseudobulk data by sample and sample class
  SCP_genes_milo<-as.data.frame(SCP_pseudobulk_milo)
  SCP_genes_milo<-as.matrix(SCP_genes_milo) # Transform into a matrix.
  
  # Normalization with DESeq2.
  sample_names <- colnames(SCP_genes_milo)
  conditions <- sapply(strsplit(sample_names, "_"), function(x) tail(x, 1)) # Sample class column
  sample_info <- data.frame(sample = sample_names,
                            condition = conditions)
  
  sample_info$condition <- factor(sample_info$condition)
  dds_milo <- DESeqDataSetFromMatrix(countData = SCP_genes_milo,
                                     colData = sample_info,
                                     design = ~ condition) # Desing matrix by condition (H or L samples)
  
  # Apply DESeq2's varianceStabilizingTranformation (log and log-derived data normalization)
  SCP_genes_norm_milo<- varianceStabilizingTransformation(dds_milo, blind = TRUE)
  SCP_genes_norm_matrix_milo <- assay(SCP_genes_norm_milo)
  SCP_genes_input_milo<-t(SCP_genes_norm_matrix_milo) # Transpose matrix
  
  # Get final train data and conditions
  set.seed(123)
  train_data <- SCP_genes_input_milo
  train_labels <- sample_info$condition
  train_data_all<-cbind(train_labels, train_data)%>%as.data.frame()
  train_data_all$train_labels<-as.factor(train_data_all$train_labels)
  levels(train_data_all$train_labels) <- c("H", "L") # Add final condition labels (debug from sample renaming)
  

### TRAIN ML MODELS WITH LOOCV (leave-one-out Cross Validation) ###

  set.seed(123)  # Ensures reproducibility of results
  
  # Define the training control settings using Leave-One-Out Cross-Validation (LOOCV)
  trct_loocv <- trainControl(
    method = "LOOCV",       # Specifies Leave-One-Out Cross-Validation (each observation is a validation set once)
    allowParallel = TRUE,   # Enables parallel computation to speed up training
    verboseIter = TRUE,     # Displays progress messages during training
    classProbs = TRUE,      # Computes class probabilities (important for classification problems)
    savePredictions = 'final',  # Saves the final predictions after training
    summaryFunction = twoClassSummary  # Uses a function to compute performance metrics for binary classification
  )
  
  set.seed(123)
  # Train a Support Vector Machine (SVM) model with a radial basis function (RBF) kernel
  svm_loocv <- train(
    train_labels ~ .,       # Formula: Predict train_labels using all available predictors
    data = train_data_all,  # The dataset used for training
    method = "svmRadial",   # Specifies the SVM with an RBF kernel
    trControl = trct_loocv, # Uses the previously defined LOOCV control settings
    metric = "ROC",         # Optimizes the model based on the ROC AUC score
    tuneLength = 10         # Specifies the number of values to try when tuning hyperparameters
  )
  
  set.seed(123)
  # Train a Random Forest (RF) model
  rf_loocv <- train(
    train_labels ~ .,       
    data = train_data_all,  
    method = "rf",          # Specifies Random Forest as the model
    trControl = trct_loocv, 
    metric = "ROC", 
    tuneLength = 10         
  )
  
  set.seed(123)
  # Train an Elastic Net model (combines Lasso and Ridge regression)
  elasticnet_loocv <- train(
    train_labels ~ .,       
    data = train_data_all,  
    method = "glmnet",      # Specifies Generalized Linear Model with Elastic Net regularization
    trControl = trct_loocv, 
    metric = "ROC", 
    tuneLength = 10         
  )
  
  set.seed(123)
  # Train a C5.0 Decision Tree model
  c50_loocv <- train(
    train_labels ~ .,       
    data = train_data_all,  
    method = "C5.0",        # Specifies the C5.0 decision tree algorithm
    trControl = trct_loocv, 
    metric = "ROC", 
    tuneLength = 10         
  )
  
  
  
  
  ### TRAIN ML MODELS WITH BOOTSTRAP ### 
  
  set.seed(123)  # Ensures reproducibility of results
  
  # Define the training control settings using Bootstrapping
  boot_control <- trainControl(
    method = "boot",       # Specifies Bootstrapping as the resampling method
    number = 100,          # Number of bootstrap resamples (default is often 25, but here we use 100)
    verboseIter = TRUE,    # Displays progress messages during training
    classProbs = TRUE,     # Computes class probabilities (useful for classification tasks)
    savePredictions = 'final',  # Saves final model predictions
    returnResamp = 'final',  # Stores the final resampling results
    summaryFunction = twoClassSummary  # Uses a function to compute performance metrics for binary classification
  )
  
  set.seed(123)
  # Train a Support Vector Machine (SVM) model with a radial basis function (RBF) kernel
  svm_boot <- train(
    train_labels ~ .,       # Formula: Predict train_labels using all available predictors
    data = train_data_all,  # The dataset used for training
    method = "svmRadial",   # Specifies the SVM with an RBF kernel
    trControl = boot_control, # Uses the previously defined Bootstrapping control settings
    metric = "ROC",         # Optimizes the model based on the ROC AUC score
    tuneLength = 10         # Specifies the number of values to try when tuning hyperparameters
  )
  
  set.seed(123)
  # Train a Random Forest (RF) model
  rf_boot <- train(
    train_labels ~ .,       
    data = train_data_all,  
    method = "rf",          # Specifies Random Forest as the model
    trControl = boot_control, 
    metric = "ROC", 
    tuneLength = 10         
  )
  
  set.seed(123)
  # Train an Elastic Net model (combines Lasso and Ridge regression)
  elasticnet_boot <- train(
    train_labels ~ .,       
    data = train_data_all,  
    method = "glmnet",      # Specifies Generalized Linear Model with Elastic Net regularization
    trControl = boot_control, 
    metric = "ROC", 
    tuneLength = 10         
  )
  
  set.seed(123)
  # Train a C5.0 Decision Tree model
  c50_boot <- train(
    train_labels ~ .,       
    data = train_data_all,  
    method = "C5.0",        # Specifies the C5.0 decision tree algorithm
    trControl = boot_control, 
    metric = "ROC", 
    tuneLength = 10         
  )  


 ### COMPUTE PERFORMANCE METRICS (LOOCV) ###
  
# We have created an acapz function to extract sensitivity and specificity from the confusion matrix 
# generated for the accuracy calculation (value matched when the metric argument of train() is set to ‘Accuracy’).
# The same matrix is used to calculate the F1 score from the extraction of precision and recall. This function also 
# allows to calculate from the predictions and observations the AUC from ROC, as well as the G-mean.
  
  calculate_loocv_metrics <- function(model) {
    # Verify that the model contains predictions
    if (is.null(model$pred)) {
      stop("The model does not contain predictions. Make sure to use savePredictions = 'all' in trainControl.")
    }
    
    # Initialize lists to store metrics
    
    # Create confusion matrix
    conf_matrix <- confusionMatrix(model$pred$pred, model$pred$obs, positive='T')
    
    # Extract metrics
    sens <- conf_matrix$byClass["Sensitivity"]  # Sensitivity
    spec <- conf_matrix$byClass["Specificity"]  # Specificity
    acc  <- conf_matrix$overall["Accuracy"]     # Accuracy
    
    # Manually calculate F1-Score
    precision <- conf_matrix$byClass["Pos Pred Value"]  # Precision
    recall <- sens  # Recall (same as Sensitivity)
    f1 <- ifelse((precision + recall) > 0, 
                 2 * (precision * recall) / (precision + recall), 
                 NA)
    
    # Compute ROC curve and AUC
    roc_curve <- roc(model$pred$obs, model$pred[['T']], 
                     levels = levels(model$pred$obs))
    auc <- auc(roc_curve)
    
    # Compute G-Mean
    gmean <- sqrt(sens * spec)
    
    # Create results list
    metrics <- c(
      Sensitivity = round(sens, 4),
      Specificity = round(spec, 4),
      Accuracy    = round(acc, 4),
      F1_Score    = round(f1, 4),
      G_Mean      = round(gmean, 4),
      AUC         = auc
    )
    
    # Return metrics
    return(metrics)
  }


  # Create a data frame to store de performance metric values of each model
  svm_metrics_loocv <- as.data.frame(calculate_loocv_metrics(svm_loocv))
  rf_metrics_loocv <- as.data.frame(calculate_loocv_metrics(rf_loocv))
  elasticnet_metrics_loocv <- as.data.frame(calculate_loocv_metrics(elasticnet_loocv))
  c50_metrics_loocv <- as.data.frame(calculate_loocv_metrics(c50_loocv))
  
  metrics_loocv<-cbind(SVM_LOOCV=svm_metrics_loocv,
                       RF_LOOCV=rf_metrics_loocv,
                       GLMNET_LOOCV=elasticnet_metrics_loocv,
                       C5.0_LOOCV=c50_metrics_loocv)



### COMPUTE PERFORMANCE METRICS (BOOTSTRAP) ###

ccalculate_boot_metrics <- function(model) {
  # Verify that the model contains predictions
  if (is.null(model$pred)) {
    stop("The model does not contain predictions. Make sure to use savePredictions = 'all' in trainControl.")
  }
  
  # Initialize lists to store metrics
  sens_list <- c()
  spec_list <- c()
  acc_list  <- c()
  f1_list   <- c()
  gmean_list <- c()
  auc_list <- c()
  
  # Loop through each Resample
  for (resample in unique(model$pred$Resample)) {
    # Filter predictions for the current iteration
    resample_preds <- model$pred %>% filter(Resample == resample)
    
    # Create confusion matrix
    conf_matrix <- confusionMatrix(resample_preds$pred, resample_preds$obs)
    
    # Extract metrics
    sens <- conf_matrix$byClass["Sensitivity"]  # Sensitivity
    spec <- conf_matrix$byClass["Specificity"]  # Specificity
    acc  <- conf_matrix$overall["Accuracy"]     # Accuracy
    
    # Manually calculate F1-Score
    precision <- conf_matrix$byClass["Pos Pred Value"]  # Precision
    recall <- sens  # Recall (same as Sensitivity)
    f1 <- ifelse((precision + recall) > 0, 
                 2 * (precision * recall) / (precision + recall), 
                 NA)
    
    # Compute ROC curve and AUC
    roc_curve <- roc(model$pred$obs, model$pred[['T']], 
                     levels = levels(model$pred$obs))
    auc <- auc(roc_curve)
    
    # Compute G-Mean
    gmean <- sqrt(sens * spec)
    
    # Store metrics
    sens_list <- c(sens_list, sens)
    spec_list <- c(spec_list, spec)
    acc_list  <- c(acc_list, acc)
    f1_list   <- c(f1_list, f1)
    gmean_list <- c(gmean_list, gmean)
    auc_list <- c(auc_list, auc)
  }
  
  # Compute average metrics
  mean_sens <- mean(sens_list, na.rm = TRUE)
  mean_spec <- mean(spec_list, na.rm = TRUE)
  mean_acc  <- mean(acc_list, na.rm = TRUE)
  mean_f1   <- mean(f1_list, na.rm = TRUE)
  mean_gmean <- mean(gmean_list, na.rm = TRUE)
  mean_auc <- mean(auc_list, na.rm = TRUE)
  
  # Create results list
  metrics <- c(
    Sensitivity = round(mean_sens, 4),
    Specificity = round(mean_spec, 4),
    Accuracy    = round(mean_acc, 4),
    F1_Score    = round(mean_f1, 4),
    G_Mean      = round(mean_gmean, 4),
    AUC         = round(mean_auc, 4)
  )
  
  # Return metrics
  return(metrics)
}

  
  # Create a data frame to store de performance metric values of each model
  svm_metrics_boot <- as.data.frame(calculate_boot_metrics(svm_boot))
  rf_metrics_boot <- as.data.frame(calculate_boot_metrics(rf_boot))
  elasticnet_metrics_boot <- as.data.frame(calculate_boot_metrics(elasticnet_boot))
  c50_metrics_boot <- as.data.frame(calculate_boot_metrics(c50_loocv))
  
  
  metrics_boot<-cbind(SVM_boot=svm_metrics_boot,
                      RF_boot=rf_metrics_boot,
                      GLMNET_boot=elasticnet_metrics_boot,
                      C5.0_boot=c50_metrics_boot)
  
  metrics_allmodels<-rbind(t(metrics_loocv), t(metrics_boot))
  
  
  metrics_plot<-read.csv('metrics_plot.csv', sep=';', stringsAsFactors = T)
  
  
  # The tables with values of the performance metrics of the models by LOOCV and bootstrap were
  # downloaded and spliced to form the table metrics_plot.csv. This file was reloaded into R for graphical representation (Figure 5a).
  
  # Reshape the 'metrics_plot' dataframe from wide format to long format 
  metrics_long <- melt(metrics_plot, 
                       id.vars = c("Model", "Validation"),  # Variables to keep unchanged
                       measure.vars = c("Accuracy", "F1.Score", "G.Mean", "AUC"))  # Metrics to be transformed into a single column

  
  # Create the dotplot dividing by validation method
  ggplot(metrics_long, aes(x = variable, y = value, color = Model)) +
    geom_point(size = 3) +
    facet_grid(. ~ Validation, ) + # facet dotplot by LOOCV and Bootstrap
    labs(
      x = "Metrics",
      y = "Values",
      color = "Model") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.title.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          strip.text = element_text(size = 14, face = "bold"))
