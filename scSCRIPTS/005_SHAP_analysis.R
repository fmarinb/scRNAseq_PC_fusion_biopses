### PACKAGE LOADING ###

# SHAP value analysis
library(kernelshap)
library(shapviz)
library(tidyverse)
library(reshape2)

### SHAP VALUE CALCULATION WITH KERNELSHAP ###

optimal_lambda <- round(elasticnet_loocv$bestTune$lambda, 5)
print(optimal_lambda) # Knwon and save the optimal lambda of the Elastic Net model (max accuracy)


# SHAP values are calculated for all Elastic Net models generated when hyperparameter tuning was performed in caret. 
# Knowing the index and the value of the Elastic Net model with the optimal lambda, the set of SHAP values corresponding 
# to the model of interest is that of the 54th combination of hyperparameter tuning and the 54th combination of hyperparameter tuning.
shap_values <- kernelshap(elasticnet_loocv$finalModel, as.matrix(train_data_all[, -1])) # Calculate SHAP values.
opt_shap_values<-shap_values$s54 # Select SHAP values from the optimal model

# Plotting top 20 informative genes based on SHAP values: beeswarm plot (Fig 5B).
shap_imp<-sv_importance(opt_shap_values, kind = "bee", max_display = 20,
                        show_numbers = T) + # Add average SHAP value for each gene.
  theme(
    legend.position = c(1.05, 0.5), 
    legend.justification = c(0, 0.5), 
    plot.margin = unit(c(1, 3, 1, 1), "cm"))
sv_importance(opt_shap_values, show_numbers = TRUE)



#### GENE EXPRESSION PLOTS (FIG 5C-D) ###  

genes_of_interest <- c("PCA3", "ENSG00000289332", "TNFSF10", "CPNE4", "SCHLAP1") # Select the top 5 genes based of the change of the average SHAP value.
train_data_long <- melt(train_data_all, id.vars = "train_labels", 
                        measure.vars = genes_of_interest, 
                        variable.name = "Gene", value.name = "Expression") # Subset the train data (pseudobulk expression of the 52 genes found DE in DA analysis).

# Boxplot of the pseudobulk data
pseudobulk_markers<-ggplot(train_data_long, aes(x = train_labels, y = Expression, fill = train_labels)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +  # Remove outliers, they already appear as observations (no repeated points).
  geom_jitter(width = 0.2, alpha = 0.6) +          
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) + # Facet by gene.
  labs(x = "Sample type",
       y = "Normalized expression levels",
       fill = "Sample type") +
  theme_minimal() +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold", size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))



# Feature plots (UMAPs) showing the expression of the genes of interest at the single cell resolution.
PCA3<-FeaturePlot(SCP_data, 
                  features = "PCA3", 
                  reduction = "umap", 
                  cols = c("lightgrey", "blue")) +  
  xlab("UMAP 1") +  
  ylab("UMAP 2") +  
  ggtitle("PCA3") +  
  theme(axis.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        axis.text = element_text(size = 12)) +
  theme(legend.position = "right") 



ENSG00000289332 <- FeaturePlot(SCP_data, 
                               features = "ENSG00000289332", 
                               reduction = "umap", 
                               cols = c("lightgrey", "blue")) +  
  xlab("UMAP 1") +  
  ylab("UMAP 2") +  
  ggtitle("ENSG00000289332") +  
  theme(axis.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        axis.text = element_text(size = 12)) +
  theme(legend.position = "right") 


TNFSF10 <- FeaturePlot(SCP_data, 
                       features = "TNFSF10", 
                       reduction = "umap", 
                       cols = c("lightgrey", "blue")) +  
  xlab("UMAP 1") +  
  ylab("UMAP 2") +  
  ggtitle("TNFSF10") +  
  theme(axis.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        axis.text = element_text(size = 12)) +
  theme(legend.position = "right") 


CPNE4 <- FeaturePlot(SCP_data, 
                     features = "CPNE4", 
                     reduction = "umap", 
                     cols = c("lightgrey", "blue")) +  
  xlab("UMAP 1") +  
  ylab("UMAP 2") +  
  ggtitle("CPNE4") +  
  theme(axis.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        axis.text = element_text(size = 12)) +
  theme(legend.position = "right") 



SCHLAP1 <- FeaturePlot(SCP_data, 
                       features = "SCHLAP1", 
                       reduction = "umap", 
                       cols = c("lightgrey", "blue")) +  
  xlab("UMAP 1") +  
  ylab("UMAP 2") +  
  ggtitle("SCHLAP1") +  
  theme(axis.title = element_text(size = 15), 
        legend.text = element_text(size = 12), 
        axis.text = element_text(size = 12)) +
  theme(legend.position = "right") 






