### PACKAGE LOADING ###

# Datasets handling and additional plots.
library(dplyr)
library(data.table)
library(parallel)
library(reticulate)
library(ggplot2)
library(plyr)
library(tidyverse)
library(Matrix)
library(patchwork)

# Differential abundance analysis.
library(miloR)
library(scran)

# Single cell analysis + differential expression analysis.
library(Seurat)
library(SingleCellExperiment)
library(scater)


# IMPORTANT: IF THE USER HAS EXCEUTED PREVIOUSLY THE COMANDS IN scATOMIC_clustering_DE_ssGSEA_(Fig1).R
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


### DIFFERENTIAL ABUNDANCE ANALYSIS ###
  
  # We'll follow the rpocedure specified in https://www.bioconductor.org/packages/release/bioc/vignettes/miloR/inst/doc/milo_demo.html
  # and  https://www.bioconductor.org/packages/release/bioc/vignettes/miloR/inst/doc/milo_gastrulation.html
  SCP_data_sce<-as.SingleCellExperiment(SCP_data)# Transform Seurat data into SingleCellExperiment.
  SCP_data_milo<-Milo(SCP_data_sce) # Run Milo with default parameters.


  set.seed(123)
  SCP_data_milo<-buildGraph(SCP_data_milo, # Milo's implementation of buildKNNGraph function in scran.
                            k=30, # Number of nearest-neighbours to consider for the graph building.
                            d=50) # The number of dimensions to use if the input is a matrix of cells X reduced dimensions.
  
  SCP_data_milo<-makeNhoods(SCP_data_milo, 
                            prop = 0.1, # A double scalar that defines what proportion of graph vertices to randomly sample (recommended value: 0.1).
                            k =30, 
                            d=50, 
                            refined = TRUE, # Refined sampling scheme,
                            refinement_scheme="graph") # Sample scheme.
  
  SCP_data_milo<-countCells(SCP_data_milo, # Quantifies the number of cells in each neighbourhood.
                            meta.data = data.frame(colData(SCP_data_milo)),
                            samples = 'Sample_Name_Group2')

  plotNhoodSizeHist(SCP_data_milo) # Histogram of the number of cells belonging to each neighbourhood.
  head(nhoodCounts(SCP_data_milo))




  set.seed(123)
  SCP_data_milo$Sample_group2<-as.factor(SCP_data_milo$Sample_group2)
  SCP_data_milo$ISUP<-as.factor(SCP_data_milo$ISUP)
  nhood_counts <- nhoodCounts(SCP_data_milo) # Returns a NxM sparse matrix of cell counts in each of N neighbourhoods with respect to the M experimental samples defined.
  
  # Prepare the design matrix for differencial abundance testing
  SCP_milo_design<-data.frame(colData(SCP_data_milo))[,c('Sample_Name_Group2', 
                                                         'Sample_group2', 
                                                         'ISUP')]
  SCP_milo_design <- distinct(SCP_milo_design)
  rownames(SCP_milo_design) <- SCP_milo_design$Sample_Name_Group2
  SCP_milo_design  <-SCP_milo_design[colnames(nhood_counts), , drop=FALSE]
  
  # Differential abundance testing (edgeR method)
  contrast<-c("Sample_group2H - Sample_group2L") # Specify the contrats H (ref) vs L samples.
  da_results <- testNhoods(SCP_data_milo, 
                           design = ~ 0 + Sample_group2,# Formula or model.matrix object describing the experimental design for differential abundance testing.
                           design.df = SCP_milo_design, model.contrasts = contrast,
                           fdr.weighting="graph-overlap", # The spatial FDR weighting scheme to use.
                           norm.method="TMM") # Trimmed mean of M-values.
  
  da_results_significative<-da_results %>%
    arrange(SpatialFDR) %>%
    filter(SpatialFDR < 0.1)
  dim(da_results_significative) # Check how many cell neighbourhoods are DA between T and L samples.


### PLOTTING BEESWARM PLOTS OF DIFF. ABUNDANCE ANALYSIS ###
 # Beeswarm plot of DA cell neigbourhoods using scATOMIC subpopulations layer.
  da_results <- annotateNhoods(SCP_data_milo, da_results, coldata_col = "scATOMIC_pred") # Annotate cell neihgbourhoods depending on its scATOMIC cell group.
  plotDAbeeswarm(da_results,
                 group.by = "scATOMIC_pred", 
                 alpha = 0.1)+ # significance level for Spatial FDR (default: 0.1).
    geom_point(size = 3) +
    theme( 
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14))
  
  # Beeswarm plot of DA cell neigbourhoods using three class layer (Cancer, CAF, No Cancer).
  # If there are cell neighbourhoods composed of cells belongig to two different groups with over 70% of sharing, the cell group is categorised as mixed.
  da_results <- annotateNhoods(SCP_data_milo, da_results, coldata_col = "CAF_can_norm")
  da_results$CAF_can_norm <- ifelse(da_results$CAF_can_norm < 0.7, "Mixed", da_results$CAF_can_norm)
  
  plotDAbeeswarm(da_results, 
                    group.by = "CAF_can_norm",
                    alpha = 0.1) + 
    geom_point(size = 3) +
    theme( 
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14))
  

  
### DIFFERENTIAL EXPRESSION ANALYSIS OF GROUPED DA NEIHBOURHOODS ###

  # Grouping the DA cancer cell neighbourhoods into a unique neighbourhood group called 'Cancer DA'.
  da_results <- da_results %>%
    mutate(NhoodGroup = case_when(
      SpatialFDR < 0.1 & CAF_can_norm == "Cancer" ~ 'Cancer DA',
      SpatialFDR > 0.1 & CAF_can_norm == "Cancer" ~ 'Cancer',
      TRUE ~ 'Normal + CAF'
    ))
  
  
  
  keep.rows <- rowSums(logcounts(SCP_data_milo)) != 0
  SCP_data_milo <- SCP_data_milo[keep.rows, ] # Filtering out genes with 0 counts in all samples.
  
  
  dec <- modelGeneVar(SCP_data_milo)
  hvgs <- getTopHVGs(dec) # Find HVGs.
  
  
  # Perform differential expression testing betweem DA Cancer and non Cancer cell neihbourhoods group .
  test2 <- findNhoodGroupMarkers( # we use the logcounts data slot (default).
    SCP_data_milo, 
    da_results, 
    subset.row = hvgs, # Use variable genes for sumamrizing over cells in neighbourhoods.
    aggregate.samples = TRUE, # cells in the same sample and neighbourhood group should be merged for DGE testing.
    subset.nhoods = da_results$NhoodGroup %in% c("Cancer", "Cancer DA"), # Compare Cancer vs Cancer DA.
    sample_col = "Sample_Name_Group2")

  markers2<-test2%>%filter(., `adj.P.Val_Cancer DA`<0.05 & abs(`logFC_Cancer DA`) > 1.5 ) # Get the resulting overexpressed DE genes.

  
### DIFFERENTIAL ABUNDANCE ANALYSIS GROUPPING BY HIHGH/LOW ISUP (NO DA NEIGHBOURHOODS FOUND) ###  
  set.seed(123)
  SCP_milo_designT<-SCP_milo_design%>%filter(., Sample_group2=='T')
  SCP_milo_designT$ISUP <- factor(SCP_milo_designT$ISUP, levels = c(1, 2, 5))
  model <- model.matrix(~ 0 + ISUP, data=SCP_milo_designT)
  contrasT<-c("ISUP5 - (ISUP1 + ISUP2)/2")
  mod.contrast<-makeContrasts(contrasts=contrasT, levels=model)

  
  

  
  
  
### DIFFERENTIAL ABUNDANCE ANALYSIS IN CANCER SUBCLUSTERS ###
  
  # As in previous analyses, the procedure followed for the differential abundance analysis, 
  # graphical representation and differential expression of clustered DA cell neighbourhoods 
  #~follows the same procedure as described above. 
  # To avoid redundancy, the comments are reduced to the essential steps
  

  selected_cancer_sce<-as.SingleCellExperiment(selected_cancer)
  selected_cancer_milo<-Milo(selected_cancer_sce)
  
  set.seed(123)
  selected_cancer_milo<-buildGraph(selected_cancer_milo,
                                   k=30,
                                   d=50)
  selected_cancer_milo<-makeNhoods(selected_cancer_milo, 
                                   prop = 0.1,
                                   k =30, 
                                   d=50, refined = TRUE, refinement_scheme="graph")
  selected_cancer_milo<-countCells(selected_cancer_milo,
                                   meta.data = data.frame(colData(selected_cancer_milo)),
                                   samples = 'Sample_Name_Group2')
  plotNhoodSizeHist(selected_cancer_milo)
  head(nhoodCounts(selected_cancer_milo))
  
  # Prepare the design matrix for differencial abundance testing.
  set.seed(123)
  selected_cancer_milo$Sample_group2<-as.factor(selected_cancer_milo$Sample_group2)
  selected_cancer_milo$ISUP<-as.factor(selected_cancer_milo$ISUP)
  nhood_counts <- nhoodCounts(selected_cancer_milo)
  selected_cancer_design<-data.frame(colData(selected_cancer_milo))[,c('Sample_Name_Group2', 
                                                                       'Sample_group2', 
                                                                       'ISUP')]
  selected_cancer_design <- distinct(selected_cancer_design)
  rownames(selected_cancer_design) <- selected_cancer_design$Sample_Name_Group2
  selected_cancer_design  <-selected_cancer_design[colnames(nhood_counts), drop=FALSE]
  

  # Differential abundance testing (edgeR method),
  set.seed(123)
  contrast<-c("Sample_group2H - Sample_group2L")
  da_results_cancer <- testNhoods(selected_cancer_milo, 
                                  design = ~ 0  + Sample_group2, 
                                  design.df = selected_cancer_design, model.contrasts = contrast,
                                  fdr.weighting="graph-overlap", 
                                  norm.method="TMM")
  
  da_results_significative_cancer<-da_results_cancer  %>%
    arrange(SpatialFDR) %>%
    filter(SpatialFDR < 0.1)
  dim(da_results_significative_cancer)
  
  
  ### PLOTTING BEESWARM PLOTS OF DIFF. ABUNDANCE ANALYSIS ###  
  
  # Annotation of cell neighbourhoods based on the subcluster each group belongs #
  da_results_cancer <- annotateNhoods(selected_cancer_milo,
                                      da_results_cancer,
                                      coldata_col = "RNA_snn_res.0.2")
  
  da_results_cancer$RNA_snn_res.0.2 <- ifelse(da_results_cancer$RNA_snn_res.0.2 < 0.7,
                                              "Mixed",
                                              da_results_cancer$RNA_snn_res.0.2)
  
  # Plot Beeswarm plot of the cell neihbourhoods based on subcluster annotation,
  ident_order <- c("1", "2", "3", "4", "Mixed")
  da_results_cancer$RNA_snn_res.0.2 <- factor(da_results_cancer$RNA_snn_res.0.2, levels=ident_order)
  
  b<-plotDAbeeswarm(da_results_cancer, group.by = "RNA_snn_res.0.2", alpha = 0.1) +
    geom_point(size = 3) +
    theme( 
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 16),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14))
  
  
  
  
  ### DIFFERENTIAL EXPRESSION ANALYSIS OF GROUPED DA NEIHBOURHOODS ###  
  
  # Since DA cell neihgbourhoods have been found in cancer subclusters 1 and 3, differential expression analysis 
  # will be aimed at finding marker genes that differentiate DA neihbourhoods from non DA neihbourhoods at the subcluster level.
  
  set.seed(123)
  selected_cancer_milo<- logNormCounts(selected_cancer_milo)
  
  da_results_cancer <- da_results_cancer %>%
    mutate(NhoodGroup = case_when(
      SpatialFDR < 0.1 & RNA_snn_res.0.2 == "1" ~ "Cluster 1 DA",
      SpatialFDR < 0.1 & RNA_snn_res.0.2 == "3" ~ "Cluster 3 DA",
      SpatialFDR > 0.1 & RNA_snn_res.0.2 == "1" ~ "Cluster 1",
      SpatialFDR > 0.1 & RNA_snn_res.0.2 == "3" ~ "Cluster 3",
      TRUE ~ "0"
    ))
  
  
  keep.rows <- rowSums(logcounts(selected_cancer_milo)) != 0
  selected_cancer_milo <- selected_cancer_milo[keep.rows, ]
  

  dec <- modelGeneVar(selected_cancer_milo)
  hvgs <- getTopHVGs(dec)
  
  
  test_cancer<-findNhoodGroupMarkers(selected_cancer_milo, da_results_cancer, subset.row = hvgs,
                                     aggregate.samples=TRUE, sample_col = "Sample_Name_Group2")
  
  # Differential expression testing between subcluster 1 DA and non DA grouped cell neighbourhoods.
  test_cancerct1<-findNhoodGroupMarkers(
    selected_cancer_milo, 
    da_results_cancer, 
    subset.row = hvgs,
    aggregate.samples = TRUE,
    subset.nhoods = da_results_cancer$NhoodGroup %in% c("Cluster 1", "Cluster 1 DA"),
    sample_col = "Sample_Name_Group2")
  
  # Differential expression testing between subcluster 3 DA and non DA grouped cell neighbourhoods.
  test_cancerct3<-findNhoodGroupMarkers(
    selected_cancer_milo, 
    da_results_cancer, 
    subset.row = hvgs,
    aggregate.samples = TRUE,
    subset.nhoods = da_results_cancer$NhoodGroup %in% c("Cluster 3", "Cluster 3 DA"),
    sample_col = "Sample_Name_Group2")
  
  
  # Get the resulting marker genes.
  markers_cancerct1<-test_cancerct1%>%filter(., `adj.P.Val_Cluster 1 DA`<0.05 & abs(`logFC_Cluster 1 DA`) > 1.5)
  markers_cancerct3<-test_cancerct3%>%filter(., `adj.P.Val_Cluster 3 DA`<0.05 & abs(`logFC_Cluster 3 DA`) > 1.5)
  
  
  # Plot the gene expression heatmaps (following guideliness in https://www.bioconductor.org/packages/release/bioc/vignettes/miloR/inst/doc/milo_gastrulation.html)
  htb<-plotNhoodExpressionGroups(selected_cancer_milo, da_results_cancer,
                                 scale=TRUE, 
                                 features = markers_cancerct1$GeneID, 
                                 grid.space = 'fixed',
                                 subset.nhoods = da_results_cancer$NhoodGroup %in% c("Cluster 1", "Cluster 1 DA")) # subset of nhoods to show in plot 
  
  htc<-plotNhoodExpressionGroups(selected_cancer_milo, da_results_cancer,
                                 scale=TRUE,
                                 features = markers_cancerct3$GeneID,
                                 grid.space = 'fixed',
                                 subset.nhoods = da_results_cancer$NhoodGroup %in% c("Cluster 3", "Cluster 3 DA") )
  
  
  
  # Differential abundance analysis separating by ISUP did not yield significant results.
  set.seed(123)
  selected_cancer_designT<-selected_cancer_design%>%filter(., Sample_group2=='T')
  selected_cancer_designT$ISUP <- factor(selected_cancer_designT$ISUP, levels = c(1, 2, 5))
  model <- model.matrix(~ 0 + ISUP, data=selected_cancer_designT)
  contrasT<-c("ISUP5 - (ISUP1 + ISUP2)/2")
  mod.contrast<-makeContrasts(contrasts=contrasT, levels=model)
  
  
  
### DIFFERENTIAL ABUNDANCE ANALYSIS IN CAF SUBCLUSTERS (NO SIGNIFICATIVE RESULTS WERE FOUND) ###  
  
  selected_caf_sce<-as.SingleCellExperiment(selected_caf)
  selected_caf_milo<-Milo(selected_caf_sce)
  
  set.seed(123)
  selected_caf_milo<-buildGraph(selected_caf_milo,
                                k=20,
                                d=50)
  selected_caf_milo<-makeNhoods(selected_caf_milo, 
                                prop = 0.1,
                                k =20, 
                                d=50, refined = TRUE, refinement_scheme="graph")
  selected_caf_milo<-countCells(selected_caf_milo,
                                meta.data = data.frame(colData(selected_caf_milo)),
                                samples = 'Sample_Name_Group2')
  plotNhoodSizeHist(selected_caf_milo)
  head(nhoodCounts(selected_caf_milo))
  
  
  set.seed(123)
  selected_caf_milo$Sample_group2<-as.factor(selected_caf_milo$Sample_group2)
  selected_caf_milo$ISUP<-as.factor(selected_caf_milo$ISUP)
  nhood_counts <- nhoodCounts(selected_caf_milo)
  
  selected_caf_design<-data.frame(colData(selected_caf_milo))[,c('Sample_Name_Group2', 
                                                                 'Sample_group2', 
                                                                 'ISUP')]
  selected_caf_design <- distinct(selected_caf_design)
  rownames(selected_caf_design) <- selected_caf_design$Sample_Name_Group2
  selected_caf_design  <-selected_caf_design[colnames(nhood_counts), , drop=FALSE]
  
  contrast<-c("Sample_group2T - Sample_group2N")
  da_results_caf <- testNhoods(selected_caf_milo, 
                               design = ~ 0 + Sample_group2, 
                               design.df = selected_caf_design , model.contrasts = contrast,
                               fdr.weighting="graph-overlap", 
                               norm.method="TMM")

