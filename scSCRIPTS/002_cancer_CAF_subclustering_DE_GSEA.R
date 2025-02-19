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


# Single cell analysis + differential expression analysis.
library(Seurat)
library(agrmt)
library(SeuratObject)
library(scran)
source('CalculateCellularDetectionRate.R')
library(MAST)
library(agrmt) # Calculate Concentration and Dispersion in Ordered Rating Scales

# ssGSEA over Seurat's object.
library(escape)
library(RColorBrewer)



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



# ORIGINAL STEPS FOR THIS SCRIPT


### CANCER RECLUSTERING ###
selected_cancer <- subset(SCP_data, 
                          cells = WhichCells(SCP_data, expression = CAF_can_norm == "Cancer")) # Extract from original Seurat object those cells annotated as Cancer.

  options(future.globals.maxSize = 1e9)
  selected_cancer <- selected_cancer %>% 
    NormalizeData(., scale.factor = 1e6, assay = "RNA", verbose = FALSE) %>% # scale.factor = 1e6 means we are using CPMs.
    FindVariableFeatures(., nfeatures = 3000, assay = "RNA", verbose = FALSE) %>% # Use top 3000 hvgs to perform scaling.
    ScaleData(., assay = "RNA", vars.to.regress = "percent.mt", verbose = FALSE) %>% # regress out % mito.
    RunPCA(., assay = "RNA", verbose = FALSE)
  
  ElbowPlot(selected_cancer, ndims = 50) # We'll select the firts 20 PCs based on the elbow plot results.


  resolutions = c(0.2, 0.4, 1) # Set the resolution parameters. The most interesting will be the lowest one as it will group larger groups of cells.
  selected_cancer <- selected_cancer %>% # Final Seurat object with filtered and reclustered cancer cells.
    FindNeighbors(., dims = 1:20, reduction = "harmony", verbose = FALSE) %>%
    FindClusters(., resolution =  resolutions, verbose = FALSE) %>%
    RunUMAP(., dims = 1:20, reduction = "harmony", verbose = FALSE)
  

### CANCER RECLUSTERED DATASET PLOTS ###
  
  DimPlot(selected_cancer, # Selected cancer Seurat object.
          reduction = "umap", 
          group.by = "RNA_snn_res.0.2", # We'll use the resolution of 0.2
          label.size = 6,
          pt.size = 1.5,
          label.color = "black") +
    xlab("UMAP 1") +  
    ylab("UMAP 2") +  
    ggtitle(NULL) +  
    theme(axis.title = element_text(size = 15), 
          legend.text = element_text(size = 12), 
          axis.text = element_text(size = 12),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3))
  
  
  # Barplot of cell percentages by cluster, three-group layer and class of sample.
  cluster_counts_cancerclust <- selected_cancer@meta.data %>%
    group_by(Sample_Name_Group2, RNA_snn_res.0.2, Sample_group2) %>%
    dplyr::summarise(count = n()) %>%
    ungroup() # Count number of cells per cancer subcluster (0.2) and per sample while keeping the sample class.
  
  cluster_percentages_cancerclust <- cluster_counts_cancerclust %>%
    group_by(Sample_Name_Group2) %>%
    mutate(total_cells = sum(count),  
           percentage = (count / total_cells) * 100) %>%
    ungroup() # Get the percentage of each cancer subcluster per sample while keeping the sample class.
  
  # Building the barplot (Figure 2)
  cluster_percentages_cancerclust$RNA_snn_res.0.2 <- factor(cluster_percentages_cancerclust$RNA_snn_res.0.2, 
                                                            levels = as.character(0:4))
  ggplot(cluster_percentages_cancerclust, aes(x = Sample_Name_Group2, y = percentage, fill = RNA_snn_res.0.2)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +  
    labs(y = "Percentage of cells", x = 'Sample') +
    facet_grid(RNA_snn_res.0.2 ~ Sample_group2, scales = "free_x") + # Barplot splitting by cancer subcluster and sample class
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, angle = 65, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 18),
          legend.title = element_blank(),
          strip.text = element_text(size = 18, face = "bold"),
          strip.text.y = element_blank())  
  
  
  # Summary table of descriptive statistics.
  
  # Get the count of cells by ISUP, cancer subcluster and sample class.
  cluster_counts_cancerclust2 <- selected_cancer@meta.data %>% 
    select(Sample_Name_Group2, RNA_snn_res.0.2, Sample_group2, ISUP)%>%
    group_by(Sample_Name_Group2, RNA_snn_res.0.2, Sample_group2, ISUP) %>%
    summarise(count = n())%>%
    ungroup()
  
  # Transform the cell count dataframe to percentage.
  cluster_percentages_cancerclust2 <- cluster_counts_cancerclust2 %>%
    group_by(Sample_Name_Group2) %>%
    mutate(total_cells = sum(count),  
           percentage = (count / total_cells) * 100) %>%
    ungroup()
  
  # Get the basic statistics.
  cluster_percentages_cancerclust2_SUMM<-cluster_percentages_cancerclust2 %>%
    group_by(RNA_snn_res.0.2, ISUP, Sample_group2)%>%
    summarise(Mean = mean(percentage, na.rm = TRUE), # Mean.
              SD = sd(percentage, na.rm = TRUE), # SD.
              Min = min(percentage, na.rm = TRUE),# Min.
              Max = max(percentage, na.rm = TRUE)) # Max.
  
  # Table transformation to publication-ready report
  cluster_percentages_cancerclust2_SUMM_gt<-cluster_percentages_cancerclust2_SUMM %>%
    mutate(
      Stats = sprintf("%.2f (%.2f) (%.1f, %.1f)", Mean, SD, Min, Max) 
    ) %>%
    select(-Mean, -SD, -Min, -Max) %>% 
    pivot_wider(names_from = ISUP, values_from = Stats)
  
  
### DIFFERENTIAL EXPRESSION ANALYSIS ###
  
  Idents(selected_cancer) = "RNA_snn_res.0.2" # Differential expression analysis based on cancer cell subclusters.
  selected_caf$Sample_group2<-as.factor(selected_caf$Sample_group2) # Factrize covariates.
  selected_caf$ISUP<-as.factor(selected_caf$ISUP)
  

  
  cancer_clust_markers<-FindAllMarkers(object = selected_cancer,
                                       test.use = "MAST", # Use MAST as DE method.
                                       latent.vars = c("orig.ident", "CDR", "Sample_group2", "ISUP"), # Correct by CDR and batch effect.
                                       assay = "RNA",
                                       slot = "data", # Use log-transformed CPMs.
                                       only.pos = FALSE, # Cell markers can only be positive markers.
                                       min.pct = 0.5, # We consider that a marker gene should be expressed in at least 50% of the cells.
                                       logfc.threshold = log2(1.5), # Filter by FC threshold.
                                       verbose = FALSE)
  
  cancer_clustall_top<-cancer_clust_markers%>%group_by(., cluster)%>%top_n(n = 50, wt = avg_log2FC) # Select top 50 overexpressed genes in each cell type.
  cancer_clustall_bottom<-cancer_clust_markers%>%slice_min(n = 50, order_by = avg_log2FC) # Select top 50 underexpressed genes in each cell type.
  subcancer_txt<-cancer_clustall_top%>%filter(., gene %in% rownames(selected_cancer@assays$RNA$scale.data)) # Select those top 50 overexpressed genes also present in the scaled data (genes hvgs, relevant).
  
  
  # Heatmap of DE genes (Fig 2)
  DoHeatmap(selected_cancer, slot = 'scale.data', assay = 'RNA', # Seurat's gene expression heatmap function.
            features = cancer_clustall_top$gene,  # Plotting top overexpressed genes per cell type.
            group.bar = TRUE,  
            label = FALSE, group.colors = cluster_colors) +
    scale_fill_gradient2(low = "#4B79A5", mid = "white", high = "#D23C44",  # Custom scale. Scaled gene exp. values = 0 are white, negative blue and positive red.
                         midpoint = 0) +
    theme(legend.text = element_text(size = 12), 
          legend.title = element_text(size = 16),
          axis.text.y = element_text(size = 6.5))


### FUNCTIONAL ANNOTATION BY ssGSEA (escape) ###

  GS.hallmark <- getGeneSets(library = "H") # Set MiSigDB gene set hallmark H (https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=H).
  Idents(SCP_data) = "RNA_snn_res.0.2" # Set cell types to compare (Cancer subclusters).
  
  selected_cancer <- runEscape(selected_cancer , 
                               method = "ssGSEA",
                               gene.sets = GS.hallmark, 
                               groups = 1000, # The number of cells to separate the enrichment calculation (default https://bioconductor.org/packages/release/bioc/vignettes/escape/inst/doc/vignette.html).
                               min.size = 5, # Minimum number of gene necessary to perform the enrichment calculation (default).
                               new.assay.name = "escape.ssGSEA") # Slot in assay section in Seurat's object to store de cell x pathway menrichmen score matrix.
  
  # performNormalization() will normalize the enrichment values by calculating the number of genes expressed in each gene set and cell.
  selected_cancer <- performNormalization(selected_cancer, 
                                          assay = "escape.ssGSEA", 
                                          gene.sets = GS.hallmark, 
                                          scale.factor = selected_cancer$nFeature_RNA)
  
  # Heatmap of mollecular pathways by normalized enrichment score and grouped by cell type. (Fig 1)
  # Same procedure as guideliness, rest of arguments set as default.
  heatmapEnrichment(selected_cancer,
                    group.by = "RNA_snn_res.0.2",
                    gene.set.use = 'all', # Selected gene sets to visualize. If "all", the heatmap will be generated across all gene sets.
                    assay = "escape.ssGSEA_normalized", # Use the normalized enrichment score matrix.
                    cluster.rows = TRUE) +
    scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) + # Set different RColorbrewer scale.
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),
          axis.text.y = element_text(size=8))
  
 


  # UMAPs of key pathways with changes of the normalized enrichment score between subclusters. also related with cancer pathways
  DefaultAssay(selected_cancer) <- "escape.ssGSEA_normalized" # We'll use the normalized data
  colorblind_vector <- hcl.colors(n=7, palette = "RdYlBu", fixup = TRUE)
  FeaturePlot(selected_cancer, pt.size = 1.5,
                   features="HALLMARK-REACTIVE-OXYGEN-SPECIES-PATHWAY") + 
    scale_color_gradientn(colors = rev(colorblind_vector)) + 
    xlab("UMAP 1") +  
    ylab("UMAP 2") +  
    ggtitle("HALLMARK-REACTIVE-OXYGEN-SPECIES") +  
    theme(axis.title = element_text(size = 18), 
          legend.text = element_text(size = 12), 
          axis.text = element_text(size = 12),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3))+
    theme(strip.text = element_text(size = 18, face = "bold"))
  
  
  FeaturePlot(selected_cancer, pt.size = 1.5,
                   features="HALLMARK-MYC-TARGETS-V1") + 
    scale_color_gradientn(colors = rev(colorblind_vector)) + 
    xlab("UMAP 1") +  
    ylab("UMAP 2") +  
    ggtitle("HALLMARK-MYC-TARGETS-V1") +  
    theme(axis.title = element_text(size = 18), 
          legend.text = element_text(size = 12), 
          axis.text = element_text(size = 12),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3))+
    theme(strip.text = element_text(size = 18, face = "bold"))
  
  
  FeaturePlot(selected_cancer, pt.size = 1.5,
                   features="HALLMARK-MTORC1-SIGNALING") + 
    scale_color_gradientn(colors = rev(colorblind_vector)) + 
    xlab("UMAP 1") +  
    ylab("UMAP 2") +  
    ggtitle("HALLMARK-MTORC1-SIGNALING") +  
    theme(axis.title = element_text(size = 18), 
          legend.text = element_text(size = 12), 
          axis.text = element_text(size = 12),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3))+
    theme(strip.text = element_text(size = 18, face = "bold"))
  
  
  
  
  
### CAF RECLUSTERING (SEURAT SUBCLUSTER ANNOTATION)###
  
  # IMPORTANT: We performed the same set of analyses as indicated for Cancer cells, but analogously for CAF cells. This is why the code documentation is omitted to avoid redundancy.  
  #selected_caf <- subset(SCP_data, cells = WhichCells(SCP_data, expression = CAF_can_norm == "CAF"))
  
  selected_caf <- selected_caf %>% 
    NormalizeData(., scale.factor = 1e6, assay = "RNA", verbose = FALSE) %>% # scale.factor = 1e6 means we are using CPMs.
    FindVariableFeatures(., nfeatures = 3000, assay = "RNA", verbose = FALSE) %>%
    ScaleData(., assay = "RNA", vars.to.regress = "percent.mt", verbose = FALSE) %>% # regress out % mito.
    RunPCA(., assay = "RNA", verbose = FALSE)
  
  ElbowPlot(selected_caf, ndims = 50)
  
  resolutions = c(0.2, 0.1)
  
  selected_caf <- selected_caf %>% 
    FindNeighbors(., dims = 1:20, reduction = "harmony", verbose = FALSE) %>%
    FindClusters(., resolution =  resolutions, verbose = FALSE) %>% 
    RunUMAP(., dims = 1:20, reduction = "harmony", verbose = FALSE)
  
  
  DimPlot(selected_caf, # Figure 3
          reduction = "umap", 
          group.by = "RNA_snn_res.0.2", 
          label.size = 6, pt.size = 1.5,
          label.color = "black") +
    xlab("UMAP 1") +  
    ylab("UMAP 2") +  
    ggtitle(NULL) +  
    theme(axis.title = element_text(size = 15), 
          legend.text = element_text(size = 12), 
          axis.text = element_text(size = 12),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3))
  
### CAF RECLUSTERED DATASET PLOTS ###
  
  cluster_counts_cafclust <- selected_caf@meta.data %>%
    group_by(Sample_Name_Group2, RNA_snn_res.0.2, Sample_group2) %>%
    summarise(count = n()) %>%
    ungroup()
  
  cluster_percentages_cafclust <- cluster_counts_cafclust %>%
    group_by(Sample_Name_Group2) %>%
    mutate(total_cells = sum(count),  
           percentage = (count / total_cells) * 100) %>%
    ungroup()
  
  cluster_percentages_cafclust$RNA_snn_res.0.2 <- factor(cluster_percentages_cafclust$RNA_snn_res.0.2, 
                                                         levels = as.character(0:3))
  
  # Building the barplot (supp Fig 2)
  ggplot(cluster_percentages_cafclust, aes(x = Sample_Name_Group2, y = percentage, fill = RNA_snn_res.0.2)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +  
    labs(y = "Percentage of cells", x = 'Sample') +
    facet_grid(RNA_snn_res.0.2 ~ Sample_group2, scales = "free_x") + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, angle = 65, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 18),
          legend.title = element_blank(),
          strip.text = element_text(size = 18, face = "bold"),
          strip.text.y = element_blank())
  
  # Summary table of descriptive statistics.
  cluster_counts_cafclust2 <- selected_caf@meta.data %>% 
    select(Sample_Name_Group2, RNA_snn_res.0.2, Sample_group2, ISUP)%>%
    group_by(Sample_Name_Group2, RNA_snn_res.0.2, Sample_group2, ISUP) %>%
    summarise(count = n())%>%
    ungroup()
  
  cluster_percentages_cafclust2 <- cluster_counts_cafclust2 %>%
    group_by(Sample_Name_Group2) %>%
    mutate(total_cells = sum(count),  
           percentage = (count / total_cells) * 100) %>%
    ungroup()
  
  cluster_percentages_cafclust2_SUMM<-cluster_percentages_cafclust2 %>%
    group_by(RNA_snn_res.0.2, ISUP, Sample_group2)%>%
    summarise(Mean = mean(percentage, na.rm = TRUE),
              SD = sd(percentage, na.rm = TRUE),
              Min = min(percentage, na.rm = TRUE),
              Max = max(percentage, na.rm = TRUE))
  
  cluster_percentages_cafclust2_SUMM_gt<-cluster_percentages_cafclust2_SUMM %>%
    mutate(
      Stats = sprintf("%.2f (%.2f) (%.1f, %.1f)", Mean, SD, Min, Max) 
    ) %>%
    select(-Mean, -SD, -Min, -Max) %>% 
    pivot_wider(names_from = ISUP, values_from = Stats)
  
  
### DIFFERENTIAL EXPRESSION ANALYSIS ###    
  
  Idents(selected_caf) = "RNA_snn_res.0.2"
  selected_caf$Sample_group2<-as.factor(selected_caf$Sample_group2)
  selected_caf$ISUP<-as.factor(selected_caf$ISUP)
  
  options(future.globals.maxSize = 1024 * 1024 * 1024)
  
  caf_clust_markers<-FindAllMarkers(object = selected_caf,
                                    test.use = "MAST", 
                                    latent.vars = c("orig.ident", "CDR", "Sample_group2", "ISUP"), 
                                    assay = "RNA",
                                    slot = "data", 
                                    only.pos = FALSE, 
                                    min.pct = 0.5, 
                                    logfc.threshold = log2(1.5), 
                                    verbose = FALSE)
  
  caf_clustall_top<-caf_clust_markers%>%group_by(., cluster)%>%top_n(n = 50, wt = avg_log2FC)
  caf_clustall_bottom<-caf_clust_markers%>%slice_min(n = 50, order_by = avg_log2FC)
  subcaf_txt<-caf_clustall_top%>%filter(., gene %in% rownames(selected_caf@assays$RNA$scale.data))
  
  DoHeatmap(selected_caf, slot = 'scale.data', assay = 'RNA', # Figure 3
            features = caf_clustall_top$gene,  
            group.bar = TRUE,  
            label = FALSE, group.colors = cluster_colors) +
    scale_fill_gradient2(low = "#4B79A5", mid = "white", high = "#D23C44", 
                         midpoint = 0) +
    theme(legend.text = element_text(size = 12), 
          legend.title = element_text(size = 16),
          axis.text.y = element_text(size = 6.5))
  
  
### CAF RECLUSTERING (scATOMIC CAF SUBPOPULATION ANNOTATION)###
  
  # We take a second approach to characterising CAF cells by taking 
  # advantage of the fact that scATOMIc classifies them into different subpopulations.
  
  ### CAF RECLUSTERED DATASET PLOTS ###
  
  DimPlot(selected_caf, # Fig 3
          reduction = "umap", 
          group.by = "scATOMIC_pred", 
          label.size = 6, pt.size = 1.5,
          label.color = "black") +
    xlab("UMAP 1") +  
    ylab("UMAP 2") +  
    ggtitle(NULL) +  
    theme(axis.title = element_text(size = 15), 
          legend.text = element_text(size = 12), 
          axis.text = element_text(size = 12),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3))
  
  cluster_counts_cafatomic <- selected_caf@meta.data %>%
    group_by(Sample_Name_Group2, scATOMIC_pred, Sample_group2) %>%
    summarise(count = n()) %>%
    ungroup()
  
  cluster_percentages_cafatomic <- cluster_counts_cafatomic %>%
    group_by(Sample_Name_Group2) %>%
    mutate(total_cells = sum(count),  # Número total de células por muestra y grupo
           percentage = (count / total_cells) * 100) %>%
    ungroup()
  
  # Barplot for supp fig 2
  ggplot(cluster_percentages_cafatomic, aes(x = Sample_Name_Group2, y = percentage, fill = scATOMIC_pred)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +  
    labs(y = "Percentage of cells", x = 'Sample') +
    facet_grid(scATOMIC_pred ~ Sample_group2, scales = "free_x") + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, angle = 65, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 18),
          legend.title = element_blank(),
          strip.text = element_text(size = 18, face = "bold"),
          strip.text.y = element_blank())
  
  # Summary table by ISUP and sample class
  
  cluster_counts_cafatomic2 <- selected_caf@meta.data %>% 
    select(Sample_Name_Group2, scATOMIC_pred, Sample_group2, ISUP)%>%
    group_by(Sample_Name_Group2, scATOMIC_pred, Sample_group2, ISUP) %>%
    summarise(count = n())%>%
    ungroup()
  
  cluster_percentages_cafatomic2 <- cluster_counts_cafatomic2 %>%
    group_by(Sample_Name_Group2) %>%
    mutate(total_cells = sum(count),  
           percentage = (count / total_cells) * 100) %>%
    ungroup()
  
  cluster_percentages_cafatomic2_SUMM<-cluster_percentages_cafatomic2 %>%
    group_by(scATOMIC_pred, ISUP, Sample_group2)%>%
    summarise(Mean = mean(percentage, na.rm = TRUE),
              SD = sd(percentage, na.rm = TRUE),
              Min = min(percentage, na.rm = TRUE),
              Max = max(percentage, na.rm = TRUE))
  
  cluster_percentages_cafatomic2_SUMM_gt<-cluster_percentages_cafatomic2_SUMM %>%
    mutate(
      Stats = sprintf("%.2f (%.2f) (%.1f, %.1f)", Mean, SD, Min, Max) 
    ) %>%
    select(-Mean, -SD, -Min, -Max) %>% 
    pivot_wider(names_from = ISUP, values_from = Stats)

### FUNCTIONAL ANNOTATION BY ssGSEA (escape) ###
  
  GS.hallmark <- getGeneSets(library = "H")
  
  
  selected_caf <- runEscape(selected_caf , 
                            method = "ssGSEA",
                            gene.sets = GS.hallmark, 
                            groups = 1000, 
                            min.size = 5,
                            new.assay.name = "escape.ssGSEA")
  
  
  selected_caf <- performNormalization(selected_caf, 
                                       assay = "escape.ssGSEA", 
                                       gene.sets = GS.hallmark, 
                                       scale.factor = selected_caf$nFeature_RNA)
  
  # Heatmap of normalized enrichmen scores of each pathway grouped by CAF subcluster (Fig 3)
  heatmapEnrichment(selected_caf,
                    group.by = "RNA_snn_res.0.2",
                    gene.set.use = 'all',
                    assay = "escape.ssGSEA_normalized", 
                    cluster.rows = TRUE) +
    scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),
          axis.text.y = element_text(size=8))
  
  # Heatmap of normalized enrichmen scores of each pathway grouped by CAF subpopulation obtained by scATOMIC (supp Fig 2)
  heatmapEnrichment(selected_caf,
                    group.by = "scATOMIC_pred",
                    gene.set.use = 'all',
                    assay = "escape.ssGSEA_normalized", 
                    cluster.rows = TRUE) +
    scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),
          axis.text.y = element_text(size=8))


  # UMAP plots of the key mollecular pathways related to cancer found in CAF subpopulations
  
  DefaultAssay(selected_caf) <- "escape.ssGSEA_normalized"
  colorblind_vector <- hcl.colors(n=7, palette = "RdYlBu", fixup = TRUE)
  FeaturePlot(selected_caf, pt.size = 1.5,
                   features="HALLMARK-TGF-BETA-SIGNALING") + 
    scale_color_gradientn(colors = rev(colorblind_vector)) + 
    xlab("UMAP 1") +  
    ylab("UMAP 2") +  
    ggtitle("HALLMARK-TGF-BETA-SIGNALING") +  
    theme(axis.title = element_text(size = 18), 
          legend.text = element_text(size = 12), 
          axis.text = element_text(size = 12),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3))+
    theme(strip.text = element_text(size = 18, face = "bold"))
  
  
  FeaturePlot(selected_caf, pt.size = 1.5,
                   features="HALLMARK-ANGIOGENESIS") + 
    scale_color_gradientn(colors = rev(colorblind_vector)) + 
    xlab("UMAP 1") +  
    ylab("UMAP 2") +  
    ggtitle("HALLMARK-ANGIOGENESIS") +  
    theme(axis.title = element_text(size = 18), 
          legend.text = element_text(size = 12), 
          axis.text = element_text(size = 12),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3))+
    theme(strip.text = element_text(size = 18, face = "bold"))
  
  
  FeaturePlot(selected_caf, pt.size = 1.5,
                   features="HALLMARK-EPITHELIAL-MESENCHYMAL-TRANSITION") + 
    scale_color_gradientn(colors = rev(colorblind_vector)) + 
    xlab("UMAP 1") +  
    ylab("UMAP 2") +  
    ggtitle("HALLMARK-EPITHELIAL-MESENCHYMAL-TRANSITION") +  
    theme(axis.title = element_text(size = 18), 
          legend.text = element_text(size = 12), 
          axis.text = element_text(size = 12),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3))+
    theme(strip.text = element_text(size = 18, face = "bold"))
