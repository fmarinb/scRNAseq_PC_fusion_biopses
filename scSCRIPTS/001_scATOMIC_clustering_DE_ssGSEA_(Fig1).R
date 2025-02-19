### PACKAGE LOADING ###

# Cancer-oriented annotation via scATOMIC + dependencies
library(scATOMIC)
library(cutoff.scATOMIC)
library(copykat)
library(caret)
library(randomForest)
library(agrmt)
library(Rmagic)
library(Matrix)

# Datasets handling and additional plots.
library(dplyr)
library(data.table)
library(parallel)
library(reticulate)
library(ggplot2)

# Single cell analysis + differential expression analysis.
library(Seurat)
library(SeuratObject)
library(scran)
source('CalculateCellularDetectionRate.R')
library(MAST)

# ssGSEA over Seurat's object.
library(escape)
library(RColorBrewer)


#### LOADING THE PRE-PROCESSED SEURAT FILE ###
    options(parallelly.fork.enable = TRUE) # Set up Seurat pararell computing.
    plan("multicore", workers = parallel::detectCores())
    SCP_data <- readRDS(file = "SC_Prostate_processed.rds") # Load seurat object with 0.4 resolution.
    SCP_data <- JoinLayers(SCP_data, assay = "RNA") # Join the per-run splitted layers.
    original_ident<-SCP_data$orig.ident
    rownm<-rownames(SCP_data@meta.data)
    SCP_data <- CalculateCDR(SCP_data) # Calculate Cell Detection Rate (CDR) with our own function (code available in a separate file).



#### CANCER-ORIENTED ANNOTATION WITH scATOMIC ###
    # Combine the datasets by normalized-unscaled counts an per run.
    count_run1 <- as.matrix(GetAssayData(SCP_data, layer = "counts.run1", assay = "RNA"))
    count_run2 <- as.matrix(GetAssayData(SCP_data, layer = "counts.run2", assay = "RNA"))
    count_run3 <- as.matrix(GetAssayData(SCP_data, layer = "counts.run3", assay = "RNA"))
    count_run4 <- as.matrix(GetAssayData(SCP_data, layer = "counts.run4", assay = "RNA"))
    combined_normalized_data <- cbind(count_run1, count_run2, count_run3,  count_run4)
    
    
    # Data pre-rpocessing following scATOMIC developer's guideliness,
    pct_mt <- colSums(combined_normalized_data[grep("^MT-", row.names(combined_normalized_data)),])/colSums(
      combined_normalized_data) * 100 # get the % of mitochondrial RNA
    nFeatureRNA <- colSums(combined_normalized_data> 0) # Filtering out genes with zero counts
    combined_normalized_data <- combined_normalized_data[, names(which(pct_mt < 25))] # Filtering out cells with a % of mitochondrial RNA >25%.
    combined_normalized_data <- combined_normalized_data[, intersect(names(which(nFeatureRNA > 500)), 
                                                                     colnames(combined_normalized_data))] # Filtering out cells with fewer than 500 unique features expressed.
    
    # Run scATOMIC and get the annotation object.
    cell_predictions <- run_scATOMIC(combined_normalized_data,
                                     mc.cores = 4,
                                     imputation = T) 
    options(future.globals.maxSize = 9 * 1024^3)  # Matrix handling.
    
    
    # Create the summary of the results obtained in cell_predictions(). 
    # Default parameters specified in https://github.com/abelson-lab/scATOMIC
    results_cap <- create_summary_matrix(
      prediction_list = cell_predictions, 
      use_CNVs = F,
      modify_results = T,
      mc.cores = 4, 
      raw_counts = combined_normalized_data, 
      min_prop = 0.5,
      known_cancer_type = "PC cell", # We know the type of cancer.
      confidence_cutoff = T,
      pan_cancer = F 
    ) 
    
    
    # Data visualization of the scATOMIC results via tree plot,
    tree_results_interactive <- scATOMICTree(predictions_list = cell_predictions,
                                             summary_matrix = results_cap, 
                                             interactive_mode = T,
                                             collapsed = T, save_results = T,
                                             height = 700, width = 700)
    
    SCP_data <- AddMetaData(SCP_data , results_cap) # Add scATOMIC matrix in Seurat's metadata slot 
    SCP_data$orig.ident<-original_ident # Maintain original ident names instead of those renamed by scATOMIC
    


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
    
    
### CELL ANNOTATION PLOTS ###
    
    # UMAP plot of the annotated cell types by Single R (supp Fig1, done during pre-processing steps).
    DimPlot(SCP_data, 
            reduction = "umap", 
            group.by = "cell.labels", 
            label = TRUE, 
            repel = TRUE, 
            label.size = 5,  
            label.color = "black") +  
      xlab("UMAP 1") +  
      ylab("UMAP 2") +  
      ggtitle(NULL) +  
      theme(
        axis.title = element_text(size = 15), 
            legend.text = element_text(size = 12), 
            axis.text = element_text(size = 12)) + 
      guides(colour = guide_legend(override.aes = list(size = 4))) 
    
    # UMAP plot of the cells classified based on Seurat's clusters (res 0.4, supp Fig1).
    DimPlot(SCP_data,
            reduction = "umap", 
            group.by = "RNA_snn_res.0.4", 
            label = TRUE, 
            repel = TRUE, 
            label.size = 6,  
            label.color = "black") +  
      xlab("UMAP 1") +  
      ylab("UMAP 2") +  
      ggtitle(NULL) +  
      theme(
        axis.title = element_text(size = 15), 
            legend.text = element_text(size = 12), 
            axis.text = element_text(size = 12)) + 
      guides(colour = guide_legend(override.aes = list(size = 4))) 
    
    # UMAP plot using the detailed annotation os scATOMIC (supp Fig 1).
    DimPlot(SCP_data, 
            reduction = "umap", 
            group.by = "scATOMIC_pred", 
            label.size = 6,  
            label.color = "black") +  
      xlab("UMAP 1") +  
      ylab("UMAP 2") +  
      ggtitle(NULL) +  
      theme(axis.title = element_text(size = 15), 
            legend.text = element_text(size = 12), 
            axis.text = element_text(size = 12)) + 
      guides(colour = guide_legend(override.aes = list(size = 4))) 
    

    # UMAP plot using the three type annotation layer (Fig 1).
    DimPlot(SCP_data, 
            reduction = "umap", 
            group.by = "CAF_can_norm", 
            label.size = 6,  
            label.color = "black",
            cols = c("#CE7D3B",  "#DD64BE", "#BBC5CC")) +    
      xlab("UMAP 1") +  
      ylab("UMAP 2") +  
      ggtitle(NULL) +  
      theme(
        text = element_text(family = "sans"),
            axis.title = element_text(size = 16), 
            legend.text = element_text(size = 12), 
            axis.text = element_text(size = 12),
            panel.border = element_rect(colour = "black", fill=NA, size=0.3))+
      theme(strip.text = element_text(size = 18, face = "bold"))
    
    # Barplot of cell percentages by cluster, three-group layer and class of sample.
    cluster_counts_cancer <- SCP_data@meta.data %>%
      group_by(Sample_Name_Group2, RNA_snn_res.0.4, CAF_can_norm, Sample_group2) %>%
      summarise(count = n())%>%
      ungroup() # Count numer of CAF, Cancer, No Cancer and NA cells per Seurat's cluster and sample while keeping the sample class.
    
    cluster_percentages_cancer <- cluster_counts_cancer %>%
      group_by(Sample_Name_Group2) %>%
      mutate(total_cells = sum(count),  
             percentage = (count / total_cells) * 100) %>%
      ungroup() # Get the percentage of each CAF, Cancer and No Cancer cell type per Seurat's cluster and sample while keeping the sample class.


    # Recalculate the cell percentage of each cell type without the SEURAT's cluster information.
    cluster_counts_cancer2 <- SCP_data@meta.data %>% 
      select(Sample_Name_Group2, CAF_can_norm, Sample_group2, ISUP)%>%
      group_by(Sample_Name_Group2, CAF_can_norm, Sample_group2, ISUP) %>%
      summarise(count = n())%>%
      ungroup()
    
    cluster_percentages_cancer2 <- cluster_counts_cancer2 %>%
      group_by(Sample_Name_Group2) %>%
      mutate(total_cells = sum(count),
             percentage = (count / total_cells) * 100) %>%
      ungroup()
    
    # Building the barplot (Figure 1)
    cluster_percentages_cancer$RNA_snn_res.0.4 <- factor(cluster_percentages_cancer$RNA_snn_res.0.4, levels = as.character(0:14))
    ggplot(cluster_percentages_cancer, aes(x = Sample_Name_Group2, y = percentage, fill = RNA_snn_res.0.4)) + # Split by Seurat cluster.
      geom_bar(stat = "identity", position = "stack", color = "black") +  
      labs(y = "Percentage of cells", x='Sample') +
      facet_grid(CAF_can_norm ~ Sample_group2, scales = "free_x") + # Barplot splitting by cell type and sample class.
      theme_bw() +
      theme( 
        text = element_text(family = "sans"),
        axis.text.x = element_text(size = 8, angle = 65, vjust = 1, hjust=1),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=18))+
      theme(legend.title = element_blank())+
      theme(strip.text = element_text(size = 18, face = "bold"))
    
    # Summary table of descriptive statistics
    cluster_percentages_cancer2_SUMM<-cluster_percentages_cancer2 %>%
      group_by(CAF_can_norm, ISUP, Sample_group2)%>% # group by CAF, Cancer and No Cancer cell type, ISUP and sample class.
      summarise(Mean = mean(percentage, na.rm = TRUE), # Mean.
                SD = sd(percentage, na.rm = TRUE), # SD.
                Min = min(percentage, na.rm = TRUE), # Min.
                Max = max(percentage, na.rm = TRUE)) # Max.
  
    # Table transformation to publication-ready report.
    cluster_percentages_cancer2_SUMM_gt<-cluster_percentages_cancer2_SUMM %>%
      mutate(
        Stats = sprintf("%.2f (%.2f) (%.1f, %.1f)", Mean, SD, Min, Max) 
      ) %>%
      select(-Mean, -SD, -Min, -Max) %>% 
      pivot_wider(names_from = ISUP, values_from = Stats)
    
    
    
### DIFFERENTIAL EXPRESSION ANALYSIS ###    
    
    Idents(SCP_data) = "CAF_can_norm" # Differential expression analysis based on Cancer vs CAF vs Non-Cancer cells.
    SCP_data <- subset(SCP_data, CAF_can_norm %in% na.omit(CAF_can_norm)) # Discard NA residual cell group.
    SCP_data$CAF_can_norm <- as.factor(SCP_data$CAF_can_norm) # Factrize covariates.
    SCP_data$Sample_Group <- as.factor(SCP_data$Sample_Group)
    
    
    cancer_caf_markers<-FindAllMarkers(object = SCP_data,
                                       test.use = "MAST", # Use MAST as DE method.
                                       latent.vars = c("orig.ident", "CDR", "Sample_group2", "ISUP"), # Correct by CDR and batch effect.
                                       assay = "RNA",
                                       slot = "data", # Use log-transformed CPMs.
                                       only.pos = FALSE, # Cell markers can only be positive markers.
                                       min.pct = 0.5, # We consider that a marker gene should be expressed in at least 50% of the cells.
                                       logfc.threshold = log2(1.5), # Filter by FC threshold.
                                       verbose = FALSE)
    
    
    cancer_markersall_top<-cancer_caf_markers%>%group_by(., cluster)%>%top_n(n = 50, wt = avg_log2FC) # Select top 50 overexpressed genes in each cell type.
    cancer_markersall_bottom<-cancer_caf_markers%>%slice_min(n = 50, order_by = avg_log2FC) # Select top 50 underexpressed genes in each cell type.
    cancer_caf_txt<-cancer_markersall_top%>%filter(., gene %in% rownames(SCP_data@assays$RNA$scale.data)) # Select those top 50 overexpressed genes also present in the scaled data (genes hvgs, relevant).
    
    # Heatmap of DE genes (Fig 1)
    DoHeatmap(SCP_data, slot = 'scale.data', assay = 'RNA', # Seurat's gene expression heatmap function.
              features = cancer_markersall_bottom$gene, # Plotting top overexpressed genes per cell type.
              group.bar = TRUE,  
              label = FALSE, group.colors = cluster_colors) +
      scale_fill_gradientn(colors = c("#4B79A5", "white", "#D23C44"), 
                           values = scales::rescale(c(-1, 0, 2))) +  # Custom scale. Scaled gene exp. values = 0 are white, negative blue and positive red.
      theme(
            legend.text = element_text(size = 12), 
            legend.title = element_text(size = 16),
            axis.text.y = element_text(size = 8))
    
### FUNCTIONAL ANNOTATION BY ssGSEA (escape) ###
    
    
    Idents(SCP_data) = "CAF_can_norm" # Set cell types to compare (Cancer, CAF and No cancer).
    GS.hallmark <- getGeneSets(library = "H") # Set MiSigDB gene set hallmark H (https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=H).
    
    # Apply runEscape function for ssGSEA analysis.
    SCP_data <- runEscape(SCP_data , 
                          method = "ssGSEA",
                          gene.sets = GS.hallmark, 
                          groups = 1000, # The number of cells to separate the enrichment calculation (default https://bioconductor.org/packages/release/bioc/vignettes/escape/inst/doc/vignette.html).
                          min.size = 5, # Minimum number of gene necessary to perform the enrichment calculation (default).
                          new.assay.name = "escape.ssGSEA") # Slot in assay section in Seurat's object to store de cell x pathway menrichmen score matrix.
    
    # performNormalization() will normalize the enrichment values by calculating the number of genes expressed in each gene set and cell.
    SCP_data <- performNormalization(SCP_data, 
                                     assay = "escape.ssGSEA", 
                                     gene.sets = GS.hallmark, 
                                     scale.factor = SCP_data$nFeature_RNA) # Number of genes as scaling factor (same as guidelines)).
                                     
    # Heatmap of mollecular pathways by normalized enrichment score and grouped by cell type. (Fig 1)
    # Same procedure as guideliness, rest of arguments set as default.
    heatmapEnrichment(SCP_data,
                      group.by = "CAF_can_norm",
                      gene.set.use = 'all', # Selected gene sets to visualize. If "all", the heatmap will be generated across all gene sets.
                      assay = "escape.ssGSEA_normalized", # Use the normalized enrichment score matrix.
                      cluster.rows = TRUE) +
      scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) + # Set different RColorbrewer scale.
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),
            axis.text.y = element_text(size=8))
    
    
    # UMAP of "HALLMARK-EPITHELIAL-MESENCHYMAL-TRANSITION pathway repsentation cell by cell
    DefaultAssay(SCP_data) <- "escape.ssGSEA_normalized" # We'll use the normalized data
    colorblind_vector <- hcl.colors(n=7, palette = "RdYlBu", fixup = TRUE) # Custom color pallete warmer colours represents higher enrochmen scores
    FeaturePlot(SCP_data, 
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
    
    
    
    
    # UMAP of "HALLMARK-ANDROGEN-RESPONSE" pathway repsentation cell by cell
    DefaultAssay(SCP_data) <- "escape.ssGSEA_normalized"
    colorblind_vector <- hcl.colors(n=7, palette = "RdYlBu", fixup = TRUE)
    FeaturePlot(SCP_data, 
                features="HALLMARK-ANDROGEN-RESPONSE") + 
      scale_color_gradientn(colors = rev(colorblind_vector)) + 
      xlab("UMAP 1") +  
      ylab("UMAP 2") +  
      ggtitle("HALLMARK-ANDROGEN-RESPONSE") +  
      theme(axis.title = element_text(size = 18), 
            legend.text = element_text(size = 12), 
            axis.text = element_text(size = 12),
            panel.border = element_rect(colour = "black", fill=NA, size=0.3))+
      theme(strip.text = element_text(size = 18, face = "bold"))