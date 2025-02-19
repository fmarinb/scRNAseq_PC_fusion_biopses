# Single-Cell transcriptomics to characterize the tumor microenvironment of prostate cancer fusion biopsies

This repository stores the R scripts corresponding to the secondary scRNAseq analyses on data already pre-processed by Genyo's Genomics and Bioinformatics Unit. The scripts contain the steps used for:

1. **Cancer-oriented cell annotation** by *scATOMIC*.
2. **Differential expression analysis** between clusters of cells.
3. **Functional enrichment analysis** by *ssGSEA*.
4. **Differential abundance analysis** to identify cells associated with high tumourgenicity.
5. **Training of machine learning models** for the classification of genes of interest in pseudobulk format based on the tumourgenicity of the samples.
6. **Study of gene importance** in the best-performing model (*SHAP values*) and **exploratory analysis of gene expression**.

## One Sentence Summary  
Single-cell analysis of fusion biopsies identified prostate cancer biomarkers and tumor-specific cell populations with distinct genetic signatures.

## Abstract  

Prostate cancer (PC) is the second most common malignancy worldwide and the leading cancer in men. Despite advances in treatment, resistance mechanisms and variability in the tumor microenvironment remain significant challenges. Single-cell RNA sequencing (scRNA-seq) enables the characterization of tumor microenvironment complexity, while transperineal fusion-guided biopsy enhances precision in PC diagnostics. However, studies using scRNA-seq from fusion-guided biopsy samples at diagnosis are limited. In the present study, we analyzed the tumor microenvironment of 31 fusion-guided biopsy samples from 18 newly diagnosed PC patients. Samples were taken from distal locations and directed to the suspected tumor site and subsequently reclassified into high- or low-tumorigenicity samples. Using scRNA-seq specific tools, we identified cancer cells and cancer-associated fibroblasts (CAFs) enriched with PC-associated signatures. Cancer cells overexpressed *TMPRSS2*, *KLK2-4*, *NKX3-1*, and *HOXB13*, while CAFs displayed high expression of *COL* family, *CALD1*, *TPM2*, and *MYL9*. Single-cell pathway enrichment analyses revealed androgen response activation in cancer cells and epithelial-to-mesenchymal transition in CAFs. Differential abundance analyses identified cancer cell subpopulations enriched in high tumorigenicity samples, characterized by overexpression of *PCA3*, *SCHLAP1*, *CRISP3*, and *ERG*. These genes were used to develop machine learning models, with Elastic Net achieving moderate to high performance. SHAP analysis of the resulting model highlighted the discriminatory power of genes like *PCA3*, *TNFSF10*, *CPNE4*, and the novel *ENSG00000289332* long non-coding RNA. Bioinformatic analysis identifies key PC biomarkers and characterizes the tumor microenvironment, overcoming challenges such as low cell viability, high tumor heterogeneity, and limited cohort size using fusion-guided biopsy.
