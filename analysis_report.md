# 🧬 Multi-Omics Project: Detailed Analysis Report

This document explains each analytical step taken to process and integrate RNA-seq and methylation data from [GSE224377](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224377), as performed during and after the ABCON 2025 Multi-Omics Workshop.

---


## 📘 Multi-Omics Analysis of Multiple Sclerosis

**Analysis Overview and Task Descriptions**  
*Prepared by:* Mohamed Hussein  
*Date:* 2025-08-10

---

### ✅ Task 0: Setup Instructions & Package Installation

- **Task Title:** Setup Instructions & Package Installation  
- **Rmd File Name:** [00_setup_instructions.Rmd](scripts/00_setup_instructions.Rmd)
- **HTML Output:** [00_setup_instructions.html](https://Mohamed-H-Hussein.github.io/Multi-Omics-Analysis-of-Multiple-Sclerosis-Lesions/00_setup_instruction.html)
- **Input:** None (packages are installed from online sources)  
- **Output:** Installation of required R and Bioconductor libraries  
- **Figures:** None  

**Description:**  
This task provides guidance for installing all R and Bioconductor packages needed to run the analysis.  
The code is initially commented out to prevent automatic execution. Users should remove the `#` to run each line manually.

> **Important Notes:**  
> - If a package is already installed, R will skip reinstalling it.  
> - `install.packages()` is used for CRAN libraries.  
> - `BiocManager::install()` is used for Bioconductor packages.

---

### ✅ Task 1: Load Required Libraries & Setup Environment

- **Task Title:** Load Required Libraries & Setup Environment  
- **Rmd File Name:** [01_libraries_setup.Rmd](scripts/01_libraries_setup.Rmd)
- **HTML Output:** [01_libraries_setup.html](https://Mohamed-H-Hussein.github.io/Multi-Omics-Analysis-of-Multiple-Sclerosis-Lesions/01_libraries_setup.html)
- **Input:** None  
- **Output:** All libraries successfully loaded into the R environment  
- **Figures:** None  

**Description:**  
This task loads all essential libraries required for multi-omics analysis, including RNA-seq processing, methylation analysis, statistical modeling, enrichment, and network construction.

**Loaded Library Categories:**  
- **RNA-seq analysis:** DESeq2, edgeR, vsn  
- **Data wrangling:** tidyverse, dplyr, genefilter  
- **Statistical testing:** multtest, DescTools  
- **Methylation analysis:** minfi, IlluminaHumanMethylation450kanno.ilmn12.hg19  
- **Network analysis:** GeneNet  
- **Functional enrichment:** clusterProfiler, org.Hs.eg.db, AnnotationDbi  
- **Visualization:** ggfortify, factoextra  

✅ The `sessionInfo()` confirms successful loading of all required libraries and environment setup.

---

### ✅ Task 2: RNA-seq Data Preprocessing and Normalization

- **Task Title:** RNA-seq Data Preprocessing and Normalization  
- **Rmd File Name:** [02_RNAseq_Preprocessing.Rmd](scripts/02_RNAseq_Preprocessing.Rmd)  
- **HTML Output:** [02_RNAseq_Preprocessing.html](https://Mohamed-H-Hussein.github.io/Multi-Omics-Analysis-of-Multiple-Sclerosis-Lesions/02_RNAseq_Preprocessing.html) (contains all code, outputs, and figures rendered)  

**Input Files:**  
- [GSE224377_raw_counts_GRCh38.p13_NCBI.tsv](data/GSE224377_raw_counts_GRCh38.p13_NCBI%20(1).tsv
) – raw gene count matrix  
- [metadata-new.csv](data/metadata-new.csv) – sample metadata (Accession, Title, Age, Gender, etc.)  

**Output Files:**  
- [Filtered_counts_less90perc_zeros.csv](results/Filtered_counts_less90perc_zeros.csv)– filtered raw counts  
- [RNA_log2_CPM_TMM_filtered.csv](results/RNA_log2_CPM_TMM_filtered.csv) – TMM-normalized log2(CPM) matrix  
- [RNA_TMM_log2_zscore.csv](results/RNA_TMM_log2_zscore.csv) – z-score normalized expression  
-  [RNA_TMM_log2_minmax_scaled.csv](results/RNA_TMM_log2_minmax_scaled.csv) – min-max scaled matrix  
- PDF plots of histograms after each normalization step  

**Figures:**  
- [histogram_0_raw_counts.pdf](figures/Histogram_0_Raw_Counts.pdf)  
- [histogram_filtered_counts_less90perc_zeros.pdf](figures/Histogram_Filtered_counts_less90perc_zeros.pdf)  
- [histogram_rna_log2_cpm_tmm_filtered.pdf](figures/Histogram_RNA_log2_CPM_TMM_filtered.pdf)  
- [histogram_rna_tmm_log2_zscore.pdf](figures/Histogram_RNA_TMM_log2_zscore.pdf)  
  

**Description:**  
This task processes raw RNA-seq count data through several normalization steps to ensure data quality and comparability across samples. The preprocessing workflow includes:  
- Loading raw counts and sample metadata  
- Matching and aligning sample identifiers  
- Filtering genes with >90% zero counts  
- TMM normalization to adjust for library size bias  
- Log2(CPM + 1) transformation  
- Z-score standardization  
- Min-max scaling to range [-1, 1]  

Each normalization stage includes quality control through histograms to visualize count distributions across samples.

> **Key Notes:**  
> - Filtering genes with many zeros removes low-quality/noisy signals  
> - TMM is a robust normalization method for RNA-seq (from edgeR)  
> - The final normalized matrices will be used in PCA, DEG analysis, and multi-omics integration.  
> - The HTML file generated from this task documents all steps and visualizations interactively.

---

### ✅ Task 3: PCA Analysis of RNA-seq Data

- **Task Title:** Principal Component Analysis (PCA) of RNA-seq Data  
- **Rmd File:** [03_PCA_Analysis.Rmd](scripts/03_PCA_Analysis.Rmd)  
- **HTML Output:** [03_PCA_Analysis.html](https://Mohamed-H-Hussein.github.io/Multi-Omics-Analysis-of-Multiple-Sclerosis-Lesions/03_PCA_Analysis.html)  

**Input Files:**  
- [RNA_TMM_log2_zscore.csv](results/RNA_TMM_log2_zscore.csv) — z-score normalized gene expression matrix (loaded at the beginning if not already in environment)

**Output Files:**  
- [pca_plots_samples.pdf](figures/PCA_plots_samples.pdf) — PC1 vs PC2, PC2 vs PC3, and PC1 vs PC3 plots for samples  
- [pca_plots_genes.pdf](figures/PCA_plots_genes.pdf) — PC1 vs PC2 plot for genes after filtering and imputation  

**Figures:**  
- Sample PCA:  
  - PC1 vs PC2  
  - PC2 vs PC3  
  - PC1 vs PC3  
- Gene PCA:  
  - PC1 vs PC2 (genes after quality filtering and mean imputation)  

**Description:**  
This task applies *Principal Component Analysis (PCA)* to explore the structure and variance in RNA-seq gene expression data.

Two PCA strategies are applied:  
1. **Sample-wise PCA:**  
   - Performed on the z-score normalized expression matrix (genes × samples).  
   - Helps visualize clustering patterns among samples (e.g., by tissue, patient, or condition).  
2. **Gene-wise PCA:**  
   - Performed on the transposed matrix (samples × genes).  
   - Genes with >20% missing values or zero variance are removed.  
   - Remaining missing values are imputed with the gene’s mean.  
   - Reveals gene clusters contributing to biological variation.

> **Key Notes:**  
> - PCA is an essential exploratory step in multi-omics analysis.  
> - It helps identify outliers, batch effects, or grouping patterns.  
> - Outputs guide downstream analyses like DEG, clustering, and multi-omics integration.

---

### ✅ Task 4: Differential Expression & Methylation Analysis

- **Task Title:** Differential Expression and Methylation Analysis  
- **Rmd File:** [04_DEG_DMR_Analysis.Rmd](scripts/04_DEG_DMR_Analysis.Rmd)  
- **HTML Output:** [04_DEG_DMR_Analysis.html](https://Mohamed-H-Hussein.github.io/Multi-Omics-Analysis-of-Multiple-Sclerosis-Lesions/04_DEG_DMR_Analysis.html)  

**Input Files:**  
- dge object (from filtered, normalized RNA-seq)  
- [metadata-new.csv](data/metadata-new.csv) (sample info)  
- [Methylation_20000_Top_expressed.csv](data/Methylation_20000_Top_expressed.csv) (preprocessed methylation matrix)  

**Output Files:**  

- **DEG:**  
  - [DEG_results_all.csv](results/DEG_results_all.csv) – All differentially expressed genes  
  - [DEG_significant_filtered.csv](results/DEG_significant_filtered.csv) – Significant genes (FDR < 0.05, |log2FC| > 1)  

- **DMR:**  
  - [dmrs_methylation_lm_based.csv](results/dmrs_methylation_lm_based.csv) – Significant DMRs via linear model  
  - [dmrs_methylation.csv](results/dmrs_methylation.csv) – Significant DMRs via t-test  
  - [dmrs_glm_filtered.csv](results/dmrs_glm_filtered.csv) – Significant DMRs via GLM  

**Description:**  
This task performs differential expression and methylation analysis to identify biologically meaningful differences between *lesion* and *NAWM* samples:

**Differential Expression (RNA-seq):**  
- Uses *edgeR* to model count data with negative binomial GLM  
- Estimates dispersions and fits model per group  
- Contrast tested: lesion vs NAWM  
- Genes filtered by FDR < 0.05 and |log2 Fold Change| > 1  

**Differential Methylation:**  
- Matrix contains beta values per CpG site  
- Groups defined as NAWM (control) vs Lesion (case)  
- Three methods applied:  
  - lm() linear regression  
  - rowttests() row-wise t-test  
  - glm() generalized linear model (Gaussian)  

> **Key Notes:**  
> - Only CpGs with complete data were analyzed.  
> - Results from all three methylation methods saved for comparison.  
> - Filter thresholds consistent (FDR < 0.05 and logFC > log2(1)).

---

### ✅ Task 5: DEG–DMR Correlation & Linear Regression (with Age + Gender)

- **Task Title:** DEG–DMR Correlation & Linear Regression (with Age + Gender)  
- **Rmd File:** [05_Integration_Regression.Rmd](scripts/05_Integration_Regression.Rmd)  
- **HTML Output:** [05_Integration_Regression.html](https://Mohamed-H-Hussein.github.io/Multi-Omics-Analysis-of-Multiple-Sclerosis-Lesions/05_Integration_Regression.html)

**Input Files:**  
- [RNA_TMM_log2_zscore.csv](results/RNA_TMM_log2_zscore.csv) – z-score normalized RNA-seq matrix  
- deg_filtered (from Task 4) – significant DEGs  
- [Methylation_20000_Top_expressed.csv](data/Methylation_20000_Top_expressed.csv) – methylation matrix  
- [dmrs_methylation_lm_based.csv](results/dmrs_methylation_lm_based.csv) – significant DMRs (lm method)  
- [metadata-new.csv](data/metadata-new.csv) – sample metadata (age, gender info)  

**Output Files:**  
- [DEG_DMR_linear_models_withAgeGender.csv](results/DEG_DMR_linear_models_withAgeGender.csv) – Regression results for each DEG–DMR pair adjusted for covariates  

**Description:**  
This task integrates differential gene expression and methylation data by performing multiple linear regression analyses adjusting for age and gender. The goal is to detect significant associations between DEGs and DMRs in matched samples.

**Processing Steps:**  
- RNA matrix filtered for significant DEGs  
- Methylation matrix filtered for DMRs  
- Samples matched across matrices and metadata  
- Random sampling of 25 genes and 100 CpGs (set.seed(42))  
- Linear model per gene–CpG pair:  
  expression ~ methylation + age + gender  
- Extract coefficients, p-values, and R²  

> **Key Notes:**  
> - Missing metadata values should be handled before regression  
> - Results contain coefficients and significance for methylation, age, gender  
> - Output can be extended or used for downstream visualization

---

### ✅ Task 6: Network Construction + Annotation

- **Task Title:** Network Construction and Annotation  
- **Rmd File:** [06_Network_Construction.Rmd](scripts/06_Network_Construction.Rmd)  
- **HTML Output:** [06_Network_Construction.html](https://Mohamed-H-Hussein.github.io/Multi-Omics-Analysis-of-Multiple-Sclerosis-Lesions/06_Network_Construction.html)  

**Input Files:**  
- [RNA_TMM_log2_zscore.csv](results/RNA_TMM_log2_zscore.csv) – Z-score normalized expression matrix for DEGs  
- [DEG_significant_filtered.csv](results/DEG_significant_filtered.csv) – Significant DEGs  
- [Methylation_20000_Top_expressed.csv](data/Methylation_20000_Top_expressed.csv) – Top methylation matrix  
- [dmrs_methylation_lm_based.csv](results/dmrs_methylation_lm_based.csv) – Significant DMRs  
- [DEG_DMR_linear_models_withAgeGender.csv](results/DEG_DMR_linear_models_withAgeGender.csv) – DEG–DMR regression results  

**Output Files:**  

- **Combined Network:**  
  - [full_multiomics_network.csv](results/full_multiomics_network.csv) – Combined edge list (Gene-Gene, CpG-CpG, CpG-Gene)  

- **Annotated Edges:**  
  - [gene_edges_annotated.csv](results/gene_edges_annotated.csv) – Gene-Gene edges  
  - [cpg_edges_annotated.csv](results/cpg_edges_annotated.csv) – CpG-CpG edges  
  - [gene-cpg_edges_annotated.csv](results/gene-cpg_edges_annotated.csv) – CpG-Gene regression edges  

**Description:**  
Constructs a comprehensive multi-omics network integrating gene expression and DNA methylation data.

- **Gene–Gene Network:**  
  - Pearson and partial correlations via GeneNet  
  - Significant edges p < 0.05  
  - Gene symbols annotated via org.Hs.eg.db  

- **CpG–CpG Network:**  
  - Partial correlations on top 1000 CpGs from DMRs  
  - Significant edges filtered  
  - CpGs annotated with IlluminaHumanMethylation450kanno.ilmn12.hg19  

- **CpG–Gene Edges (Regression):**  
  - From DEG ~ methylation + age + gender regression  
  - Edges weighted by -log10(p-value)  
  - Annotations combine CpG and gene symbols  

> **Key Notes:**  
> - Final network includes Gene-Gene, CpG-CpG, and CpG-Gene edges  
> - Annotated for biological interpretation  
> - Visualizable in Cytoscape or igraph  

---

### ✅ Task 7: GO Enrichment Analysis of DEGs and DMR-Associated Genes

- **Task Title:** Functional Enrichment Analysis  
- **Rmd File:** [07_GO_Enrichment_Analysis.Rmd](scripts/07_GO_Enrichment_Analysis.Rmd)  
- **HTML Output:** [07_GO_Enrichment_Analysis.html](https://Mohamed-H-Hussein.github.io/Multi-Omics-Analysis-of-Multiple-Sclerosis-Lesions/07_Enrichment_In-terpretation.html)  

**Input Files:**  
- [DEG_significant_filtered.csv](results/DEG_significant_filtered.csv) – DEGs  
- [dmrs_methylation_lm_based.csv](results/dmrs_methylation_lm_based.csv) – DMRs  
- Illumina 450k annotation

**Output Files:**  
- Dot plots for DEG and DMR-associated gene GO enrichment (in HTML report)

- **Figures:**  
  - [deg_go_enrichment.png](figures/deg_go_enrichment.png)  
  - [dmr_associated_gene_go_enrichment.png](figures/dmr_associated_gene_go_enrichment.png)
  

**Description:**  
Performs GO enrichment analysis focusing on Biological Process (BP) ontology for DEGs and DMR-associated genes.

- **DEG Enrichment:**  
  - enrichGO() with ENTREZ IDs  
  - FDR adjustment  
  - Visualized with dotplot()  

- **DMR Gene Mapping + Enrichment:**  
  - Annotate CpGs to genes using Illumina 450k annotation  
  - Convert gene symbols to ENTREZ IDs  
  - GO enrichment and dot plot visualization  

> **Key Notes:**  
> - DEG and DMR-associated gene sets analyzed separately  
> - BP category selected for function/pathways  
> - Dot plots summarize top enriched terms  
> - Supports biological interpretation  

---

### ✅ Task 8: Multi-Omics Integration with DIABLO

- **Task Title:** Multi-Omics DIABLO Modeling  
- **Rmd File:** [08_DIABLO-Model.Rmd](scripts/08_DIABLO-Model.Rmd)  
- **HTML Output:** [08_DIABLO-Model.html](https://Mohamed-H-Hussein.github.io/Multi-Omics-Analysis-of-Multiple-Sclerosis-Lesions/08_DIABLO-Model.html)  

**Input Files:**  
- [RNA_TMM_log2_zscore.csv](results/RNA_TMM_log2_zscore.csv) – Gene expression  
- [Methylation_20000_Top_expressed.csv](data/Methylation_20000_Top_expressed.csv) – Methylation  
- [metadata-new.csv](data/metadata-new.csv) – Sample metadata  

**Output Files:**  
- [selected_genes_component1.csv](results/selected_genes_component1.csv), [selected_genes_component2.csv](results/selected_genes_component2.csv)  
- [selected_cpgs_component1.csv](results/selected_cpgs_component1.csv), [selected_cpgs_component2.csv](results/selected_cpgs_component2.csv)  
- DIABLO plots (sample projections, circos plots, loading plots)in html report  
- Model performance plots (AUC, confusion matrix, error rate) in html report 

**Figures:**  
- [diablo_plot.png](figures/diablo_plot.png)  
- [diablo_circos_plot.png](figures/diablo_circos_plot.png)  
- [diablo_plot_rna_loadings.png](figures/diablo_plot_rna_loadings.png)  
- [diablo_plot_methylation_loadings.png](figures/diablo_plot_methylation_loadings.png)  
- [roc_curve_1.png](figures/roc_curve_1.png)  
- [roc_curve_2.png](figures/roc_curve_2.png)  
- [peformance_diablo.png](figures/performance_diablo.png)  

**Description:**  
Implements DIABLO multi-omics integration to link gene expression and methylation features distinguishing lesion vs NAWM samples.

**Preprocessing:**  
- Filter, normalize, transform RNA and methylation data  
- Match samples across omics  

**Modeling:**  
- Construct data blocks  
- Outcome: sample.type (lesion vs NAWM)  
- Remove near-zero variance features  
- Hyperparameter tuning for keepX  
- Final model with 2 components  

**Visualizations:**  
- ROC curves, sample projections, circos, loading plots  

**Model Performance:**  
- Accuracy: 0.85  
- Confusion matrix:  

```
Predicted  lesion  NAWM
Actual
lesion     7     0
NAWM       2     5

```


**Biological Insight:**  
- Identified gene and CpG markers driving lesion vs NAWM  
- Useful for enrichment and network analysis  

> **Key Notes:**  
> - `set.seed(42)` for reproducibility

