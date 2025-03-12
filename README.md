# GSE65194 HER2 & TNBC Breast Cancer Analysis

## About

This project focuses on analyzing HER2+ and Triple-Negative Breast Cancer (TNBC) gene expression data using the GSE65194 dataset. It involves differential expression analysis (DEA), identification of key biomarkers, pathway enrichment analysis and machine learning models for classification. The goal is to uncover potential therapeutic targets and gain insights into drug resistance mechanisms in these aggressive breast cancer subtypes.

## Project Objectives

- Perform **differential gene expression analysis (DEA)** to identify significantly upregulated and downregulated genes in HER2+ and TNBC samples.
- Conduct **principal component analysis (PCA)** to explore gene expression variance.
- Identify key **biomarkers** associated with HER2+ and TNBC.
- Perform **pathway enrichment analysis** to understand affected biological pathways.
- Implement **machine learning models** (Random Forest, SVM, Logistic Regression, Neural Networks) to classify breast cancer subtypes based on gene expression.
- Develop visualization techniques such as **volcano plots, heatmaps, and PCA plots**.
- Future goals include integrating multi-omics data and deploying a **bioinformatics web application**.

## Dataset

- **GSE65194**: A publicly available gene expression dataset from the Gene Expression Omnibus (GEO), containing transcriptomic profiles of HER2+ and TNBC breast cancer samples.

## Tools & Technologies Used

- **Programming Language:** R (Primary)
- **Libraries & Packages:**
  - `Bioconductor` (DESeq2, edgeR, limma, clusterProfiler, TCGAbiolinks, ggplot2)
  - `caret`, `randomForest`, `e1071` (for machine learning models)
  - `pheatmap`, `ggpubr`, `ggfortify` (for visualization)
- **Environment:** RStudio

## Analysis Workflow

1. **Data Preprocessing:**
   - Load and normalize gene expression data.
   - Remove low-expression genes and normalize counts.
2. **Exploratory Data Analysis (EDA):**
   - PCA for variance analysis.
   - Visualize sample distribution.
3. **Differential Expression Analysis (DEA):**
   - Identify differentially expressed genes (DEGs) using DESeq2.
   - Generate volcano plots and heatmaps.
4. **Pathway Enrichment Analysis:**
   - Gene Ontology (GO) and KEGG pathway analysis.
5. **Machine Learning Models:**
   - Train classifiers (Random Forest, SVM, etc.) to predict HER2+ vs TNBC.
6. **Result Interpretation & Future Directions:**
   - Identify key biomarkers.
   - Discuss clinical relevance and future improvements.

## Future Work

- Expand ML models with deep learning approaches.
- Integrate clinical metadata for better predictive modeling.
- Develop a web-based interactive tool using Flask/Django.
- Scale the project with Docker/Kubernetes for deployment.

## Repository Structure

```
📂 HER2_TNBC_Breast_Cancer_Analysis
├── 📁 data/                     # Raw & processed datasets
├── 📁 scripts/                  # R scripts for analysis
├── 📁 results/                  # Outputs (plots, tables, models)
├── 📄 README.md                 # Project documentation
├── 📄 analysis.R                 # Main analysis script
```

## How to Run

1. Clone the repository:
   ```bash
   git clone https://github.com/tahagill/Bioinformatics_GSE65194_Breast_Cancer_Resistance.git
   ```
2. Open RStudio and install required packages:
   ```r
   install.packages(c("ggplot2", "pheatmap", "caret", "randomForest", "e1071"))
   BiocManager::install(c("DESeq2", "edgeR", "limma", "clusterProfiler", "TCGAbiolinks"))
   ```
3. Run the analysis script:
   ```r
   source("scripts/optimus.R")
   ```

##
