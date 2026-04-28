# Single-cell RNA-seq Analysis of CD8+ T cells and Monocytes/Macrophages

## Overview

This repository contains the code used for the analysis of single-cell RNA-seq data comparing **WT vs Casp1-/- conditions**, focusing on:

- **CD8+ T cells**
- **Monocytes/Macrophages (Mo/Mac)**

This code is associated with a **manuscript currently in preparation / under submission**.

The pipeline includes:

- Quality control and filtering
- Doublet removal (scDblFinder)
- Data integration (Seurat v5)
- Cell type annotation using canonical markers
- Subsetting of CD8+ and Mo/Mac populations
- Pseudobulk differential expression analysis (DESeq2)
- Functional enrichment analysis (GSEA, GO Biological Processes)
- Visualization (UMAP, violin plots, heatmaps, MA plots, Venn diagrams)

---

## Requirements

- **R 4.5.0**
- Platform: x86_64-pc-linux-gnu  
- OS: Ubuntu 20.04.6 LTS  

---

## R package versions (critical for reproducibility)

### Single-cell
- Seurat (5.4.0)  
- SeuratObject (5.0.2)  
- SingleCellExperiment (1.26.0)  
- scDblFinder  

### Differential expression
- DESeq2 (1.44.0)  
- apeglm  

### Functional enrichment
- clusterProfiler  
- enrichplot (1.24.0)  
- org.Mm.eg.db  
- GO.db (3.19.1)  

### Data manipulation
- dplyr (1.2.1)  
- tidyr (1.3.1)  
- tibble (3.2.1)  
- readr  

### Visualization
- ggplot2 (3.5.1)  
- patchwork (1.3.0)  
- ggpubr  
- ggrepel (0.9.6)  
- pheatmap  
- ggVennDiagram  
- circlize  
- RColorBrewer  

---

## Workflow summary

1. Preprocessing and quality control  
2. Doublet removal (per sample)  
3. Data integration (Seurat)  
4. Clustering and annotation  
5. Subsetting CD8+ and Mo/Mac populations  
6. Pseudobulk differential expression (DESeq2)  
7. GSEA (GO Biological Processes)  
8. Visualization and figure generation  

---

## Usage

Run the main analysis script:

```r
source("script/analysis.R")
```

Input data paths must be adapted to the local environment.

---

## Reproducibility

All analyses were performed using the R session described above.

The full session information is available in:

- `sessionInfo.txt`

---

## Output

The pipeline generates:

- UMAP visualizations  
- Differential expression results  
- Heatmaps (CD8+ and Mo/Mac)  
- GSEA results (GO BP)  
- Publication-ready figures  

---

## Citation

If you use this code, please cite:

**Huang et al. (2026). SCaspase-1 regulates TBK1 and PPARα to suppress steatosis-associated
hepatocarcinogenesis. Manuscript in preparation.**

The repository will be updated with publication details upon acceptance.

---

## Notes

- Raw data is not included  
- File paths should be adapted to local environment  
- Scripts require preprocessed Seurat objects  

---

 


