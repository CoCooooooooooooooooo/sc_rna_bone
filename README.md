# Single-cell RNA-seq Analysis: Bone Defect Repair Model

## Overview

Single-cell transcriptomic analysis of rat bone defect tissue across three treatment groups and two time points, performed using Python/Scanpy.

| | Group A | Group B | Group C |
|---|---|---|---|
| **Description** | Defect control | Reference material | Experimental material |
| **Role** | Negative control | Positive control | Focus of study |

**Time points**: 1 week and 2 weeks post-implantation  
**Sequencing platform**: Singleron scRNA-seq  
**Total cells analyzed**: 86,233  
**Genes detected**: 25,559  

---

## Repository Structure

```
├── analysis_fixed.py                   # Main preprocessing pipeline
├── annotate.py                         # Cell type annotation
├── marker_genes.py                     # Marker gene identification
├── comparison.py                       # Cross-group comparison analysis
├── download_fixed.sh                   # Data download script
├── cell_proportions.csv                # Cell type proportions per group
├── proportion_differences.csv          # Pairwise proportion differences
├── marker_genes.csv                    # Top marker genes per cluster
├── marker_genes_full.csv               # Full marker gene statistics
├── qc_stats.csv                        # QC metrics summary
├── cluster_sample_composition.csv      # Cluster-sample composition
└── figures/                            # All visualizations
    ├── violin_qc.png
    ├── pca_variance_ratio_pca.png
    ├── umap_clusters.png
    ├── umap_samples.png
    ├── umap_annotated.png
    ├── umap_annotated_legend.png
    ├── dotplot__markers_dotplot.png
    ├── heatmap__markers_heatmap.png
    └── comparison/
        ├── heatmap_all.png
        ├── stacked_all.png
        ├── CvsA_1w.png
        ├── CvsA_2w.png
        ├── CvsA_combined.png
        ├── CvsA_stacked.png
        ├── CvsB_1w.png
        ├── CvsB_2w.png
        ├── CvsB_combined.png
        ├── CvsB_stacked.png
        ├── C_temporal.png
        └── C_temporal_stacked.png
```

---

## Analysis Pipeline

### Step 1 — Preprocessing (`analysis_fixed.py`)

- Read 10x-format matrices for 6 samples (A, B, C, A_2w, B_2w, C_2w)
- Merge into a single AnnData object
- **QC filtering**: min_genes=200, max_genes=6000, mitochondrial gene pct < 20%
- **Normalization**: normalize_total (target=10,000) + log1p
- **Highly variable gene selection** with batch_key='sample'
- **Dimensionality reduction**: PCA (50 PCs) → KNN graph → UMAP
- **Clustering**: Leiden algorithm, resolution=0.5 → 20 clusters

### Step 2 — Cell Type Annotation (`annotate.py`)

- Wilcoxon rank-sum test to identify marker genes per cluster
- Manual annotation based on known rat bone tissue markers
- 19 cell types identified:

| Category | Cell Types |
|----------|-----------|
| Stromal | Fibroblast, Stromal/Perivascular cell, Endothelial cell |
| Bone/Cartilage | Chondrocyte, Chondrocyte-related |
| Myeloid immune | Neutrophil, Mature Neutrophil, Neutrophil (defensin+), Neutrophil/Monocyte, Macrophage/DC, Antigen-presenting cell, Dendritic cell, Mast cell |
| Lymphoid | T cell, B cell/DC, B cell precursor |
| Erythroid | Erythrocyte, Erythroid precursor |
| Proliferating | Proliferating cell |

### Step 3 — Marker Gene Analysis (`marker_genes.py`)

- Wilcoxon rank-sum test, top 50 genes per cluster
- Outputs: dotplot, heatmap, CSV tables

### Step 4 — Cross-group Comparison (`comparison.py`)

- Cell type proportion calculation per sample
- Three comparisons:
  - **C vs A**: Experimental material vs defect control
  - **C vs B**: Experimental material vs reference material
  - **C temporal**: C group 1w → 2w dynamics
- Outputs: stacked bar plots, bar charts, heatmaps, proportion difference tables

---

## Key Findings

### C vs B — Biocompatibility

Group B showed **46% neutrophil infiltration** at 1 week vs only **5.3% in Group C**, indicating significantly superior biocompatibility of the experimental material. By 2 weeks, B group neutrophils remained at 19.2% vs 2.0% in C.

### C vs A — Early Repair Activation

Chondrocyte-related cells were **2× higher in C at 1 week** (10.8% vs 5.3%), suggesting earlier recruitment of bone progenitor cells. Proliferating cells were also elevated in C (6.3% vs 3.9%).

### C Temporal Dynamics (1w → 2w)

| Cell Type | 1 week | 2 weeks | Trend |
|-----------|--------|---------|-------|
| Neutrophil | 5.3% | 2.0% | ↓ inflammation resolving |
| Macrophage/DC | 6.4% | 3.5% | ↓ inflammation resolving |
| Fibroblast | 17.9% | 23.4% | ↑ matrix remodeling |
| Erythroid precursor | 3.8% | 6.0% | ↑ vascularization |

The temporal pattern is consistent with normal bone healing progression: acute inflammation at 1w transitioning to matrix remodeling at 2w.

---

## Environment Setup

```bash
python3 -m venv scanpy_env
source scanpy_env/bin/activate
pip install scanpy leidenalg python-igraph matplotlib seaborn anndata
```

**Key package versions used**:
- Python 3.11
- Scanpy 1.9.x
- AnnData 0.9.x

---

## Data Availability

Raw data (.loom, matrix.tar files) and processed AnnData objects (.h5ad, >500MB) are not included in this repository due to file size limits and data sharing restrictions (unpublished).

**Data source**: Singleron VE platform, Project P26032704  
**Species**: Rattus norvegicus  
**Tissue**: Bone defect repair tissue  
