# ScRNA-seq T Cell Trajectory Analysis

This repository presents a detailed and reproducible pipeline for analyzing single-cell RNA-seq data to study **T cell differentiation** using R. It leverages the power of **Seurat** for data preprocessing and clustering and **Monocle3** for trajectory inference and pseudotime analysis. The analysis walks through quality control, clustering, cell type annotation, and dynamic modeling of cell state transitions.

## Key Features

- Preprocessing and normalization of UMI count data
- Quality control filtering based on RNA counts, feature counts, and mitochondrial content
- Dimensionality reduction (PCA, UMAP)
- Clustering of cell populations using graph-based methods
- Annotation of clusters using known marker genes
- Conversion of Seurat object to Monocle3 format
- Trajectory graph construction and pseudotime inference
- Visualization of gene expression trends along pseudotime


## Repository Structure

```
 Tcell-Trajectory-Analysis/
├── ScRNA_Tcell_Trajectory_Analysis.Rmd    # Complete analysis pipeline (R Markdown)
├── ABC_umi_matrix.csv                     # UMI count matrix (genes × cells)
├── ABC_Meta.txt                           # Cell-level metadata (nCount_RNA, nFeature_RNA, cluster info)
├── ABC_Marker.txt                         # Marker genes used for manual annotation
├── README.md                              # Project documentation
```

##  Analysis Workflow

### 1. Load and Prepare Data

- Load expression matrix, cell metadata, and marker list
- Transpose expression matrix and create Seurat object

### 2. Quality Control

- Filter cells based on total counts, gene features, and mitochondrial gene percentage using `PercentageFeatureSet()`
- Example:
```r
seu.obj.filtered <- subset(seu.obj, subset = nCount_RNA > 800 & nFeature_RNA > 500 & mitopercent < 20)
```

### 3. Normalization and Feature Selection

- Normalize data using `NormalizeData()`
- Identify highly variable features with `FindVariableFeatures()`

### 4. Dimensionality Reduction and Clustering

- Perform PCA: `RunPCA()`
- Run UMAP: `RunUMAP()`
- Construct SNN graph and cluster cells with `FindNeighbors()` and `FindClusters()`

### 5. Cell Type Annotation

- Match cluster marker genes with known gene signatures from `ABC_Marker.txt`
- Assign cell type labels based on marker gene expression

### 6. Trajectory Inference with Monocle3

- Convert Seurat object to Monocle3 CDS using `as.cell_data_set()`
- Reduce dimensions with `reduce_dimension()`
- Learn trajectory graph and order cells with `learn_graph()` and `order_cells()`

### 7. Visualization

- UMAP with cluster and cell type labels
- Pseudotime trajectories
- Gene expression trends along pseudotime
- Feature plots for selected marker genes

## Sample Outputs

- `UMAP`: Cluster and cell-type labeled plots
- `Trajectory`: Lineage tree showing progression and branches
- `Pseudotime Heatmap`: Dynamic changes in gene expression across cell states
- `Gene Trends`: Expression of lineage-driving genes over pseudotime
