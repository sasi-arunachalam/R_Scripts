
# load libraries
#install.packages('Seurat')
#install.packages("hdf5r")
library(Seurat)
library(hdf5r)
library(tidyverse)

# Load the test dataset
test.sparse.m <- Read10X_h5(filename = 'testdata.h5')
str(test.sparse.m)
cts <-  test.sparse.m$`Gene Expression`
# Initialize the Seurat object with the raw (non-normalized data).
test.seurat.obj <- CreateSeuratObject(counts = cts, project = "test", min.cells = 3, min.features = 200)
str(test.seurat.obj)
test.seurat.obj
# 29552 features across 42081 samples

# 1. QC -------
View(test.seurat.obj@meta.data)
# % MT reads
test.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(test.seurat.obj, pattern = "^MT-")
View(test.seurat.obj@meta.data)

VlnPlot(test.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(test.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# 2. Filtering -----------------
test.seurat.obj <- subset(test.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)
# 3. Normalize data ----------
#test.seurat.obj <- NormalizeData(test.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
test.seurat.obj <- NormalizeData(test.seurat.obj)
str(test.seurat.obj)


# 4. Identify highly variable features --------------
test.seurat.obj <- FindVariableFeatures(test.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(test.seurat.obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(test.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


# 5. Scaling -------------
all.genes <- rownames(test.seurat.obj)
test.seurat.obj <- ScaleData(test.seurat.obj, features = all.genes)

str(test.seurat.obj)

# 6. Perform Linear dimensionality reduction --------------
test.seurat.obj <- RunPCA(test.seurat.obj, features = VariableFeatures(object = test.seurat.obj))

# visualize PCA results
print(test.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(test.seurat.obj, dims = 1, cells = 500, balanced = TRUE)

# determine dimensionality of the data
ElbowPlot(test.seurat.obj)

# 7. Clustering ------------
test.seurat.obj <- FindNeighbors(test.seurat.obj, dims = 1:15)

# understanding resolution
test.seurat.obj <- FindClusters(test.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(test.seurat.obj@meta.data)

DimPlot(test.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)

# setting identity of clusters
Idents(test.seurat.obj)
Idents(test.seurat.obj) <- "RNA_snn_res.0.1"
Idents(test.seurat.obj)

# non-linear dimensionality reduction --------------
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
test.seurat.obj <- RunUMAP(test.seurat.obj, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(test.seurat.obj, reduction = "umap")
