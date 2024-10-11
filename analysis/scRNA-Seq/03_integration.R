setwd("#PATH#/data/hcc/10xdropseq-inte")
library(dplyr)
library(patchwork)
library(ggplot2)
library(Seurat)

# load pre-processed data from each sample
hcc1 <- readRDS("#PATH#/data/hcc/hcc1/10x/raw.filter.rds")
hcc2 <- readRDS("#PATH#/data/hcc/hcc2/10x/raw.filter.rds")
hcc3 <- readRDS("#PATH#/data/hcc/hcc3/drop-seq/raw.filter.rds")
hcc4 <- readRDS("#PATH#/data/hcc/hcc4/drop-seq/raw.filter.rds")
hcc5 <- readRDS("#PATH#/data/hcc/hcc5/drop-seq/raw.filter.rds")
hcc6 <- readRDS("#PATH#/data/hcc/hcc6/drop-seq/raw.filter.rds")
hcc7 <- readRDS("#PATH#/data/hcc/hcc7/drop-seq/raw.filter.rds")
hcc8 <- readRDS("#PATH#/data/hcc/hcc8/10x/raw.filter.rds")
hcc9 <- readRDS("#PATH#/data/hcc/hcc9/10x/raw.filter.rds")

# merge all the data
spl.name <- c("HCC1", "HCC2", "HCC3", "HCC4", "HCC5", "HCC6", "HCC7", "HCC8", "HCC9")
data <- merge(x = hcc1, y = c(hcc2, hcc3, hcc4, hcc5, hcc6, hcc7, hcc8, hcc9), add.cell.ids = spl.name, project = "hcc-inte")

# data normalization for each sample
data.list <- SplitObject(data, split.by = "patient")
rm(hcc1, hcc2, hcc3, hcc4, hcc5, hcc6, hcc7, hcc8, hcc9, data)
data.list <- lapply(X = data.list, FUN = function(x){
  x <- NormalizeData(object = x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(object = x, selection.method = "vst", nfeatures = 3000)
})
# add patient information
for(i in 1:length(spl.name)){
  data.list[[i]]@meta.data$patient <- spl.name[i]
}

# data integration using 10x data as reference
ref_dataset <- which(names(data.list) %in% c("HCC1", "HCC2", "HCC8", "HCC9"))
data.anchors <- FindIntegrationAnchors(data.list, normalization.method = "LogNormalize", anchor.features = 3000, reference = ref_dataset, dims = 1:40)
data.combined <- IntegrateData(anchorset = data.anchors, new.assay.name = "integrated", normalization.method = "LogNormalize", dims = 1:40)

# run the standard seurat integration workflow using CCA method
var.genes <- VariableFeatures(data.combined)
var.genes <- var.genes[!(grepl(pattern = "^MT-", x = var.genes))]
data.combined <- ScaleData(data.combined, features = var.genes, vars.to.regress = "percent.mt")
data.combined <- RunPCA(data.combined, features = var.genes, npcs = 50)
rm(data.list, ref_dataset, data.anchors)

# non-linear dimensional reduction
ElbowPlot(data.combined, ndims = 50)
data.combined <- RunUMAP(data.combined, dims = 1:30, reduction = "pca")
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindClusters(data.combined, resolution = 1)

# assign new idents for each cluster
new.cluster.ids <- c("CD4+ memory", "Macrophage", "Macrophage", "CD8+ exhausted",
                     "Macrophage", "Endothelial cell", "Dendritic cell", "CD8+ memory",
                     "NK", "CD4+ Treg", "HPC", "Unknown", "NK", "HPC", "Dendritic cell", 
                     "B cell", "HPC", "CD8+ cytotoxic", "Endothelial cell", "Fibroblast", 
                     "Proliferative T", "CD8+ cytotoxic", "Endothelial cell", "HPC", 
                     "CD4+ memory", "Macrophage", "PlasmaB cell", "CD8+ memory", "HPC",
                     "Unknown", "Dendritic cell", "Neutrophil", "Unknown", "Dendritic cell",
                     "Unknown", "Dendritic cell", "Mast cell", "Unknown")
names(new.cluster.ids) <- levels(data.combined)
data.combined <- RenameIdents(data.combined, new.cluster.ids)

# umap showing clustering of all cell types
DimPlot(data.combined, reduction = "umap", label = TRUE, repel = TRUE)

## note
# 1. remove all the Unknown cells and non-immune cells(HPC, Fibroblast, Endothelial cell) with expression of PTPRC > 1
# 2. repeat data processing above for the second round of clustering
