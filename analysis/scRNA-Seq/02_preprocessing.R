setwd("#PATH#/data/hcc/hcc1/10x")
library(ggplot2)
library(Seurat)
library(dplyr)
library(Matrix)
library(DoubletFinder)

# sample information
sample <- "HCC1"#example
region <- c("NT", "PT2", "PT4", "PT5")
group <- c("NT", "P", "P", "S")
data.dir <- "#PATH#/data/hcc/hcc1/10x/outs/count/filtered_feature_bc_matrix"

# load raw count data
raw.matrix <- Read10X(data.dir = data.dir, gene.column = 2, unique.features = TRUE)
# remove cells with high umi count
total.umi <- colSums(raw.matrix)
umi.cutoff <- mean(total.umi) + 3*sd(total.umi)
raw.matrix <- raw.matrix[, total.umi <= umi.cutoff]
# remove cells with low umi count
raw.matrix <- raw.matrix[, colSums(raw.matrix) >= 800]
# create the Seurat object
raw <- CreateSeuratObject(counts = raw.matrix, project = sample, assay = "RNA", min.cells = 0, min.features = 0)
rm(raw.matrix)
# add region information
raw$orig.ident <- as.character(raw$orig.ident)
for(i in 1:length(region)){
  raw@meta.data[grepl(pattern = i, x = Cells(raw)), "orig.ident"] <- region[i]
}
# add group information
raw$group <- "NT"
for(i in 1:length(region)){
  raw@meta.data[raw$orig.ident == region[i], "group"] <- group[i]
}
# percentage of mitochondrial genes
raw$percent.mt <- PercentageFeatureSet(raw, pattern = "^MT-")

# filter low quality cells (genes and percent.mt)
raw.filter <- subset(raw, nFeature_RNA >= 300 & percent.mt <= 20)
# filter genes
gene.min.cell <- 30 # for drop-seq data: gene.min.cell <- 10
raw.filter <- raw.filter[rowSums(raw[['RNA']]@counts > 0) >= gene.min.cell,]

# add library and patient information
raw.filter$lib.method <- "10x"
raw.filter$patient <- sample
rm(raw)

# normalization
data.norm <- NormalizeData(raw.filter, normalization.method = "LogNormalize", scale.factor = 10000)
# identification of highly variable features
data.norm <- FindVariableFeatures(data.norm, selection.method = "vst", nfeatures = 2000)
# Scaling
data.scale <- ScaleData(data.norm, features = VariableFeatures(data.norm))
# perform linear dimensional reduction
data <- RunPCA(data.scale, features = VariableFeatures(data.scale), npcs = 50)
# decide number of PCs
ElbowPlot(data,ndims = 50)
# non-linear dimensional reduction
data <- FindNeighbors(data, dims = 1:25)
data <- FindClusters(data, resolution = 0.8)
data <- RunUMAP(data, dims = 1:25, reduction = "pca")
rm(data.norm, data.scale)
# remove doublets
set.seed(123456)
doublet.rate <- 0.046 # defined by platforms of scRNA-seq libraries and total cells
sweep.res.list <- paramSweep_v3(data, PCs = 1:25, sct = FALSE, num.cores = 20)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
all.pk <- as.numeric(as.character(bcmvn$pK))
pk.best <- all.pk[bcmvn$BCmetric == max(bcmvn$BCmetric)]
nExp_poi <- round(doublet.rate*ncol(data))
data <- doubletFinder_v3(data, PCs = 1:25, pN = 0.25, pK = pk.best, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
colnames(data@meta.data)[ncol(data@meta.data)] <- "doublet"

# remove detected doublets, remove cells with genes > 5000, and save data for data integration
data <- subset(data, doublet == 'Singlet' & nFeature_RNA <= 5000)
data@meta.data <- data@meta.data[, 1:7]
saveRDS(data, "raw.filter.rds")
