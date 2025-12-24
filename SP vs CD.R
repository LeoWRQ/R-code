library(Seurat)
library(dplyr)
library(mindr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggrepel)

## Section 1 - Load all data sets, perform basic quality control, and rename cell id

## 1.1 - Read all data sets

# From CD paper
CD1.data <- Read10X(data.dir = "~/Desktop/GSE197547_RAW/LPS")
CD1 <- CreateSeuratObject(counts = SP1.data, project = "CD_LPS")
CD2.data <- Read10X(data.dir = "~/Desktop/GSE197547_RAW/Saline")
CD2 <- CreateSeuratObject(counts = SP2.data, project = "CD_Ctrl")

# Our LPS vs PBS data
SP3.data <- Read10X(data.dir = "~/Desktop/SP3")
SP3 <- CreateSeuratObject(counts = SP3.data, project = "LPS1")
SP4.data <- Read10X(data.dir = "~/Desktop/SP4")
SP4 <- CreateSeuratObject(counts = SP4.data, project = "Ctrl3")

## 1.2 - Quality control
mito.features <- grep(pattern = "^mt-", x = rownames(x = CD1), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = CD1, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = CD1, slot = 'counts'))
CD1[['percent.mito']] <- percent.mito
CD1.QC <- subset(x = CD1, subset = nFeature_RNA > 1000 & nFeature_RNA < 20000 & percent.mito < 0.10)

mito.features <- grep(pattern = "^mt-", x = rownames(x = CD2), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = CD2, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = CD2, slot = 'counts'))
CD2[['percent.mito']] <- percent.mito
CD2.QC <- subset(x = CD2, subset = nFeature_RNA > 1000 & nFeature_RNA < 20000 & percent.mito < 0.10)

mito.features <- grep(pattern = "^mt-", x = rownames(x = SP3), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = SP3, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = SP3, slot = 'counts'))
SP3[['percent.mito']] <- percent.mito
SP3.QC <- subset(x = SP3, subset = nFeature_RNA > 1000 & nFeature_RNA < 20000 & percent.mito < 0.10)

mito.features <- grep(pattern = "^mt-", x = rownames(x = SP4), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = SP4, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = SP4, slot = 'counts'))
SP4[['percent.mito']] <- percent.mito
SP4.QC <- subset(x = SP4, subset = nFeature_RNA > 1000 & nFeature_RNA < 20000 & percent.mito < 0.10)

## 1.3 - Rename the cell id
CD1.QC <- RenameCells(CD1.QC, add.cell.id = "CD1")
CD2.QC <- RenameCells(CD2.QC, add.cell.id = "CD2")
SP3.QC <- RenameCells(SP3.QC, add.cell.id = "SP3")
SP4.QC <- RenameCells(SP4.QC, add.cell.id = "SP4")


CD1.QC$Stim = CD1.QC$orig.ident
head(CD1.QC[[]])
levels(CD1.QC$Stim) <- c(levels(CD1.QC$Stim), "LPS")
CD1.QC$Stim[cells = WhichCells(CD1.QC)] <- "LPS"
head(CD1.QC[[]])

CD2.QC$Stim = CD2.QC$orig.ident
levels(CD2.QC$Stim) <- c(levels(CD2.QC$Stim), "Ctrl")
CD2.QC$Stim[cells = WhichCells(CD2.QC)] <- "Ctrl"
head(CD2.QC[[]])

SP3.QC$Stim = SP3.QC$orig.ident
levels(SP3.QC$Stim) <- c(levels(SP3.QC$Stim), "LPS")
SP3.QC$Stim[cells = WhichCells(SP3.QC)] <- "LPS"
head(SP3.QC[[]])

SP4.QC$Stim = SP4.QC$orig.ident
levels(SP4.QC$Stim) <- c(levels(SP4.QC$Stim), "Ctrl")
SP4.QC$Stim[cells = WhichCells(SP4.QC)] <- "Ctrl"
head(SP4.QC[[]])

## 1.4 - Normalization
SPCD.QC <- list(CD1.QC, CD2.QC, SP3.QC, SP4.QC)
SPCD.QC <- lapply(X = SPCD.QC, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = SPCD.QC)

## 1.5 - Combine all data sets
SPCD.anchors <- FindIntegrationAnchors(object.list = SPCD.QC, anchor.features = features)
SPCD.combined <- IntegrateData(anchorset = SPCD.anchors)

## Section 2 - Perform an integrated analysis
# 2.1 - Scaling and dimension deduction
DefaultAssay(SPCD.combined) <- "integrated"
SPCD.combined <- ScaleData(SPCD.combined, verbose = FALSE)
SPCD.combined <- RunPCA(object = SPCD.combined, verbose = FALSE, dims = 50)
SPCD.combined <- RunUMAP(object = SPCD.combined, dims = 1:50)
SPCD.combined <- FindNeighbors(object = SPCD.combined, dims = 1:50)
SPCD.combined <- FindClusters(object = SPCD.combined, resolution = 0.1)
FeaturePlot(SPCD.combined, features = c("nFeature_RNA", "nCount_RNA"))

# 2.2 - Dimplot for samples and conditions
DefaultAssay(SPCD.combined) <- 'RNA'
Idents(SPCD.combined) <- "seurat_clusters"
DimPlot(object = SPCD.combined, reduction = 'umap', label = TRUE)
DimPlot(object = SPCD.combined, reduction = 'umap', label = FALSE, group.by ='Stim')
table(SPCD.combined$seurat_clusters)

Idents(SPCD.combined) <- "orig.ident"
DimPlot(object = SPCD.combined, reduction = 'umap', label = FALSE, group.by = "orig.ident")


## Section 3 - Explore cell clusters

# 3.0 - Visualization of cell clusters
DimPlot(object = SPCD.combined, reduction = 'umap', label = TRUE)
# Ependymal Cells 
FeaturePlot(SPCD.combined, features = c("Cfap43", "Dnah12", "Dcdc2a"))
FeaturePlot(SPCD.combined, features = c("Dcdc2a", "Tppp3", "S100b", "Foxj1", "Vim"))
FeaturePlot(SPCD.combined, features = c("Pifo", "Ccdc153", "Rabl2",  "Vim"))

Idents(SPCD.combined) <- "seurat_clusters"
SPCD.combined <- RenameIdents(SPCD.combined, `0` = "Epend")

## Section 4 - Ependymal cells

Idents(SPCD.combined) <- "seurat_clusters"
SPCD.Epend <- subset(SPCD.combined, idents = "Epend")
DefaultAssay(SPCD.Epend) <- "integrated"
SPCD.Epend <- ScaleData(SPCD.Epend, verbose = FALSE)
SPCD.Epend <- RunPCA(object = SPCD.Epend, verbose = FALSE, dims = 50)
SPCD.Epend <- RunUMAP(object = SPCD.Epend, dims = 1:50)
SPCD.Epend <- FindNeighbors(object = SPCD.Epend, dims = 1:50)
SPCD.Epend <- FindClusters(object = SPCD.Epend, resolution = 0.05)
table(SPCD.Epend$seurat_clusters)

Idents(SPCD.Epend) <- "orig.ident"
DimPlot(object = SPCD.Epend, reduction = 'umap', label = FALSE, group.by = "orig.ident")

DimPlot(object = SPCD.Epend, reduction = 'umap', label = FALSE, group.by ='Stim')

DefaultAssay(SPCD.Epend) <- 'RNA'
Idents(SPCD.Epend) <- "Stim"
VlnPlot(SPCD.Epend, features = c("Fos"), pt.size=0)

Idents(SPCD.Epend) <- "orig.ident"
VlnPlot(SPCD.Epend, features = c("Fos"), pt.size=0, split.by = "orig.ident")
FeaturePlot(SPCD.Epend, features = c("Fos"), split.by = "orig.ident")

## Only Control

SPCD.Ctrl.QC <- list(CD2.QC, SP4.QC)
SPCD.Ctrl.QC <- lapply(X = SPCD.Ctrl.QC, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = SPCD.Ctrl.QC)

SPCD.Ctrl.anchors <- FindIntegrationAnchors(object.list = SPCD.Ctrl.QC, anchor.features = features)
SPCD.Ctrl <- IntegrateData(anchorset = SPCD.Ctrl.anchors)

DefaultAssay(SPCD.Ctrl) <- "integrated"
SPCD.Ctrl <- ScaleData(SPCD.Ctrl, verbose = FALSE)
SPCD.Ctrl <- RunPCA(object = SPCD.Ctrl, verbose = FALSE, dims = 50)
SPCD.Ctrl <- RunUMAP(object = SPCD.Ctrl, dims = 1:50)
SPCD.Ctrl <- FindNeighbors(object = SPCD.Ctrl, dims = 1:50)
SPCD.Ctrl <- FindClusters(object = SPCD.Ctrl, resolution = 0.1)
FeaturePlot(SPCD.Ctrl, features = c("nFeature_RNA", "nCount_RNA"))

DimPlot(object = SPCD.Ctrl, reduction = 'umap', label = TRUE)
Idents(SPCD.Ctrl) <- "orig.ident"
DimPlot(object = SPCD.Ctrl, reduction = 'umap', label = FALSE, group.by = "orig.ident")

# Ependymal Cells 
DefaultAssay(SPCD.Ctrl) <- 'RNA'
FeaturePlot(SPCD.Ctrl, features = c("Cfap43", "Dnah12", "Dcdc2a"))
FeaturePlot(SPCD.Ctrl, features = c("Dcdc2a", "Tppp3", "S100b", "Foxj1", "Vim"))
FeaturePlot(SPCD.Ctrl, features = c("Pifo", "Ccdc153", "Rabl2",  "Vim"))

Idents(SPCD.Ctrl) <- "seurat_clusters"
SPCD.Ctrl <- RenameIdents(SPCD.Ctrl, `0` = "Epend")
SPCD.Epend.Ctrl <- subset(SPCD.Ctrl, idents = "Epend")
DefaultAssay(SPCD.Epend.Ctrl) <- "integrated"
SPCD.Epend.Ctrl <- ScaleData(SPCD.Epend.Ctrl, verbose = FALSE)
SPCD.Epend.Ctrl <- RunPCA(object = SPCD.Epend.Ctrl, verbose = FALSE, dims = 50)
SPCD.Epend.Ctrl <- RunUMAP(object = SPCD.Epend.Ctrl, dims = 1:50)
SPCD.Epend.Ctrl <- FindNeighbors(object = SPCD.Epend.Ctrl, dims = 1:50)
SPCD.Epend.Ctrl <- FindClusters(object = SPCD.Epend.Ctrl, resolution = 0.05)
table(SPCD.Epend.Ctrl$seurat_clusters)

Idents(SPCD.Epend.Ctrl) <- "orig.ident"
DefaultAssay(SPCD.Epend.Ctrl) <- 'RNA'
DimPlot(object = SPCD.Epend.Ctrl, reduction = 'umap', label = FALSE, group.by = "orig.ident")

Epend.brainspine <- FindAllMarkers(SPCD.Epend.Ctrl, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25)
write.table(Epend.brainspine, "~/Desktop/Epend.brainspine.txt", sep="\t")

FeaturePlot(SPCD.Epend.Ctrl, features = c("Acta2"), split.by = "orig.ident")

Epend.brainspine.2 <- FindMarkers(SPCD.Epend.Ctrl, ident.1 = "CD_Ctrl", ident.2 = "Ctrl3", verbose = FALSE)
write.table(Epend.brainspine.2, "~/Desktop/Epend.brainspine.2.txt", sep="\t")
Epend.brainspine.3 <- FindMarkers(SPCD.Epend.Ctrl, ident.1 = "Ctrl3", ident.2 = "CD_Ctrl", verbose = FALSE)









