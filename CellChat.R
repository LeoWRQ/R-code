## Section 0 - Package preparation
# 0.1 Required packages for CellChat
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biobase")
BiocManager::install("BiocNeighbors")

source('http://renozao.github.io/repotools/install.R')
library(repotools)

install.packages('NMF', force = TRUE, quiet = TRUE)
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")

# 0.2 Install CellChat
suppressMessages(if(!require(CellChat))devtools::install_github("sqjin/CellChat"))

# 0.3 Install other supportive packages
suppressMessages(if(!require(Seurat))install.packages('Seurat'))
suppressMessages(reticulate::py_install(packages = 'umap-learn'))
suppressMessages(if(!require(ggplot2))install.packages('ggplot2'))
suppressMessages(if(!require(patchwork))install.packages('patchwork') )
suppressMessages(if(!require(ggalluvial))install.packages('ggalluvial'))
suppressMessages(if(!require(igraph))install.packages('igraph'))
suppressMessages(if(!require(dplyr))install.packages('dplyr'))
suppressMessages(options(stringsAsFactors = FALSE))
suppressMessages(options(futrue.globlas.Maxsize=2*1024**3))
suppressWarnings(suppressMessages(future::plan("multiprocess", workers = 8)))

# 0.4 Load packages
library(CellChat)
library(Seurat)
library(patchwork)

## Section 1 - Load data and database

# 1.1 Read Seurat data
SP.combined <- readRDS("~/ondemand/data/SP.combined.rds")

#cellchat <- createCellChat(object = SP.combined, group.by = "celltype") # Simply load data

# If select cell types
data.input <- SP.combined[["RNA"]]@data
meta <- SP.combined@meta.data
cell.use <- rownames(meta)[meta$celltype == c("CSF-cN", "Epend", "GABAN", "GluN", "Micro", "Astro")]
data.input <- data.input[, cell.use]
meta = meta[cell.use, ]
meta <- as.data.frame(meta)

# 1.2 Create CellChat
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
groupSize <- as.numeric(table(cellchat@idents))

# 1.2 Load database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Select database categories
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use

## Section 2 - Process the data
# 2.1 Pre-processing
cellchat <- subsetData(cellchat,features = NULL) 
cellchat <- identifyOverExpressedGenes(cellchat)

# 2.2 Identify over expressed pathways
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = T) #Default type = "truncatedMean"
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)
DT::datatable(df.net)
write.csv(df.net,'SP.df.net.csv')

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Section 3 - Visualization
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 targets.use = 'CSF-cN')

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 targets.use = 'Epend')





