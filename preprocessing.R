library(Seurat)
library(tidyverse)

# setwd("C:\\Users\\atyrlik\\Documents\\hackathon\\")
setwd("D:\\Users\\tyrli\\Documents\\Things\\R\\biosciencehackathon2024")
load("Hackathon2024.RData")


dat <- CreateSeuratObject(Hackathon2024.RNA, meta.data = Hackathon2024.Meta, assay = "RNA")
dat[["ATAC"]] <- CreateAssay5Object(counts = Hackathon2024.ATAC)

dat <- PercentageFeatureSet(dat, pattern = "^MT-", col.name = "mitofrac") #makes sure to start at the start of the string
VlnPlot(dat, features = c("nCount_ATAC", "nFeature_ATAC", "mitofrac"), group.by = "orig.ident")
dat <- subset(dat, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & mitofrac < 20 & nFeature_ATAC > 500)

dat <- dat %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(dat)

dat <- dat %>%
  FindNeighbors(dims = 1:50) %>%
  FindClusters(resolution = 5, algorithm = 2) %>%
  RunUMAP(dims = 1:50)

dat <- dat %>%
  NormalizeData(assay = "ATAC") %>%
  FindVariableFeatures(assay = "ATAC")# %>%
#   ScaleData(assay = "ATAC") %>%
#   RunPCA(assay = "ATAC")
# 
# ElbowPlot(dat, reduction = "pca_ATAC")
# 
# dat <- dat %>%
#   FindNeighbors(dims = 1:50, reduction = "pca_ATAC") %>%
#   FindClusters(., resolution = 0.2, graph.name = "ATAC_snn") %>%
#   RunUMAP(dims = 1:50, reduction = "pca_ATAC", reduction.name = "umap_ATAC")

cRNA <- AggregateExpression(dat, assays = "RNA", group.by = "RNA_snn_res.4", normalization.method = "RC")$RNA
colnames(cRNA) <- paste0("RNA_c", 1:ncol(cRNA))
cATAC <- AggregateExpression(dat, assays = "ATAC", group.by = "RNA_snn_res.4", normalization.method = "RC")$ATAC
colnames(cATAC) <- paste0("ATAC_c", 1:ncol(cATAC))
