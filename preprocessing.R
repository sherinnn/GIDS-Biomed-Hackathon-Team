library(Seurat)
library(tidyverse)

setwd("D:\\Users\\tyrli\\Documents\\Things\\R\\biosciencehackathon2024\\")

load("Hackathon2024.RData")


dat <- CreateSeuratObject(Hackathon2024.RNA, meta.data = Hackathon2024.Meta, assay = "RNA")
dat[["ATAC"]] <- CreateAssay5Object(counts = Hackathon2024.ATAC)

dat <- PercentageFeatureSet(dat, pattern = "^MT-", col.name = "mitofrac") #makes sure to start at the start of the string
VlnPlot(dat, features = c("nCount_RNA", "nFeature_RNA", "mitofrac"))
dat <- subset(dat, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & mitofrac < 20)

dat <- dat %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(dat)

dat <- dat %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.3) %>%
  RunUMAP(dims = 1:20)

dat <- dat %>%
  NormalizeData(assay = "ATAC") %>%
  FindVariableFeatures(assay = "ATAC") %>%
  ScaleData(assay = "ATAC") %>%
  RunPCA(assay = "ATAC", reduction.name = "pca_ATAC")

ElbowPlot(dat, reduction = "pca_ATAC")

dat <- dat %>%
  FindNeighbors(dims = 1:14, reduction = "pca_ATAC") %>%
  FindClusters(resolution = 0.3) %>%
  RunUMAP(dims = 1:14, reduction = "pca_ATAC", reduction.name = "umap_ATAC")

