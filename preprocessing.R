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
  FindNeighbors(dims = 1:14) %>%
  FindClusters(resolution = 0.3) %>%
  RunUMAP(dims = 1:14)

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

norm <- VST(Hackathon2024.RNA)

vars <- apply(Hackathon2024.RNA, 1, var)
mean <- apply(Hackathon2024.RNA, 1, mean)
norm <- VST(Hackathon2024.RNA)

df <- data.frame(vars = vars, mean = mean, stvars = norm$variance.standardized)

rna <- Hackathon2024.RNA[rownames(Hackathon2024.RNA) == "YPEL5",]
atac <- Hackathon2024.ATAC[rownames(Hackathon2024.ATAC) == "chr2-30144160-30169627",]

df <- data.frame(rna = rna, atac = atac)
ggplot(df, aes(x = rna, y = atac)) + geom_point() +
  geom_smooth(method = "lm")

train <- Hackathon2024.Training.Set.Peak2Gene.Pairs

apply(Hackathon2024.RNA, 1, function(x){
  this <- Hackathon2024.Training.Set.Peak2Gene.Pairs[Hackathon2024.Training.Set.Peak2Gene.Pairs$gene == rownames(x),]$peak
  y <- Hackathon2024.ATAC[rownames(Hackathon2024.ATAC) == this,]
  # cor(x, y)
})
