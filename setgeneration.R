library(Seurat)
library(tidyverse)


nrommat <- dat@assays$RNA$data
atacnormmat <- dat@assays$ATAC$data

nrommat <- nrommat[rowSums(nrommat) > 0,]
atacnormmat <- atacnormmat[rowSums(atacnormmat) > 0,]


mean <- rowMeans(nrommat)
sd <- apply(nrommat, 1, sd)                
zmat <-(nrommat-mean)/sd

mean <- rowMeans(atacnormmat)
sd <- apply(atacnormmat, 1, sd)               
ataczmat <-(atacnormmat-mean)/sd

genesets <- split(zmat, rownames(zmat))
genesets <- lapply(genesets, function(x){ 
  names(x) <-colnames(zmat)
  x[x > 0]
})
names(genesets) <- rownames(nrommat)


atacsets <- split(ataczmat, rownames(ataczmat))
atacsets <- lapply(atacsets, function(x){ 
  names(x) <- colnames(ataczmat)
  x[x > 0]
})
names(atacsets) <- rownames(ataczmat)


# df <- as.data.frame(nrommat[rownames(nrommat) == "JAK1",])
# colnames(df) <- "values"
# ggplot(df, aes(x = values)) + geom_histogram()
