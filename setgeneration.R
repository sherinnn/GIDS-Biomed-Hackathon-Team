library(Seurat)
library(tidyverse)
library(biomaRt)

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

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
martout <- getBM(attributes=c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
      filters= 'hgnc_symbol',
      values= toupper(Hackathon2024.Training.Set.Peak2Gene.Pairs$gene),
      mart=ensembl) %>%
  filter((chromosome_name %in% as.character(1:23)))

train <- Hackathon2024.Training.Set.Peak2Gene.Pairs %>%
  separate_wider_delim(col = peak, names = c("ch", "start", "end"), cols_remove = FALSE, delim = "-") %>%
  left_join(martout, by = join_by("gene" == "hgnc_symbol")) %>%
  mutate(chromosome_name = paste0("chr", chromosome_name),
         start = as.numeric(start),
         end = as.numeric(end),
         start_position = as.numeric(start_position),
         end_position = as.numeric(end_position),
         dist = case_when(
         ch == chromosome_name & pmax(start, end) < pmin(start_position, end_position) ~ pmin(start_position, end_position)-pmax(start, end),
         ch == chromosome_name & pmin(start, end) > pmax(start_position, end_position) ~ pmin(start, end)-pmax(start_position, end_position),
         ch == chromosome_name ~ 0,
         TRUE ~ NA
         )
  )

ggplot(train, aes(x = dist)) + geom_histogram() + facet_wrap(~Peak2Gene)


