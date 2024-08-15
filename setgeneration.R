library(Seurat)
library(tidyverse)
library(biomaRt)

nrommat <- dat@assays$RNA$data
atacnormmat <- dat@assays$ATAC$data

nrommat <- nrommat[rowSums(nrommat) > 0,]
atacnormmat <- atacnormmat[rowSums(atacnormmat) > 0,]



pcrna <- prcomp(dat@assays$RNA$scale.data)
pcatac <- prcomp(dat@assays$ATAC$scale.data)


pcatac$x[,1:100]

train <- as.data.frame(pcrna$x[,1:100]) %>%
  rename_with(.fn = ~ paste0("RNA_", .x)) %>% 
  rownames_to_column() %>%
  inner_join(Hackathon2024.Training.Set.Peak2Gene.Pairs, ., by = join_by("gene" == "rowname"))

train <- as.data.frame(pcatac$x[,1:1000]) %>%
  rename_with(.fn = ~ paste0("ATAC_", .x)) %>% 
  rownames_to_column() %>%
  inner_join(train, ., by = join_by("peak" == "rowname"))



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

subgenesets <- genesets[names(genesets) %in% Hackathon2024.Training.Set.Peak2Gene.Pairs$gene]
subatacsets <- atacsets[names(atacsets) %in% Hackathon2024.Training.Set.Peak2Gene.Pairs$peak]


temp <- lapply(subgenesets, function(x){
  out <- lapply(subatacsets, function(y){
    rnadf <- as.data.frame(x) %>%
      rownames_to_column()
    colnames(rnadf) <- c("rname","RNA")
    atacdf <- as.data.frame(y) %>%
      rownames_to_column()
    colnames(atacdf) <- c("rname","ATAC")
    df <- full_join(rnadf, atacdf, by = join_by("rname"))%>%
      mutate(RNA = replace_na(RNA, 0),
             ATAC = replace_na(ATAC, 0))
    summary(lm(RNA~ATAC, data=df))$r.squared
  })
})
summary(lm(RNA~ATAC, data=df))$r.squared
tempmat <- do.call(rbind, temp)

scores <- apply(train, 1, function(x){
  tempmat[as.character(x[["gene"]]),x[["peak"]]]
})
train$corr <- unlist(scores)

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
         ), 
         overlap = case_when(
          dist == 0 & between(start, pmin(start_position, end_position), pmax(start_position, end_position)) & between(start, pmin(start_position, end_position), pmax(start_position, end_position)) ~ abs(start-end),
          dist == 0 & between(pmin(start, end), pmin(start_position, end_position), pmax(start_position, end_position)) ~ pmax(start_position, end_position) -pmin(start, end),
          dist == 0 & between(pmax(start, end), pmin(start_position, end_position), pmax(start_position, end_position)) ~ pmax(start, end) - pmin(start_position, end_position),
          dist > 0 ~ 0,
          TRUE ~ NA
         )
  )

prcomp(nrommat, scale = TRUE)


ggplot(train, aes(x = corr)) + geom_histogram() + facet_wrap(~Peak2Gene)

write_csv(train, "D:\\Users\\tyrli\\Documents\\Hackathon2024.Training.Set.Peak2Gene.Pairs.pcloadings.csv")

