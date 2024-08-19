library(Seurat)
library(tidyverse)
library(biomaRt)
library(GSVA)
library(mclust)
library(mltools)

trainset <- read_csv("C:\\Users\\atyrlik\\Downloads\\train_biomed_hackathon.csv") %>%
  .[, 1:4]
valset <- read_csv("C:\\Users\\atyrlik\\Downloads\\val_biomed_hackathon.csv") %>%
  .[, 1:4]


nrommat <- dat@assays$RNA$data
atacnormmat <- dat@assays$ATAC$data

nrommat <- nrommat[rowSums(nrommat) > 0,]
atacnormmat <- atacnormmat[rowSums(atacnormmat) > 0,]

traintest <- Hackathon2024.Testing.Set.Peak2Gene.Pairs

#PC embeddings
pcrna <- nrommat[rownames(nrommat) %in% c(VariableFeatures(dat), as.character(traintest$gene)), ] %>%
  .[,colSums(.) > 0] %>%
  prcomp(., scale. = TRUE, center = TRUE)

pcatac <- atacnormmat[rownames(atacnormmat) %in% c(VariableFeatures(dat, assay = "ATAC"), as.character(traintest$peak)), ] %>%
  .[,colSums(.) > 0] %>%
  prcomp(., scale. = TRUE, center = TRUE)

train <- as.data.frame(pcrna$x[,1:356]) %>%
  rename_with(.fn = ~ paste0("RNA_", .x)) %>% 
  rownames_to_column() %>%
  left_join(traintest, ., by = join_by("gene" == "rowname"))

train <- as.data.frame(pcatac$x[,1:671]) %>%
  rename_with(.fn = ~ paste0("ATAC_", .x)) %>% 
  rownames_to_column() %>%
  left_join(train, ., by = join_by("peak" == "rowname"))
rm(pcrna, pcatac)
#correlation

# zmat <-(nrommat-rowMeans(nrommat))/apply(nrommat, 1, sd)
ataczmat <-(atacnormmat-rowMeans(atacnormmat))/apply(atacnormmat, 1, sd)


atacsets <- split(ataczmat, rownames(ataczmat))
atacsets <- lapply(atacsets, function(x){ 
  names(x) <- colnames(ataczmat)
  x[x > 0]
})
names(atacsets) <- rownames(ataczmat)

# subgenesets <- genesets[names(genesets) %in% traintest$gene]
subatacsets <- atacsets[names(atacsets) %in% traintest$peak]

#gsva

namedsubatacsets <- lapply(subatacsets, names)

param <- gsvaParam(t(nrommat[rownames(nrommat) %in% traintest$gene,]), namedsubatacsets, minSize = 1)
gsvaout <- gsva(param, verbose=TRUE)

gsvascores <- apply(train, 1, function(x){
  if(x[["peak"]] %in% rownames(gsvaout) & x[["gene"]] %in% colnames(gsvaout)){gsvaout[x[["peak"]], x[["gene"]]]}
  else{NA}
})

train$gsva <- gsvascores
rm(ataczmat, param, gsvaout)
#physical distance
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
martout <- getBM(attributes=c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
      filters= 'hgnc_symbol',
      values= toupper(traintest$gene),
      mart=ensembl) %>%
  filter(chromosome_name %in% as.character(1:23)) %>%
  distinct(hgnc_symbol, .keep_all = TRUE)

out <- traintest %>%
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
          dist == 0 & between(start_position, pmin(start, end), pmax(start, end)) & between(end_position, pmin(start, end), pmax(start, end)) ~ abs(start_position-end_position),
          dist == 0 & between(pmin(start, end), pmin(start_position, end_position), pmax(start_position, end_position)) ~ pmax(start_position, end_position) -pmin(start, end),
          dist == 0 & between(pmax(start, end), pmin(start_position, end_position), pmax(start_position, end_position)) ~ pmax(start, end) - pmin(start_position, end_position),
          dist > 0 ~ 0,
          TRUE ~ NA
         ),
         reloverlap = if_else(dist == 0, overlap/abs(end_position-start_position), 0)
  )





train$dist <- out$dist
train$overlap <- out$overlap
train <- mutate(train, overlap = replace_na(overlap, 0))

cRNA <- AggregateExpression(dat, assays = "RNA", group.by = "RNA_snn_res.5", normalization.method = "RC")$RNA
colnames(cRNA) <- paste0("RNA_c", 1:ncol(cRNA))
cATAC <- AggregateExpression(dat, assays = "ATAC", group.by = "RNA_snn_res.5", normalization.method = "RC")$ATAC
colnames(cATAC) <- paste0("ATAC_c", 1:ncol(cATAC))

caRNA <- AggregateExpression(dat, assays = "RNA", group.by = "ATAC_snn_res.3", normalization.method = "RC")$RNA
colnames(caRNA) <- paste0("RNA_ca", 1:ncol(caRNA))
caATAC <- AggregateExpression(dat, assays = "ATAC", group.by = "ATAC_snn_res.3", normalization.method = "RC")$ATAC
colnames(caATAC) <- paste0("ATAC_ca", 1:ncol(caATAC))


train <- train %>%
  left_join(as.data.frame(cRNA) %>% rownames_to_column(), by = join_by("gene" == "rowname")) %>%
  left_join(as.data.frame(cATAC) %>% rownames_to_column(), by = join_by("peak" == "rowname")) %>%
  left_join(as.data.frame(caRNA) %>% rownames_to_column(), by = join_by("gene" == "rowname")) %>%
  left_join(as.data.frame(caATAC) %>% rownames_to_column(), by = join_by("peak" == "rowname"))
  
  
mattrain <- train %>%
  filter(Pair %in% trainset$Pair) %>%
  dplyr::select(matches("PC"))


# these <- train %>% 
#   filter(is.na(dist)) %>% 
#   magrittr::extract2("gene")

trainmodel <- MclustDA(as.matrix(mattrain), c(rep(TRUE, 20), rep(FALSE, 20)))
cout <- predict(trainmodel, newdata = as.matrix(train %>% filter(Pair %in% valset$Pair) %>% dplyr::select(matches("PC"))))

# check <- Hackathon2024.Training.Set.Peak2Gene.Pairs %>%
#   filter(!(gene %in% these))
check <- Hackathon2024.Training.Set.Peak2Gene.Pairs %>%
  filter(Pair %in% valset$Pair)
# check$classification <- out$classification
check$classification_clust <- cout$classification

# check$gsva <- train$gsva
# check$scoreoverlap <-out$z[,"TRUE"]
# check$scoreclust <-cout$z[,"TRUE"]

# check %>%
#   group_by(Peak2Gene, classification, classification_clust) %>%
#   summarise(n = n(), gsva = mean(gsva))


mcc(as.logical(check$Peak2Gene), as.logical(check$classification_clust))





ggplot(out, aes(x = reloverlap)) + geom_histogram() + facet_wrap(~Peak2Gene)

write_csv(train, "Hackathon2024.Testing.Set.Peak2Gene.Pairs.everything.csv")
# train <- read_csv("Hackathon2024.Training.Set.Peak2Gene.Pairs.everything.csv")

