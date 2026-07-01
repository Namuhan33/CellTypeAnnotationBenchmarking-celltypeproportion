setwd("R:/ANNCOM2023-Q6025")
library(Seurat)
library(dplyr)

pbmc.68k <- readRDS(file = "08_StartOver_NewStory/05_interdatasets/01_pbmc/01_datasets/seuobj/pbmc_68k.rds")
pbmc.covid <- readRDS(file = "08_StartOver_NewStory/05_interdatasets/01_pbmc/01_datasets/seuobj/pbmc_covid.rds")

pbmc.68k@meta.data$Index <- 1:nrow(pbmc.68k@meta.data)
pbmc.covid@meta.data$Index <- 1:nrow(pbmc.covid@meta.data)

table(pbmc.covid$cell_type_1)/sum(table(pbmc.covid$cell_type_1))*100
table(pbmc.covid$cell_type_2)/sum(table(pbmc.covid$cell_type_2))*100
table(pbmc.covid$cell_type_3)/sum(table(pbmc.covid$cell_type_3))*100

table(pbmc.68k$cell_type_1)/sum(table(pbmc.68k$cell_type_1))*100
table(pbmc.68k$cell_type_2)/sum(table(pbmc.68k$cell_type_2))*100
table(pbmc.68k$cell_type_3)/sum(table(pbmc.68k$cell_type_3))*100

# cell count each iteration (depends on myeloid in reference) is 5238*(1+0.1)=5761 
# Session 1: query counts = reference counts
covid.meta <- pbmc.covid@meta.data
health.meta <- pbmc.68k@meta.data

query.index.s1 <- covid.meta %>% group_by(cell_type_1) %>% slice_sample(prop=(5761+5)*5/nrow(covid.meta))
index <- rBayesianOptimization::KFold(target = query.index.s1$cell_type_1, nfolds = 5, stratified = TRUE, seed = 233)
query.index.s1 <- lapply(index, function(ind){query.index.s1[ind,]$Index})

# P1 = 10:1 (lymphoid:myeloid)
reference.index.p1 <- lapply(1:5, function(i){
  p=10
  lymphoid <- health.meta[health.meta$cell_type_3=="Lymphoid",]
  myeloid <- health.meta[health.meta$cell_type_3=="Myeloid",]
  lym <- lymphoid %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(5761+5)/nrow(lymphoid))
  mye <- myeloid %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(5761+5)/nrow(myeloid))
  df <- rbind(lym, mye)
  df$Index
})

table(health.meta[reference.index.p1[[1]],]$cell_type_3)

index.p1 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], reference.index.p1[[i]])
})

# P2 = 5:1 (lymphoid:myeloid)
reference.index.p2 <- lapply(1:5, function(i){
  p=5
  lymphoid <- health.meta[health.meta$cell_type_3=="Lymphoid",]
  myeloid <- health.meta[health.meta$cell_type_3=="Myeloid",]
  lym <- lymphoid %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(5761+5)/nrow(lymphoid))
  mye <- myeloid %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(5761+5)/nrow(myeloid))
  df <- rbind(lym, mye)
  df$Index
})

table(health.meta[reference.index.p2[[1]],]$cell_type_3)

index.p2 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], reference.index.p2[[i]])
})

# P3 = 2:1 (lymphoid:myeloid)
reference.index.p3 <- lapply(1:5, function(i){
  p=2
  lymphoid <- health.meta[health.meta$cell_type_3=="Lymphoid",]
  myeloid <- health.meta[health.meta$cell_type_3=="Myeloid",]
  lym <- lymphoid %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(5761+5)/nrow(lymphoid))
  mye <- myeloid %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(5761+5)/nrow(myeloid))
  df <- rbind(lym, mye)
  df$Index
})

table(health.meta[reference.index.p3[[1]],]$cell_type_3)

index.p3 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], reference.index.p3[[i]])
})

# P4 = 1:1 (lymphoid:myeloid)
reference.index.p4 <- lapply(1:5, function(i){
  p=1
  lymphoid <- health.meta[health.meta$cell_type_3=="Lymphoid",]
  myeloid <- health.meta[health.meta$cell_type_3=="Myeloid",]
  lym <- lymphoid %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(5761+5)/nrow(lymphoid))
  mye <- myeloid %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(5761+5)/nrow(myeloid))
  df <- rbind(lym, mye)
  df$Index
})

table(health.meta[reference.index.p4[[1]],]$cell_type_3)

index.p4 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], reference.index.p4[[i]])
})

# P5 = 0.5:1 (lymphoid:myeloid)
reference.index.p5 <- lapply(1:5, function(i){
  p=0.5
  lymphoid <- health.meta[health.meta$cell_type_3=="Lymphoid",]
  myeloid <- health.meta[health.meta$cell_type_3=="Myeloid",]
  lym <- lymphoid %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(5761+5)/nrow(lymphoid))
  mye <- myeloid %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(5761+5)/nrow(myeloid))
  df <- rbind(lym, mye)
  df$Index
})

table(health.meta[reference.index.p5[[1]],]$cell_type_3)

index.p5 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], reference.index.p5[[i]])
})

# P6 = 0.2:1 (lymphoid:myeloid)
reference.index.p6 <- lapply(1:5, function(i){
  p=0.2
  lymphoid <- health.meta[health.meta$cell_type_3=="Lymphoid",]
  myeloid <- health.meta[health.meta$cell_type_3=="Myeloid",]
  lym <- lymphoid %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(5761+5)/nrow(lymphoid))
  mye <- myeloid %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(5761+5)/nrow(myeloid))
  df <- rbind(lym, mye)
  df$Index
})

table(health.meta[reference.index.p6[[1]],]$cell_type_3)

index.p6 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], reference.index.p6[[i]])
})

# P7 = 0.1:1 (lymphoid:myeloid)
reference.index.p7 <- lapply(1:5, function(i){
  p=0.1
  lymphoid <- health.meta[health.meta$cell_type_3=="Lymphoid",]
  myeloid <- health.meta[health.meta$cell_type_3=="Myeloid",]
  lym <- lymphoid %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(5761+5)/nrow(lymphoid))
  mye <- myeloid %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(5761+5)/nrow(myeloid))
  df <- rbind(lym, mye)
  df$Index
})

table(health.meta[reference.index.p7[[1]],]$cell_type_3)

index.p7 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], reference.index.p7[[i]])
})

# Session 2: all cells in query will be used
covid.meta <- pbmc.covid@meta.data
health.meta <- pbmc.68k@meta.data

index <- rBayesianOptimization::KFold(target = covid.meta$cell_type_1, nfolds = 5, stratified = TRUE, seed = 233)
query.index.s2 <- lapply(index, function(ind){covid.meta[ind,]$Index})
# p8 = 10:1 (lymphoid:myeloid), all query cells are invited
index.p8 <- lapply(1:5, function(i){
  n=8
  ref <- get(paste0("reference.index.p", n-7))
  list(query.index.s2[[i]], ref[[i]])
})

# p9 = 5:1 (lymphoid:myeloid), all query cells are invited
index.p9 <- lapply(1:5, function(i){
  n=9
  ref <- get(paste0("reference.index.p", n-7))
  list(query.index.s2[[i]], ref[[i]])
})

# p10 = 2:1 (lymphoid:myeloid), all query cells are invited
index.p10 <- lapply(1:5, function(i){
  n=10
  ref <- get(paste0("reference.index.p", n-7))
  list(query.index.s2[[i]], ref[[i]])
})

# p11 = 1:1 (lymphoid:myeloid), all query cells are invited
index.p11 <- lapply(1:5, function(i){
  n=11
  ref <- get(paste0("reference.index.p", n-7))
  list(query.index.s2[[i]], ref[[i]])
})

# p12 = 0.5:1 (lymphoid:myeloid), all query cells are invited
index.p12 <- lapply(1:5, function(i){
  n=12
  ref <- get(paste0("reference.index.p", n-7))
  list(query.index.s2[[i]], ref[[i]])
})

# p13 = 0.2:1 (lymphoid:myeloid), all query cells are invited
index.p13 <- lapply(1:5, function(i){
  n=13
  ref <- get(paste0("reference.index.p", n-7))
  list(query.index.s2[[i]], ref[[i]])
})

# p14 = 0.1:1 (lymphoid:myeloid), all query cells are invited
index.p14 <- lapply(1:5, function(i){
  n=14
  ref <- get(paste0("reference.index.p", n-7))
  list(query.index.s2[[i]], ref[[i]])
})

for (i in 1:14) {
  saveRDS(get(paste0("index.p", i)), file = paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/01_datasets/index/Index_P", i, ".rds"))
}


lengths(Index_P7[[1]])
lengths(Index_P14[[1]])


