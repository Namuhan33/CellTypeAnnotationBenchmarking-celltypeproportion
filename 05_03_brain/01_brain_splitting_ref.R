setwd("R:/ANNCOM2023-Q6025")

library(Seurat)
library(SeuratObject)
library(Matrix)
library(tidyverse)

seaAD <- readRDS("08_StartOver_NewStory/05_interdatasets/03_brain/01_datasets/seuobj/seaAD.rds")
humanMTG <- readRDS("08_StartOver_NewStory/05_interdatasets/03_brain/01_datasets/seuobj/humanMTG.rds")

# splitting datasets ####
# seaAD always be the query and humanMTG always be the reference
seaAD@meta.data$index <- 1:ncol(seaAD)
humanMTG@meta.data$index <- 1:ncol(humanMTG)

## Session 1: counts in query = counts in reference #####
# 10:1, 5:1, 2:1, 1:1, 0.5:1, 0.2:1, 0.1:1. Total cell count is 1411*1.1=1552
seaAD.meta <- seaAD@meta.data
humanMTG.meta <- humanMTG@meta.data

seaAD.meta <- seaAD.meta %>% group_by(cell_type_1) %>% slice_sample(prop = (1552*5+20)/nrow(seaAD.meta))
index.list.Q <- rBayesianOptimization::KFold(target = seaAD.meta$cell_type_1, nfolds = 5, stratified = TRUE, seed = 233)
index.list.Q <- lapply(index.list.Q, function(ind){
  seaAD.meta[ind,]$index
})

ref.neu <- humanMTG.meta[humanMTG.meta$cell_type_4=="Neuronal",]
ref.non <- humanMTG.meta[humanMTG.meta$cell_type_4=="NonNeuronal",]

# P1 = 10:1, count Q= count R
index.list.R.p1 <- lapply(1:5, function(i){
  p=10
  ref.index.neu <- ref.neu %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(1552+5)/nrow(ref.neu))
  ref.index.non <- ref.non %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(1552+5)/nrow(ref.non))
  ref <- rbind(ref.index.neu, ref.index.non)
  ref$index
})

index.list.p1 <- lapply(1:5, function(i){
  list(index.list.Q[[i]], index.list.R.p1[[i]])
})

# P2 = 5:1, count Q = count R
index.list.R.p2 <- lapply(1:5, function(i){
  p=5
  ref.index.neu <- ref.neu %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(1552+5)/nrow(ref.neu))
  ref.index.non <- ref.non %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(1552+5)/nrow(ref.non))
  ref <- rbind(ref.index.neu, ref.index.non)
  ref$index
})

index.list.p2 <- lapply(1:5, function(i){
  list(index.list.Q[[i]], index.list.R.p2[[i]])
})

table(humanMTG.meta[index.list.p2[[1]][[2]],]$cell_type_4)

# P3 = 2:1, count Q = count R
index.list.R.p3 <- lapply(1:5, function(i){
  p=2
  ref.index.neu <- ref.neu %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(1552+5)/nrow(ref.neu))
  ref.index.non <- ref.non %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(1552+5)/nrow(ref.non))
  ref <- rbind(ref.index.neu, ref.index.non)
  ref$index
})

index.list.p3 <- lapply(1:5, function(i){
  list(index.list.Q[[i]], index.list.R.p3[[i]])
})

table(humanMTG.meta[index.list.p3[[1]][[2]],]$cell_type_4)

# P4 = 1:1, count Q = count R
index.list.R.p4 <- lapply(1:5, function(i){
  p=1
  ref.index.neu <- ref.neu %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(1552+5)/nrow(ref.neu))
  ref.index.non <- ref.non %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(1552+5)/nrow(ref.non))
  ref <- rbind(ref.index.neu, ref.index.non)
  ref$index
})

index.list.p4 <- lapply(1:5, function(i){
  list(index.list.Q[[i]], index.list.R.p4[[i]])
})

table(humanMTG.meta[index.list.p4[[1]][[2]],]$cell_type_4)

# P5 = 1:1, count Q = count R
index.list.R.p5 <- lapply(1:5, function(i){
  p=0.5
  ref.index.neu <- ref.neu %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(1552+5)/nrow(ref.neu))
  ref.index.non <- ref.non %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(1552+5)/nrow(ref.non))
  ref <- rbind(ref.index.neu, ref.index.non)
  ref$index
})

index.list.p5 <- lapply(1:5, function(i){
  list(index.list.Q[[i]], index.list.R.p5[[i]])
})

table(humanMTG.meta[index.list.p5[[1]][[2]],]$cell_type_4)

# P6 = 0.2:1, count Q = count R
index.list.R.p6 <- lapply(1:5, function(i){
  p=0.2
  ref.index.neu <- ref.neu %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(1552+5)/nrow(ref.neu))
  ref.index.non <- ref.non %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(1552+5)/nrow(ref.non))
  ref <- rbind(ref.index.neu, ref.index.non)
  ref$index
})

index.list.p6 <- lapply(1:5, function(i){
  list(index.list.Q[[i]], index.list.R.p6[[i]])
})

table(humanMTG.meta[index.list.p6[[1]][[2]],]$cell_type_4)

# P7 = 0.1:1, count Q = count R
index.list.R.p7 <- lapply(1:5, function(i){
  p=0.1
  ref.index.neu <- ref.neu %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(1552+5)/nrow(ref.neu))
  ref.index.non <- ref.non %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(1552+5)/nrow(ref.non))
  ref <- rbind(ref.index.neu, ref.index.non)
  ref$index
})

index.list.p7 <- lapply(1:5, function(i){
  list(index.list.Q[[i]], index.list.R.p7[[i]])
})

table(humanMTG.meta[index.list.p7[[1]][[2]],]$cell_type_4)

## Session 2: counts in query are the originals #####
# 10:1, 5:1, 2:1, 1:1, 0.5:1, 0.2:1, 0.1:1. Total cell count in reference is 1411*1.1=1552
seaAD.meta <- seaAD@meta.data
humanMTG.meta <- humanMTG@meta.data

index.list.Q.S2 <- rBayesianOptimization::KFold(target = seaAD.meta$cell_type_1, nfolds = 5, stratified = TRUE, seed = 233)
index.list.Q.S2 <- lapply(index.list.Q.S2, function(ind){
  seaAD.meta[ind,]$index
})

# P8 = 10:1, reference = 1552, query = 137303/5 = 27460
index.list.p8 <- lapply(1:5, function(i){
  pn=8
  list.index <- get(paste0("index.list.R.p", pn-7))
  list(index.list.Q.S2[[i]], list.index[[i]])
})

# P9 = 5:1, reference = 1552, query = 137303/5 = 27460
index.list.p9 <- lapply(1:5, function(i){
  pn=9
  list.index <- get(paste0("index.list.R.p", pn-7))
  list(index.list.Q.S2[[i]], list.index[[i]])
})

# P10 = 2:1, reference = 1552, query = 137303/5 = 27460
index.list.p10 <- lapply(1:5, function(i){
  pn=10
  list.index <- get(paste0("index.list.R.p", pn-7))
  list(index.list.Q.S2[[i]], list.index[[i]])
})

# P11 = 1:1, reference = 1552, query = 137303/5 = 27460
index.list.p11 <- lapply(1:5, function(i){
  pn=11
  list.index <- get(paste0("index.list.R.p", pn-7))
  list(index.list.Q.S2[[i]], list.index[[i]])
})

# P12 = 0.5:1, reference = 1552, query = 137303/5 = 27460
index.list.p12 <- lapply(1:5, function(i){
  pn=12
  list.index <- get(paste0("index.list.R.p", pn-7))
  list(index.list.Q.S2[[i]], list.index[[i]])
})

# P13 = 0.2:1, reference = 1552, query = 137303/5 = 27460
index.list.p13 <- lapply(1:5, function(i){
  pn=13
  list.index <- get(paste0("index.list.R.p", pn-7))
  list(index.list.Q.S2[[i]], list.index[[i]])
})

# P14 = 0.1:1, reference = 1552, query = 137303/5 = 27460
index.list.p14 <- lapply(1:5, function(i){
  pn=14
  list.index <- get(paste0("index.list.R.p", pn-7))
  list(index.list.Q.S2[[i]], list.index[[i]])
})

for (i in 1:14) {
  saveRDS(get(paste0("index.list.p", i)), file = paste0("08_StartOver_NewStory/05_interdatasets/03_brain/01_datasets/index/Index_P", i, ".rds"))
}







