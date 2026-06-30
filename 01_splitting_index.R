setwd("R:/ANNCOM2023-Q6025")
library(Seurat)
library(tidyverse)
library(dplyr)

# Session 1 ######
# cell-type proportions change in query and reference together and cell counts remains the same each iteration
# 2794 cells each iteration = 1863*(0.5+1), so total cells will be 2794*5=13970
# P1=10:1, P2=5:1, P3=2:1, P4=1:1, P5=0.5:1
combined <- readRDS("R:/ANNCOM2023-Q6025/08_StartOver_NewStory/01_datasets/seuobj/combined.p1.rds")
table(combined$bm_cell_type)
table(combined$bm_general)
combined$index <- 1:ncol(combined)
meta.ori <- combined@meta.data[,c("index", "bm_cell_type", "bm_general")]

## P1 ####
meta <- combined@meta.data[,c("index", "bm_cell_type", "bm_general")]
epi <- meta[meta$bm_general=="EPITHELIAL",]
imm <- meta[meta$bm_general=="IMMUNE",]
# set.seed(233) # maybe remove this or change every time
epi <- epi %>% group_by(bm_cell_type) %>% slice_sample(prop = 10*nrow(imm)/nrow(epi))
meta <- rbind(epi, imm)
# set.seed(233)
testindex.p1 <- meta %>% group_by(bm_cell_type) %>% slice_sample(prop = 13970/nrow(meta))
index <- rBayesianOptimization::KFold(target = testindex.p1$bm_cell_type, nfolds = 5, stratified = TRUE, seed = 233)

real.index.p1.query <- lapply(index, function(ind){
  testindex.p1[ind,]$index
})

real.index.p1.ref <- lapply(seq_along(index), function(i){
  rest <- meta.ori[meta.ori$index %in% testindex.p1[-index[[i]],]$index,]
  yy <- rest %>% group_by(bm_cell_type) %>% slice_sample(prop = (nrow(testindex.p1[index[[i]],])+5)/nrow(rest))
  yy$index
})

real.index.p1 <- lapply(seq_along(real.index.p1.query), function(i){
  list(real.index.p1.query[[i]], real.index.p1.ref[[i]])
})

## P2 #####
meta <- combined@meta.data[,c("index", "bm_cell_type", "bm_general")]
epi <- meta[meta$bm_general=="EPITHELIAL",]
imm <- meta[meta$bm_general=="IMMUNE",]
# set.seed(233) # maybe remove this or change every time
epi <- epi %>% group_by(bm_cell_type) %>% slice_sample(prop = 5*nrow(imm)/nrow(epi)) # change here for every proportion
meta <- rbind(epi, imm)
# set.seed(233)
testindex.p2 <- meta %>% group_by(bm_cell_type) %>% slice_sample(prop = 13970/nrow(meta))
index <- rBayesianOptimization::KFold(target = testindex.p2$bm_cell_type, nfolds = 5, stratified = TRUE, seed = 233)

real.index.p2.query <- lapply(index, function(ind){
  testindex.p2[ind,]$index
})

real.index.p2.ref <- lapply(seq_along(index), function(i){
  rest <- meta.ori[meta.ori$index %in% testindex.p2[-index[[i]],]$index,]
  yy <- rest %>% group_by(bm_cell_type) %>% slice_sample(prop = (nrow(testindex.p2[index[[i]],])+5)/nrow(rest))
  yy$index
})

real.index.p2 <- lapply(seq_along(real.index.p2.query), function(i){
  list(real.index.p2.query[[i]], real.index.p2.ref[[i]])
})


## P3 ######
meta <- combined@meta.data[,c("index", "bm_cell_type", "bm_general")]
epi <- meta[meta$bm_general=="EPITHELIAL",]
imm <- meta[meta$bm_general=="IMMUNE",]
# set.seed(233) # maybe remove this or change every time
epi <- epi %>% group_by(bm_cell_type) %>% slice_sample(prop = 2*nrow(imm)/nrow(epi)) # change here for every proportion
meta <- rbind(epi, imm)
# set.seed(233)
testindex.p3 <- meta %>% group_by(bm_cell_type) %>% slice_sample(prop = 13970/nrow(meta))
index <- rBayesianOptimization::KFold(target = testindex.p3$bm_cell_type, nfolds = 5, stratified = TRUE, seed = 233)

real.index.p3.query <- lapply(index, function(ind){
  testindex.p3[ind,]$index
})

real.index.p3.ref <- lapply(seq_along(index), function(i){
  rest <- meta.ori[meta.ori$index %in% testindex.p3[-index[[i]],]$index,]
  yy <- rest %>% group_by(bm_cell_type) %>% slice_sample(prop = (nrow(testindex.p3[index[[i]],])+5)/nrow(rest))
  yy$index
})

real.index.p3 <- lapply(seq_along(real.index.p3.query), function(i){
  list(real.index.p3.query[[i]], real.index.p3.ref[[i]])
})


## P4 ######
meta <- combined@meta.data[,c("index", "bm_cell_type", "bm_general")]
epi <- meta[meta$bm_general=="EPITHELIAL",]
imm <- meta[meta$bm_general=="IMMUNE",]
# set.seed(233) # maybe remove this or change every time
epi <- epi %>% group_by(bm_cell_type) %>% slice_sample(prop = 1*nrow(imm)/nrow(epi)) # change here for every proportion
meta <- rbind(epi, imm)
# set.seed(233)
testindex.p4 <- meta %>% group_by(bm_cell_type) %>% slice_sample(prop = 13970/nrow(meta))
index <- rBayesianOptimization::KFold(target = testindex.p4$bm_cell_type, nfolds = 5, stratified = TRUE, seed = 233)

real.index.p4.query <- lapply(index, function(ind){
  testindex.p4[ind,]$index
})

real.index.p4.ref <- lapply(seq_along(index), function(i){
  rest <- meta.ori[meta.ori$index %in% testindex.p4[-index[[i]],]$index,]
  yy <- rest %>% group_by(bm_cell_type) %>% slice_sample(prop = (nrow(testindex.p4[index[[i]],])+5)/nrow(rest))
  yy$index
})

real.index.p4 <- lapply(seq_along(real.index.p4.query), function(i){
  list(real.index.p4.query[[i]], real.index.p4.ref[[i]])
})


## P5 ######
meta <- combined@meta.data[,c("index", "bm_cell_type", "bm_general")]
epi <- meta[meta$bm_general=="EPITHELIAL",]
imm <- meta[meta$bm_general=="IMMUNE",]
# set.seed(233) # maybe remove this or change every time
epi <- epi %>% group_by(bm_cell_type) %>% slice_sample(prop = 0.5*nrow(imm)/nrow(epi)) # change here for every proportion
meta <- rbind(epi, imm)
# set.seed(233)
testindex.p5 <- meta %>% group_by(bm_cell_type) %>% slice_sample(prop = 13970/nrow(meta))
index <- rBayesianOptimization::KFold(target = testindex.p5$bm_cell_type, nfolds = 5, stratified = TRUE, seed = 233)

real.index.p5.query <- lapply(index, function(ind){
  testindex.p5[ind,]$index
})

real.index.p5.ref <- lapply(seq_along(index), function(i){
  rest <- meta.ori[meta.ori$index %in% testindex.p5[-index[[i]],]$index,]
  yy <- rest %>% group_by(bm_cell_type) %>% slice_sample(prop = (nrow(testindex.p5[index[[i]],])+5)/nrow(rest))
  yy$index
})

real.index.p5 <- lapply(seq_along(real.index.p5.query), function(i){
  list(real.index.p5.query[[i]], real.index.p5.ref[[i]])
})

for (i in 1:5) {
  saveRDS(get(paste0("real.index.p", i)), file = paste0("08_StartOver_NewStory/01_datasets/index/Index_P", i, ".rds"))
}


# Session 2: ####################
# keep 1:1 in query and varying proportions in reference but keep counts of query and reference and every settings the same
combined <- readRDS("R:/ANNCOM2023-Q6025/08_StartOver_NewStory/01_datasets/seuobj/combined.p1.rds")
table(combined$bm_cell_type)
table(combined$bm_general)
combined$index <- 1:ncol(combined)

# set query
meta <- combined@meta.data[,c("index", "bm_cell_type", "bm_general")]
imm <- meta[meta$bm_general=="IMMUNE",]
epi <- meta[meta$bm_general=="EPITHELIAL",]
epi <- epi %>% group_by(bm_cell_type) %>% slice_sample(prop = nrow(imm)/nrow(epi))
meta.Q <- rbind(epi, imm)
index.temp <- rBayesianOptimization::KFold(target = meta.Q$bm_cell_type, nfolds = 5, stratified = TRUE, seed = 233)

index.Q.real <- lapply(index.temp, function(ind){
  meta.Q[ind,]$index
})

index.R.pool <- lapply(index.Q.real, function(ind){
  meta[-ind,]$index
})

# 3731 3731 3731 3731 3709

## P6 #####
# P6 = 10:1
index.R.p6 <- lapply(index.R.pool, function(ind){
  p=10
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*(3731+5)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*(3731+5)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p6[[1]],]$bm_general)

real.index.p6 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p6[[ind]])
})

## P7 #####
# P7 = 5:1
index.R.p7 <- lapply(index.R.pool, function(ind){
  p=5
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*(3731+5)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*(3731+5)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p7[[1]],]$bm_general)

real.index.p7 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p7[[ind]])
})

## P8 #####
# P8 = 2:1
index.R.p8 <- lapply(index.R.pool, function(ind){
  p=2
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*(3731+5)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*(3731+5)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p8[[1]],]$bm_general)

real.index.p8 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p8[[ind]])
})

## P9 #####
# P9 = 1:1
index.R.p9 <- lapply(index.R.pool, function(ind){
  p=1
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*(3731+5)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*(3731+5)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p9[[1]],]$bm_general)

real.index.p9 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p9[[ind]])
})

## P10 #####
# P10 = 0.5:1
index.R.p10 <- lapply(index.R.pool, function(ind){
  p=0.5
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*(3731+5)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*(3731+5)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p10[[1]],]$bm_general)

real.index.p10 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p10[[ind]])
})

## P11 #####
# P11 = 0.2:1
index.R.p11 <- lapply(index.R.pool, function(ind){
  p=0.2
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*(3731+5)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*(3731+5)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p11[[1]],]$bm_general)

real.index.p11 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p11[[ind]])
})

## P12 #####
# P12 = 0.1:1
index.R.p12 <- lapply(index.R.pool, function(ind){
  p=0.1
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*(3731+5)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*(3731+5)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p12[[1]],]$bm_general)

real.index.p12 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p12[[ind]])
})

for (i in 6:12) {
  saveRDS(get(paste0("real.index.p", i)), file = paste0("08_StartOver_NewStory/01_datasets/index/Index_P", i, ".rds"))
}

# Session 3:#####
# increase counts in reference:
# P13-P19 1.5 times of Query, P20-P26 2 times of Query

combined <- readRDS("R:/ANNCOM2023-Q6025/08_StartOver_NewStory/01_datasets/seuobj/combined.p1.rds")
table(combined$bm_cell_type)
table(combined$bm_general)
combined$index <- 1:ncol(combined)

# set query
meta <- combined@meta.data[,c("index", "bm_cell_type", "bm_general")]
imm <- meta[meta$bm_general=="IMMUNE",]
epi <- meta[meta$bm_general=="EPITHELIAL",]
epi <- epi %>% group_by(bm_cell_type) %>% slice_sample(prop = nrow(imm)/nrow(epi))
meta.Q <- rbind(epi, imm)
index.temp <- rBayesianOptimization::KFold(target = meta.Q$bm_cell_type, nfolds = 5, stratified = TRUE, seed = 233)

index.Q.real <- lapply(index.temp, function(ind){
  meta.Q[ind,]$index
})

index.R.pool <- lapply(index.Q.real, function(ind){
  meta[-ind,]$index
})

# P13 = 10:1, 1.5 times
index.R.p13 <- lapply(index.R.pool, function(ind){
  p=10
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*((3731+4)*1.5)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*((3731+4)*1.5)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p13[[1]],]$bm_general)

real.index.p13 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p13[[ind]])
})


# P14 = 5:1, 1.5 times
index.R.p14 <- lapply(index.R.pool, function(ind){
  p=5
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*((3731+4)*1.5)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*((3731+4)*1.5)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p14[[1]],]$bm_general)

real.index.p14 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p14[[ind]])
})


# P15 = 2:1, 1.5 times
index.R.p15 <- lapply(index.R.pool, function(ind){
  p=2
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*((3731+4)*1.5)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*((3731+4)*1.5)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p15[[1]],]$bm_general)

real.index.p15 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p15[[ind]])
})

# P16 = 1:1, 1.5 times
index.R.p16 <- lapply(index.R.pool, function(ind){
  p=1
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*((3731+4)*1.5)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*((3731+4)*1.5)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p16[[1]],]$bm_general)

real.index.p16 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p16[[ind]])
})

# P17 = 0.5:1, 1.5 times
index.R.p17 <- lapply(index.R.pool, function(ind){
  p=0.5
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*((3731+4)*1.5)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*((3731+4)*1.5)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p17[[1]],]$bm_general)

real.index.p17 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p17[[ind]])
})

# P18 = 0.2:1, 1.5 times
index.R.p18 <- lapply(index.R.pool, function(ind){
  p=0.2
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*((3731+4)*1.5)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*((3731+4)*1.5)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p18[[1]],]$bm_general)

real.index.p18 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p18[[ind]])
})


# P19 = 0.1:1, 1.5 times
index.R.p19 <- lapply(index.R.pool, function(ind){
  p=0.1
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*((3731+4)*1.5)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*((3731+4)*1.5)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p19[[1]],]$bm_general)

real.index.p19 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p19[[ind]])
})

# P20 = 10:1, 2 times
index.R.p20 <- lapply(index.R.pool, function(ind){
  p=10
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*((3731+4)*2)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*((3731+4)*2)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p20[[1]],]$bm_general)

real.index.p20 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p20[[ind]])
})

# P21 = 5:1, 2 times
index.R.p21 <- lapply(index.R.pool, function(ind){
  p=5
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*((3731+4)*2)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*((3731+4)*2)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p21[[1]],]$bm_general)

real.index.p21 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p21[[ind]])
})

# P22 = 2:1, 2 times
index.R.p22 <- lapply(index.R.pool, function(ind){
  p=2
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*((3731+4)*2)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*((3731+4)*2)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p22[[1]],]$bm_general)

real.index.p22 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p22[[ind]])
})

# P23 = 1:1, 2 times
index.R.p23 <- lapply(index.R.pool, function(ind){
  p=1
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*((3731+4)*2)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*((3731+4)*2)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p23[[1]],]$bm_general)

real.index.p23 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p23[[ind]])
})

# P24 = 0.5:1, 2 times
index.R.p24 <- lapply(index.R.pool, function(ind){
  p=0.5
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*((3731+4)*2)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*((3731+4)*2)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p24[[1]],]$bm_general)

real.index.p24 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p24[[ind]])
})

# P25 = 0.2:1, 2 times
index.R.p25 <- lapply(index.R.pool, function(ind){
  p=0.2
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*((3731+4)*2)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*((3731+4)*2)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p25[[1]],]$bm_general)

real.index.p25 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p25[[ind]])
})

# P26 = 0.1:1, 2 times
index.R.p26 <- lapply(index.R.pool, function(ind){
  p=0.1
  ref.pool <- meta[ind,]
  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]
  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = p/(p+1)*((3731+4)*2)/nrow(ref.pool.epi))
  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1/(p+1))*((3731+4)*2)/nrow(ref.pool.imm))
  ref <- rbind(ref.epi, ref.imm)
  index <- ref$index
})

table(meta[index.R.p26[[1]],]$bm_general)

real.index.p26 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p26[[ind]])
})

for (i in 13:26) {
  saveRDS(get(paste0("real.index.p", i)), file = paste0("08_StartOver_NewStory/01_datasets/index/Index_P", i, ".rds"))
}

# Session 4: ######
# Subtype ratios varying
combined <- readRDS("R:/ANNCOM2023-Q6025/08_StartOver_NewStory/01_datasets/seuobj/combined.p1.rds")
table(combined$bm_cell_type)
table(combined$bm_general)
combined$index <- 1:ncol(combined)

# set query
meta <- combined@meta.data[,c("index", "bm_cell_type", "bm_general")]
imm <- meta[meta$bm_general=="IMMUNE",]
epi <- meta[meta$bm_general=="EPITHELIAL",]
epi <- epi %>% group_by(bm_cell_type) %>% slice_sample(prop = nrow(imm)/nrow(epi))
meta.Q <- rbind(epi, imm)
index.temp <- rBayesianOptimization::KFold(target = meta.Q$bm_cell_type, nfolds = 5, stratified = TRUE, seed = 233)

index.Q.real <- lapply(index.temp, function(ind){
  meta.Q[ind,]$index
})

index.R.pool <- lapply(index.Q.real, function(ind){
  meta[-ind,]$index
})

yy.e <- c(76652, 19604, 503, 493)
yy.i <- c(6235, 1427, 903, 626, 64, 62)
names(yy.e) <- c("LC2", "LC1", "basal", "cyc")
names(yy.i) <- c("mac", "T", "neutro", "DC", "NK", "B")

yy.e/97252*1865
yy.i/9317*1866

# so counts for cycling=9 and basal=10 when ratios changing within epithelial
# counts for neutrophil=181, DC=125, NK=13, B=12
# LC2:LC1 = 3.910018, mac:T = 4.369306
# LC1+LC2 = 1846, mac+T = 1535

## P27 = 4:1 (LC2:LC1) and rest keep consistent in immune ######
index.R.p27 <- lapply(index.R.pool, function(ind){
  p=4
  ref.pool <- meta[ind,]

  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]

  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1866+3)/nrow(ref.pool.imm))

  ref.basal <- ref.pool.epi[ref.pool.epi$bm_cell_type=="Basal",] %>% slice_sample(n=10)
  ref.cyc <- ref.pool.epi[ref.pool.epi$bm_cell_type=="G2MS",] %>% slice_sample(n=11)
  ref.LC2 <- ref.pool.epi[ref.pool.epi$bm_cell_type=="LC2",] %>% slice_sample(n=round(1846*p/(p+1)))
  ref.LC1 <- ref.pool.epi[ref.pool.epi$bm_cell_type=="LC1",] %>% slice_sample(n=round(1846*1/(p+1)))
  ref <- do.call("rbind", list(ref.imm, ref.basal, ref.cyc, ref.LC2, ref.LC1))
  index <- ref$index

})

table(meta[index.R.p27[[1]],]$bm_general)

real.index.p27 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p27[[ind]])
})

## P28 = 2.5:1 (LC2:LC1) and rest keep consistent in immune ######
index.R.p28 <- lapply(index.R.pool, function(ind){
  p=2.5
  ref.pool <- meta[ind,]

  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]

  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1866+3)/nrow(ref.pool.imm))

  ref.basal <- ref.pool.epi[ref.pool.epi$bm_cell_type=="Basal",] %>% slice_sample(n=10)
  ref.cyc <- ref.pool.epi[ref.pool.epi$bm_cell_type=="G2MS",] %>% slice_sample(n=11)
  ref.LC2 <- ref.pool.epi[ref.pool.epi$bm_cell_type=="LC2",] %>% slice_sample(n=round(1846*p/(p+1)))
  ref.LC1 <- ref.pool.epi[ref.pool.epi$bm_cell_type=="LC1",] %>% slice_sample(n=round(1846*1/(p+1)))
  ref <- do.call("rbind", list(ref.imm, ref.basal, ref.cyc, ref.LC2, ref.LC1))
  index <- ref$index

})

table(meta[index.R.p28[[1]],]$bm_general)

real.index.p28 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p28[[ind]])
})

## p29 = 2:1 (LC2:LC1) and rest keep consistent in immune ######
index.R.p29 <- lapply(index.R.pool, function(ind){
  p=2
  ref.pool <- meta[ind,]

  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]

  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1866+3)/nrow(ref.pool.imm))

  ref.basal <- ref.pool.epi[ref.pool.epi$bm_cell_type=="Basal",] %>% slice_sample(n=10)
  ref.cyc <- ref.pool.epi[ref.pool.epi$bm_cell_type=="G2MS",] %>% slice_sample(n=11)
  ref.LC2 <- ref.pool.epi[ref.pool.epi$bm_cell_type=="LC2",] %>% slice_sample(n=round(1846*p/(p+1)))
  ref.LC1 <- ref.pool.epi[ref.pool.epi$bm_cell_type=="LC1",] %>% slice_sample(n=round(1846*1/(p+1)))
  ref <- do.call("rbind", list(ref.imm, ref.basal, ref.cyc, ref.LC2, ref.LC1))
  index <- ref$index

})

table(meta[index.R.p29[[1]],]$bm_general)

real.index.p29 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p29[[ind]])
})

## p30 = 1:1 (LC2:LC1) and rest keep consistent in immune ######
index.R.p30 <- lapply(index.R.pool, function(ind){
  p=1
  ref.pool <- meta[ind,]

  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]

  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1866+3)/nrow(ref.pool.imm))

  ref.basal <- ref.pool.epi[ref.pool.epi$bm_cell_type=="Basal",] %>% slice_sample(n=10)
  ref.cyc <- ref.pool.epi[ref.pool.epi$bm_cell_type=="G2MS",] %>% slice_sample(n=11)
  ref.LC2 <- ref.pool.epi[ref.pool.epi$bm_cell_type=="LC2",] %>% slice_sample(n=round(1846*p/(p+1)))
  ref.LC1 <- ref.pool.epi[ref.pool.epi$bm_cell_type=="LC1",] %>% slice_sample(n=round(1846*1/(p+1)))
  ref <- do.call("rbind", list(ref.imm, ref.basal, ref.cyc, ref.LC2, ref.LC1))
  index <- ref$index

})

table(meta[index.R.p30[[1]],]$bm_general)

real.index.p30 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p30[[ind]])
})

## p31 = 0.5:1 (LC2:LC1) and rest keep consistent in immune ######
index.R.p31 <- lapply(index.R.pool, function(ind){
  p=0.5
  ref.pool <- meta[ind,]

  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]

  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1866+3)/nrow(ref.pool.imm))

  ref.basal <- ref.pool.epi[ref.pool.epi$bm_cell_type=="Basal",] %>% slice_sample(n=10)
  ref.cyc <- ref.pool.epi[ref.pool.epi$bm_cell_type=="G2MS",] %>% slice_sample(n=11)
  ref.LC2 <- ref.pool.epi[ref.pool.epi$bm_cell_type=="LC2",] %>% slice_sample(n=round(1846*p/(p+1)))
  ref.LC1 <- ref.pool.epi[ref.pool.epi$bm_cell_type=="LC1",] %>% slice_sample(n=round(1846*1/(p+1)))
  ref <- do.call("rbind", list(ref.imm, ref.basal, ref.cyc, ref.LC2, ref.LC1))
  index <- ref$index

})

table(meta[index.R.p31[[1]],]$bm_general)

real.index.p31 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p31[[ind]])
})

## p32 = 0.4:1 (LC2:LC1) and rest keep consistent in immune ######
index.R.p32 <- lapply(index.R.pool, function(ind){
  p=0.4
  ref.pool <- meta[ind,]

  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]

  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1866+3)/nrow(ref.pool.imm))

  ref.basal <- ref.pool.epi[ref.pool.epi$bm_cell_type=="Basal",] %>% slice_sample(n=10)
  ref.cyc <- ref.pool.epi[ref.pool.epi$bm_cell_type=="G2MS",] %>% slice_sample(n=11)
  ref.LC2 <- ref.pool.epi[ref.pool.epi$bm_cell_type=="LC2",] %>% slice_sample(n=round(1846*p/(p+1)))
  ref.LC1 <- ref.pool.epi[ref.pool.epi$bm_cell_type=="LC1",] %>% slice_sample(n=round(1846*1/(p+1)))
  ref <- do.call("rbind", list(ref.imm, ref.basal, ref.cyc, ref.LC2, ref.LC1))
  index <- ref$index

})

table(meta[index.R.p32[[1]],]$bm_general)

real.index.p32 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p32[[ind]])
})

## p33 = 0.25:1 (LC2:LC1) and rest keep consistent in immune ######
index.R.p33 <- lapply(index.R.pool, function(ind){
  p=0.25
  ref.pool <- meta[ind,]

  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]

  ref.imm <- ref.pool.imm %>% group_by(bm_cell_type) %>% slice_sample(prop = (1866+3)/nrow(ref.pool.imm))

  ref.basal <- ref.pool.epi[ref.pool.epi$bm_cell_type=="Basal",] %>% slice_sample(n=10)
  ref.cyc <- ref.pool.epi[ref.pool.epi$bm_cell_type=="G2MS",] %>% slice_sample(n=11)
  ref.LC2 <- ref.pool.epi[ref.pool.epi$bm_cell_type=="LC2",] %>% slice_sample(n=round(1846*p/(p+1)))
  ref.LC1 <- ref.pool.epi[ref.pool.epi$bm_cell_type=="LC1",] %>% slice_sample(n=round(1846*1/(p+1)))
  ref <- do.call("rbind", list(ref.imm, ref.basal, ref.cyc, ref.LC2, ref.LC1))
  index <- ref$index

})

table(meta[index.R.p33[[1]],]$bm_general)

real.index.p33 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p33[[ind]])
})

## P34 = 4:1 (Mac:T) and rest keep consistent in epithelial ######
index.R.p34 <- lapply(index.R.pool, function(ind){
  p=4
  ref.pool <- meta[ind,]

  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]

  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = (1866+2)/nrow(ref.pool.epi))

  ref.mac <- ref.pool.imm[ref.pool.imm$bm_cell_type=="Macrophage",] %>% slice_sample(n=round(1535*p/(1+p)))
  ref.T <- ref.pool.imm[ref.pool.imm$bm_cell_type=="T_cell",] %>% slice_sample(n=round(1535*1/(1+p)))

  ref.neu <- ref.pool.imm[ref.pool.imm$bm_cell_type=="Neutrophil",] %>% slice_sample(n=181)
  ref.B <- ref.pool.imm[ref.pool.imm$bm_cell_type=="B_cell",] %>% slice_sample(n=12)
  ref.DC <- ref.pool.imm[ref.pool.imm$bm_cell_type=="DC",] %>% slice_sample(n=125)
  ref.NK <- ref.pool.imm[ref.pool.imm$bm_cell_type=="NK",] %>% slice_sample(n=13)

  ref <- do.call("rbind", list(ref.epi, ref.mac, ref.T, ref.neu, ref.B, ref.DC, ref.NK))
  index <- ref$index

})

table(meta[index.R.p34[[1]],]$bm_general)

real.index.p34 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p34[[ind]])
})


## p35 = 2.5:1 (Mac:T) and rest keep consistent in epithelial ######
index.R.p35 <- lapply(index.R.pool, function(ind){
  p=2.5
  ref.pool <- meta[ind,]

  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]

  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = (1866+2)/nrow(ref.pool.epi))

  ref.mac <- ref.pool.imm[ref.pool.imm$bm_cell_type=="Macrophage",] %>% slice_sample(n=round(1535*p/(1+p)))
  ref.T <- ref.pool.imm[ref.pool.imm$bm_cell_type=="T_cell",] %>% slice_sample(n=round(1535*1/(1+p)))

  ref.neu <- ref.pool.imm[ref.pool.imm$bm_cell_type=="Neutrophil",] %>% slice_sample(n=181)
  ref.B <- ref.pool.imm[ref.pool.imm$bm_cell_type=="B_cell",] %>% slice_sample(n=12)
  ref.DC <- ref.pool.imm[ref.pool.imm$bm_cell_type=="DC",] %>% slice_sample(n=125)
  ref.NK <- ref.pool.imm[ref.pool.imm$bm_cell_type=="NK",] %>% slice_sample(n=13)

  ref <- do.call("rbind", list(ref.epi, ref.mac, ref.T, ref.neu, ref.B, ref.DC, ref.NK))
  index <- ref$index

})

table(meta[index.R.p35[[1]],]$bm_general)

real.index.p35 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p35[[ind]])
})

## p36 = 2:1 (Mac:T) and rest keep consistent in epithelial ######
index.R.p36 <- lapply(index.R.pool, function(ind){
  p=2
  ref.pool <- meta[ind,]

  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]

  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = (1866+2)/nrow(ref.pool.epi))

  ref.mac <- ref.pool.imm[ref.pool.imm$bm_cell_type=="Macrophage",] %>% slice_sample(n=round(1535*p/(1+p)))
  ref.T <- ref.pool.imm[ref.pool.imm$bm_cell_type=="T_cell",] %>% slice_sample(n=round(1535*1/(1+p)))

  ref.neu <- ref.pool.imm[ref.pool.imm$bm_cell_type=="Neutrophil",] %>% slice_sample(n=181)
  ref.B <- ref.pool.imm[ref.pool.imm$bm_cell_type=="B_cell",] %>% slice_sample(n=12)
  ref.DC <- ref.pool.imm[ref.pool.imm$bm_cell_type=="DC",] %>% slice_sample(n=125)
  ref.NK <- ref.pool.imm[ref.pool.imm$bm_cell_type=="NK",] %>% slice_sample(n=13)

  ref <- do.call("rbind", list(ref.epi, ref.mac, ref.T, ref.neu, ref.B, ref.DC, ref.NK))
  index <- ref$index

})

table(meta[index.R.p36[[1]],]$bm_general)

real.index.p36 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p36[[ind]])
})

## p37 = 1:1 (Mac:T) and rest keep consistent in epithelial ######
index.R.p37 <- lapply(index.R.pool, function(ind){
  p=1
  ref.pool <- meta[ind,]

  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]

  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = (1866+2)/nrow(ref.pool.epi))

  ref.mac <- ref.pool.imm[ref.pool.imm$bm_cell_type=="Macrophage",] %>% slice_sample(n=round(1535*p/(1+p)))
  ref.T <- ref.pool.imm[ref.pool.imm$bm_cell_type=="T_cell",] %>% slice_sample(n=round(1535*1/(1+p)))

  ref.neu <- ref.pool.imm[ref.pool.imm$bm_cell_type=="Neutrophil",] %>% slice_sample(n=181)
  ref.B <- ref.pool.imm[ref.pool.imm$bm_cell_type=="B_cell",] %>% slice_sample(n=12)
  ref.DC <- ref.pool.imm[ref.pool.imm$bm_cell_type=="DC",] %>% slice_sample(n=125)
  ref.NK <- ref.pool.imm[ref.pool.imm$bm_cell_type=="NK",] %>% slice_sample(n=13)

  ref <- do.call("rbind", list(ref.epi, ref.mac, ref.T, ref.neu, ref.B, ref.DC, ref.NK))
  index <- ref$index

})

table(meta[index.R.p37[[1]],]$bm_general)

real.index.p37 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p37[[ind]])
})

## p38 = 0.5:1 (Mac:T) and rest keep consistent in epithelial ######
index.R.p38 <- lapply(index.R.pool, function(ind){
  p=0.5
  ref.pool <- meta[ind,]

  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]

  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = (1866+2)/nrow(ref.pool.epi))

  ref.mac <- ref.pool.imm[ref.pool.imm$bm_cell_type=="Macrophage",] %>% slice_sample(n=round(1535*p/(1+p)))
  ref.T <- ref.pool.imm[ref.pool.imm$bm_cell_type=="T_cell",] %>% slice_sample(n=round(1535*1/(1+p)))

  ref.neu <- ref.pool.imm[ref.pool.imm$bm_cell_type=="Neutrophil",] %>% slice_sample(n=181)
  ref.B <- ref.pool.imm[ref.pool.imm$bm_cell_type=="B_cell",] %>% slice_sample(n=12)
  ref.DC <- ref.pool.imm[ref.pool.imm$bm_cell_type=="DC",] %>% slice_sample(n=125)
  ref.NK <- ref.pool.imm[ref.pool.imm$bm_cell_type=="NK",] %>% slice_sample(n=13)

  ref <- do.call("rbind", list(ref.epi, ref.mac, ref.T, ref.neu, ref.B, ref.DC, ref.NK))
  index <- ref$index

})

table(meta[index.R.p38[[1]],]$bm_general)

real.index.p38 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p38[[ind]])
})

## p39 = 0.4:1 (Mac:T) and rest keep consistent in epithelial ######
index.R.p39 <- lapply(index.R.pool, function(ind){
  p=0.4
  ref.pool <- meta[ind,]

  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]

  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = (1866+2)/nrow(ref.pool.epi))

  ref.mac <- ref.pool.imm[ref.pool.imm$bm_cell_type=="Macrophage",] %>% slice_sample(n=round(1535*p/(1+p)))
  ref.T <- ref.pool.imm[ref.pool.imm$bm_cell_type=="T_cell",] %>% slice_sample(n=round(1535*1/(1+p)))

  ref.neu <- ref.pool.imm[ref.pool.imm$bm_cell_type=="Neutrophil",] %>% slice_sample(n=181)
  ref.B <- ref.pool.imm[ref.pool.imm$bm_cell_type=="B_cell",] %>% slice_sample(n=12)
  ref.DC <- ref.pool.imm[ref.pool.imm$bm_cell_type=="DC",] %>% slice_sample(n=125)
  ref.NK <- ref.pool.imm[ref.pool.imm$bm_cell_type=="NK",] %>% slice_sample(n=13)

  ref <- do.call("rbind", list(ref.epi, ref.mac, ref.T, ref.neu, ref.B, ref.DC, ref.NK))
  index <- ref$index

})

table(meta[index.R.p39[[1]],]$bm_general)

real.index.p39 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p39[[ind]])
})


## p40 = 0.25:1 (Mac:T) and rest keep consistent in epithelial ######
index.R.p40 <- lapply(index.R.pool, function(ind){
  p=0.25
  ref.pool <- meta[ind,]

  ref.pool.epi <- ref.pool[ref.pool$bm_general=="EPITHELIAL",]
  ref.pool.imm <- ref.pool[ref.pool$bm_general=="IMMUNE",]

  ref.epi <- ref.pool.epi %>% group_by(bm_cell_type) %>% slice_sample(prop = (1866+2)/nrow(ref.pool.epi))

  ref.mac <- ref.pool.imm[ref.pool.imm$bm_cell_type=="Macrophage",] %>% slice_sample(n=round(1535*p/(1+p)))
  ref.T <- ref.pool.imm[ref.pool.imm$bm_cell_type=="T_cell",] %>% slice_sample(n=round(1535*1/(1+p)))

  ref.neu <- ref.pool.imm[ref.pool.imm$bm_cell_type=="Neutrophil",] %>% slice_sample(n=181)
  ref.B <- ref.pool.imm[ref.pool.imm$bm_cell_type=="B_cell",] %>% slice_sample(n=12)
  ref.DC <- ref.pool.imm[ref.pool.imm$bm_cell_type=="DC",] %>% slice_sample(n=125)
  ref.NK <- ref.pool.imm[ref.pool.imm$bm_cell_type=="NK",] %>% slice_sample(n=13)

  ref <- do.call("rbind", list(ref.epi, ref.mac, ref.T, ref.neu, ref.B, ref.DC, ref.NK))
  index <- ref$index

})

table(meta[index.R.p40[[1]],]$bm_general)

real.index.p40 <- lapply(seq_along(index.Q.real), function(ind){
  list(index.Q.real[[ind]], index.R.p40[[ind]])
})

for (i in 27:40) {
  saveRDS(get(paste0("real.index.p", i)), file = paste0("08_StartOver_NewStory/01_datasets/index/Index_P", i, ".rds"))
}

















