setwd("R:/ANNCOM2023-Q6025")
library(Seurat)
library(dplyr)
library(tidyverse)

seg <- readRDS(file="08_StartOver_NewStory/05_interdatasets/02_pancreas/01_datasets/seuobj/pancreas_seg.rds")
baron <- readRDS(file="08_StartOver_NewStory/05_interdatasets/02_pancreas/01_datasets/seuobj/pancreas_baron.rds")

seg@meta.data <- seg@meta.data %>% mutate(cell_type_3 = case_when(cell_type_2 == "Endocrine" ~ "Endocrine",
                                                                  cell_type_2 != "Endocrine" ~ "Non-Endocrine"))

baron@meta.data <- baron@meta.data %>% mutate(cell_type_3 = case_when(cell_type_2 == "Endocrine" ~ "Endocrine",
                                                                      cell_type_2 != "Endocrine" ~ "Non-Endocrine"))

seg@meta.data$index <- 1:ncol(seg)
baron@meta.data$index <- 1:ncol(baron)

saveRDS(seg, file = "08_StartOver_NewStory/05_interdatasets/02_pancreas/01_datasets/seuobj/pancreas_seg.rds")
saveRDS(baron, file = "08_StartOver_NewStory/05_interdatasets/02_pancreas/01_datasets/seuobj/pancreas_baron.rds")

# baron health will be reference and seg T2D will be query
# query = 2209 (1554 + 655), reference =  8569 (5725 + 2844)

table(seg$cell_type_1)/sum(table(seg$cell_type_1))*100
table(seg$cell_type_2)/sum(table(seg$cell_type_2))*100
table(seg$cell_type_3)/sum(table(seg$cell_type_3))*100


table(baron$cell_type_1)/sum(table(baron$cell_type_1))*100
table(baron$cell_type_2)/sum(table(baron$cell_type_2))*100
table(baron$cell_type_3)/sum(table(baron$cell_type_3))*100

# Session 1: query count = reference count = 2209/5 = 441 per iteration
t2d <- seg@meta.data
health <- baron@meta.data

index <- rBayesianOptimization::KFold(target = t2d$cell_type_1, nfolds = 5, stratified = TRUE, seed = 233)
query.index.s1 <- lapply(index, function(ind){
  t2d[ind,]$index
})

# P1 = 10:1
ref.index.p1 <- lapply(1:5, function(i){
  p=10
  endo <- health[health$cell_type_3=="Endocrine",]
  nonendo <- health[health$cell_type_3=="Non-Endocrine",]
  endo.ref <- endo %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(441+5)/nrow(endo))
  nonendo.ref <- nonendo %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(441+5)/nrow(nonendo))
  ref <- rbind(endo.ref, nonendo.ref)
  ref$index
})
table(health[ref.index.p1[[1]],]$cell_type_3)

index.p1 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], ref.index.p1[[i]])
})

# P2 = 5:1
ref.index.p2 <- lapply(1:5, function(i){
  p=5
  endo <- health[health$cell_type_3=="Endocrine",]
  nonendo <- health[health$cell_type_3=="Non-Endocrine",]
  endo.ref <- endo %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(441+5)/nrow(endo))
  nonendo.ref <- nonendo %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(441+5)/nrow(nonendo))
  ref <- rbind(endo.ref, nonendo.ref)
  ref$index
})
table(health[ref.index.p2[[1]],]$cell_type_3)

index.p2 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], ref.index.p2[[i]])
})

# P3 = 2:1
ref.index.p3 <- lapply(1:5, function(i){
  p=2
  endo <- health[health$cell_type_3=="Endocrine",]
  nonendo <- health[health$cell_type_3=="Non-Endocrine",]
  endo.ref <- endo %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(441+5)/nrow(endo))
  nonendo.ref <- nonendo %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(441+5)/nrow(nonendo))
  ref <- rbind(endo.ref, nonendo.ref)
  ref$index
})
table(health[ref.index.p3[[1]],]$cell_type_3)

index.p3 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], ref.index.p3[[i]])
})

# P4 = 1:1
ref.index.p4 <- lapply(1:5, function(i){
  p=1
  endo <- health[health$cell_type_3=="Endocrine",]
  nonendo <- health[health$cell_type_3=="Non-Endocrine",]
  endo.ref <- endo %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(441+5)/nrow(endo))
  nonendo.ref <- nonendo %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(441+5)/nrow(nonendo))
  ref <- rbind(endo.ref, nonendo.ref)
  ref$index
})
table(health[ref.index.p4[[1]],]$cell_type_3)

index.p4 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], ref.index.p4[[i]])
})

# P5 = 0.5:1
ref.index.p5 <- lapply(1:5, function(i){
  p=0.5
  endo <- health[health$cell_type_3=="Endocrine",]
  nonendo <- health[health$cell_type_3=="Non-Endocrine",]
  endo.ref <- endo %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(441+5)/nrow(endo))
  nonendo.ref <- nonendo %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(441+5)/nrow(nonendo))
  ref <- rbind(endo.ref, nonendo.ref)
  ref$index
})
table(health[ref.index.p5[[1]],]$cell_type_3)

index.p5 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], ref.index.p5[[i]])
})

# P6 = 0.2:1
ref.index.p6 <- lapply(1:5, function(i){
  p=0.2
  endo <- health[health$cell_type_3=="Endocrine",]
  nonendo <- health[health$cell_type_3=="Non-Endocrine",]
  endo.ref <- endo %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(441+5)/nrow(endo))
  nonendo.ref <- nonendo %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(441+5)/nrow(nonendo))
  ref <- rbind(endo.ref, nonendo.ref)
  ref$index
})
table(health[ref.index.p6[[1]],]$cell_type_3)

index.p6 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], ref.index.p6[[i]])
})


# P7 = 0.1:1
ref.index.p7 <- lapply(1:5, function(i){
  p=0.1
  endo <- health[health$cell_type_3=="Endocrine",]
  nonendo <- health[health$cell_type_3=="Non-Endocrine",]
  endo.ref <- endo %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(441+5)/nrow(endo))
  nonendo.ref <- nonendo %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(441+5)/nrow(nonendo))
  ref <- rbind(endo.ref, nonendo.ref)
  ref$index
})
table(health[ref.index.p7[[1]],]$cell_type_3)

index.p7 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], ref.index.p7[[i]])
})

# Session 2: query count not equal to reference count, all reference cells invited, 3128 cells each iteration
# P8 = 10:1
ref.index.p8 <- lapply(1:5, function(i){
  p=10
  endo <- health[health$cell_type_3=="Endocrine",]
  nonendo <- health[health$cell_type_3=="Non-Endocrine",]
  endo.ref <- endo %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(3128+7)/nrow(endo))
  nonendo.ref <- nonendo %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(3128+7)/nrow(nonendo))
  ref <- rbind(endo.ref, nonendo.ref)
  ref$index
})
table(health[ref.index.p8[[1]],]$cell_type_3)

index.p8 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], ref.index.p8[[i]])
})

# P9 = 5:1
ref.index.p9 <- lapply(1:5, function(i){
  p=5
  endo <- health[health$cell_type_3=="Endocrine",]
  nonendo <- health[health$cell_type_3=="Non-Endocrine",]
  endo.ref <- endo %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(3128+7)/nrow(endo))
  nonendo.ref <- nonendo %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(3128+7)/nrow(nonendo))
  ref <- rbind(endo.ref, nonendo.ref)
  ref$index
})
table(health[ref.index.p9[[1]],]$cell_type_3)

index.p9 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], ref.index.p9[[i]])
})

# P10 = 2:1
ref.index.p10 <- lapply(1:5, function(i){
  p=2
  endo <- health[health$cell_type_3=="Endocrine",]
  nonendo <- health[health$cell_type_3=="Non-Endocrine",]
  endo.ref <- endo %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(3128+7)/nrow(endo))
  nonendo.ref <- nonendo %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(3128+7)/nrow(nonendo))
  ref <- rbind(endo.ref, nonendo.ref)
  ref$index
})
table(health[ref.index.p10[[1]],]$cell_type_3)

index.p10 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], ref.index.p10[[i]])
})

# P11 = 1:1
ref.index.p11 <- lapply(1:5, function(i){
  p=1
  endo <- health[health$cell_type_3=="Endocrine",]
  nonendo <- health[health$cell_type_3=="Non-Endocrine",]
  endo.ref <- endo %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(3128+7)/nrow(endo))
  nonendo.ref <- nonendo %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(3128+7)/nrow(nonendo))
  ref <- rbind(endo.ref, nonendo.ref)
  ref$index
})
table(health[ref.index.p11[[1]],]$cell_type_3)

index.p11 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], ref.index.p11[[i]])
})

# P12 = 0.5:1
ref.index.p12 <- lapply(1:5, function(i){
  p=0.5
  endo <- health[health$cell_type_3=="Endocrine",]
  nonendo <- health[health$cell_type_3=="Non-Endocrine",]
  endo.ref <- endo %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(3128+7)/nrow(endo))
  nonendo.ref <- nonendo %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(3128+7)/nrow(nonendo))
  ref <- rbind(endo.ref, nonendo.ref)
  ref$index
})
table(health[ref.index.p12[[1]],]$cell_type_3)

index.p12 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], ref.index.p12[[i]])
})

# P13 = 0.2:1
ref.index.p13 <- lapply(1:5, function(i){
  p=0.2
  endo <- health[health$cell_type_3=="Endocrine",]
  nonendo <- health[health$cell_type_3=="Non-Endocrine",]
  endo.ref <- endo %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(3128+7)/nrow(endo))
  nonendo.ref <- nonendo %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(3128+7)/nrow(nonendo))
  ref <- rbind(endo.ref, nonendo.ref)
  ref$index
})
table(health[ref.index.p13[[1]],]$cell_type_3)

index.p13 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], ref.index.p13[[i]])
})

# P14 = 0.1:1
ref.index.p14 <- lapply(1:5, function(i){
  p=0.1
  endo <- health[health$cell_type_3=="Endocrine",]
  nonendo <- health[health$cell_type_3=="Non-Endocrine",]
  endo.ref <- endo %>% group_by(cell_type_1) %>% slice_sample(prop=p/(1+p)*(3128+7)/nrow(endo))
  nonendo.ref <- nonendo %>% group_by(cell_type_1) %>% slice_sample(prop=1/(1+p)*(3128+7)/nrow(nonendo))
  ref <- rbind(endo.ref, nonendo.ref)
  ref$index
})
table(health[ref.index.p14[[1]],]$cell_type_3)

index.p14 <- lapply(1:5, function(i){
  list(query.index.s1[[i]], ref.index.p14[[i]])
})


for (i in 1:14) {
  saveRDS(get(paste0("index.p", i)), file = paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/01_datasets/index/Index_P", i, ".rds"))
}





