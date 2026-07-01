setwd("R:/ANNCOM2023-Q6025")
library(Seurat)
library(dplyr)
library(tidyverse)

pbmc.68k <- readRDS(file = "08_StartOver_NewStory/05_interdatasets/01_pbmc/01_datasets/seuobj/pbmc_68k.rds")
pbmc.covid <- readRDS(file = "08_StartOver_NewStory/05_interdatasets/01_pbmc/01_datasets/seuobj/pbmc_covid.rds")

# scmap-cluster #####
library(scmap)
library(SingleCellExperiment)
benchmarking_pbmc_scmapcluster <- function(query, reference, prop){
  set.seed(233)
  index <- readRDS(paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/01_datasets/index/Index_", prop,".rds"))
  running_time_df <- data.frame(Iteration = 1:length(index), Time  = NA)
  Ann.level = "cell_type_1"
  
  for (i in 1:length(index)) {
    gc()
    
    testing_index <- index[[i]][[1]]
    ref_index <- index[[i]][[2]]
    trial <- query[,testing_index]
    ref <- reference[,ref_index]
    
    trial  <- SCTransform(trial, verbose = TRUE, vst.flavor="v2", return.only.var.genes = FALSE)
    ref <- SCTransform(ref, verbose = TRUE, vst.flavor="v2", return.only.var.genes = FALSE)
    commonGenes = intersect(rownames(trial), rownames(ref))
    
    trial <- trial[commonGenes,]
    ref <- ref[commonGenes,]
    
    start_time <- Sys.time()
    
    trial <- DietSeurat(trial, assays = "SCT")
    ref <- DietSeurat(ref, assays = "SCT")
    
    trial <- as.SingleCellExperiment(trial, assay = "SCT")
    ref <- as.SingleCellExperiment(ref, assay = "SCT")
    
    rowData(ref)$feature_symbol <- rownames(ref)
    colData(ref)$cell_type1 <- ref$cell_type_1
    rowData(trial)$feature_symbol <- rownames(trial)
    
    ref <- selectFeatures(ref, suppress_plot = FALSE)
    
    table(rowData(ref)$scmap_features)
    
    ref <- indexCluster(ref)
    
    scmapCluster_results <- scmapCluster(
      projection = trial, 
      index_list = list(
        ref = metadata(ref)$scmap_cluster_index
      ))
    
    pred_results <- data.frame(Barcodes = colnames(trial), True_label = colData(trial)[,Ann.level], Prediction = as.character(scmapCluster_results$scmap_cluster_labs))
    
    result.saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/scmapcluster/", prop,"/prediction_results_", i, "_5.csv")
    
    write.csv(pred_results, file = result.saving.path)
    
    end_time <- Sys.time()
    running_time <- difftime(end_time, start_time, units = 'mins')
    
    running_time_df[i,"Time"] <- as.vector(running_time)
    
    rm(ref, trial, scmapCluster_results, end_time, start_time, running_time, pred_results)
    print(paste0("done_", i))
  }
  write.csv(running_time_df, paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/scmapcluster/", prop,"/running_time.csv"))
  print("Done!")
  return()
}
benchmarking_pbmc_scmapcluster(query = pbmc.covid, reference = pbmc.68k, prop="P1")

# singleCellNet #####
library(singleCellNet)
benchmarking_pbmc_singleCellNet <- function(query, reference, prop){
  set.seed(233)
  index <- readRDS(paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/01_datasets/index/Index_", prop,".rds"))
  running_time_df <- data.frame(Iteration = 1:length(index), Time  = NA)
  Ann.level = "cell_type_1"
  
  for (i in 1:length(index)) {
    gc()
    
    testing_index <- index[[i]][[1]]
    ref_index <- index[[i]][[2]]
    trial <- query[,testing_index]
    ref <- reference[,ref_index]
    
    trial  <- SCTransform(trial, verbose = TRUE, vst.flavor="v2", return.only.var.genes = FALSE)
    ref <- SCTransform(ref, verbose = TRUE, vst.flavor="v2", return.only.var.genes = FALSE)
    commonGenes = intersect(rownames(trial), rownames(ref))
    
    trial <- trial[commonGenes,]
    ref <- ref[commonGenes,]
    
    start_time <- Sys.time()
    
    trial <- DietSeurat(trial, assays = "SCT")
    ref <- DietSeurat(ref, assays = "SCT")
    
    ref@meta.data$cell <- colnames(ref)
    seuratfile <- extractSeurat(ref, exp_slot_name = "data")
    
    sampTab <- seuratfile$sampTab
    expDat <- seuratfile$expDat
    rm(seuratfile)
    
    super_rare = names(which(table(sampTab$cell_type_1)==1))
    
    if(length(super_rare) != 0){
      ref <- subset(ref, subset = cell_type_1 %in% setdiff(unique(ref$cell_type_1), super_rare))
      seuratfile <- extractSeurat(ref, exp_slot_name = "data")
      sampTab <- seuratfile$sampTab
      expDat <- seuratfile$expDat
      rm(seuratfile)
    }
    
    if(min(table(sampTab$cell_type_1)) < 3) {
      cells_reserved = 2
    } else {
      cells_reserved = 3
    }
    
    stList <- splitCommon(sampTab=sampTab, dLevel=Ann.level, cells_reserved = cells_reserved)
    stTrain <- stList[[1]]
    expTrain <- expDat[,rownames(stTrain)]
    class_info <- scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 100, nRand = 70, nTrees = 1000, nTopGenePairs = 250, dLevel = Ann.level, colName_samp = "cell")
    stTestList = splitCommon(sampTab=stList[[2]], dLevel=Ann.level, cells_reserved = cells_reserved) 
    stTest = stTestList[[1]]
    expTest = expDat[,rownames(stTest)]
    classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 50)
    trial@meta.data$cell <- colnames(trial)
    seuratfile.trial <- extractSeurat(trial, exp_slot_name = "data")
    rm(trial)
    sampTab.trial <- seuratfile.trial$sampTab
    expDat.trial <- seuratfile.trial$expDat
    rm(seuratfile.trial)
    crParkall <- scn_predict(class_info[['cnProc']], expDat.trial, nrand=50)
    sampTab.trial <- get_cate(classRes = crParkall, sampTab = sampTab.trial, dLevel = Ann.level, sid = "cell", nrand = 50)
    
    pred_results <- data.frame(Barcodes = sampTab.trial[, "cell"], True_label = sampTab.trial[,Ann.level], Prediction = sampTab.trial[,"category"])
    
    result.saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/SingleCellNet/", prop,"/prediction_results_", i, "_5.csv")
    write.csv(pred_results, file = result.saving.path)
    end_time <- Sys.time()
    running_time <- difftime(end_time, start_time, units = 'mins')
    running_time_df[i,"Time"] <- as.vector(running_time)
    
    rm(start_time, testing_index, sampTab, expDat, stList, stTrain, expTrain, class_info, stTestList, stTest, expTest, classRes_val_all, sampTab.trial, expDat.trial, crParkall, running_time, end_time, result.saving.path)
    print(paste0("done_", i))
  }
  write.csv(running_time_df, paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/SingleCellNet/",prop,"/running_time.csv"))
  print("Done!")
  return()
}
benchmarking_pbmc_singleCellNet(query = pbmc.covid, reference = pbmc.68k, prop="P1")

# CHETAH #####
library(CHETAH)
benchmarking_pbmc_CHETAH <- function(query, reference, prop){
  set.seed(233)
  index <- readRDS(paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/01_datasets/index/Index_", prop,".rds"))
  running_time_df <- data.frame(Iteration = 1:length(index), Time  = NA)
  ribo <- read.table("02_labelling/CHETAH/ribosomal.txt", header = FALSE, sep = '\t')
  Ann.level = "cell_type_1"
  
  for (i in 1:length(index)) {
    gc()
    
    testing_index <- index[[i]][[1]]
    ref_index <- index[[i]][[2]]
    trial <- query[,testing_index]
    ref <- reference[,ref_index]
    
    trial  <- SCTransform(trial, verbose = TRUE, vst.flavor="v2", return.only.var.genes = FALSE)
    ref <- SCTransform(ref, verbose = TRUE, vst.flavor="v2", return.only.var.genes = FALSE)
    commonGenes = intersect(rownames(trial), rownames(ref))
    
    trial <- trial[commonGenes,]
    ref <- ref[commonGenes,]
    
    start_time <- Sys.time()
    
    trial <- DietSeurat(trial, assays = "SCT")
    ref <- DietSeurat(ref, assays = "SCT")
    
    trial <- as.SingleCellExperiment(trial, assay = "SCT")
    ref <- as.SingleCellExperiment(ref, assay = "SCT")
    
    ref <- ref[!rownames(ref) %in% ribo[,1], ]
    trial <- CHETAHclassifier(input = trial,
                              ref_cells = ref,
                              ref_ct = Ann.level)
    
    pred_results <- data.frame(Barcodes=rownames(trial@colData), True_label=trial@colData[,Ann.level], Prediction_ori=trial@colData$celltype_CHETAH)
    
    pred_results <- pred_results %>% mutate(Prediction = case_when(Prediction_ori %in% unique(ref@colData[,Ann.level]) ~ Prediction_ori,
                                                                   .default = "unassinged"))
    
    result.saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/CHETAH/", prop,"/prediction_results_", i, "_", length(index), ".csv")
    write.csv(pred_results, file = result.saving.path)
    end_time <- Sys.time()
    running_time <- difftime(end_time, start_time, units = 'mins')
    running_time_df[i,"Time"] <- as.vector(running_time)
    
    rm(start_time, testing_index, running_time, end_time, result.saving.path, trial, ref, pred_results)
    print(paste0("done_", i))
  }
  write.csv(running_time_df, paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/CHETAH/", prop, "/running_time.csv"))
  print("Done!")
  return()
}
benchmarking_pbmc_CHETAH(query = pbmc.covid, reference = pbmc.68k, prop="P1")

# SingleR #####
library(SingleR)
library(SingleCellExperiment)
benchmarking_pbmc_SingleR <- function(query, reference, prop){
  set.seed(233)
  index <- readRDS(paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/01_datasets/index/Index_", prop,".rds"))
  running_time_df <- data.frame(Iteration = 1:length(index), Time  = NA)
  Ann.level = "cell_type_1"
  
  for (i in 1:length(index)) {
    gc()
    
    
    testing_index <- index[[i]][[1]]
    ref_index <- index[[i]][[2]]
    trial <- query[,testing_index]
    ref <- reference[,ref_index]
    
    trial  <- SCTransform(trial, verbose = TRUE, vst.flavor="v2", return.only.var.genes = FALSE)
    ref <- SCTransform(ref, verbose = TRUE, vst.flavor="v2", return.only.var.genes = FALSE)
    commonGenes = intersect(rownames(trial), rownames(ref))
    
    trial <- trial[commonGenes,]
    ref <- ref[commonGenes,]
    
    start_time <- Sys.time()
    
    trial <- DietSeurat(trial, assays = "SCT")
    ref <- DietSeurat(ref, assays = "SCT")
    
    ref <- as.SingleCellExperiment(ref, assay = "SCT")
    DefaultAssay(trial) <- "SCT"
    sce <- as.SingleCellExperiment(DietSeurat(trial))
    singler.annotation <- SingleR(test = sce, assay.type.test = 1, ref = ref, labels = ref@colData[, Ann.level])
    trial@meta.data$prediction <- singler.annotation$pruned.labels
    pred_results <- data.frame(Barcodes = rownames(trial@meta.data), True_label = trial@meta.data[, Ann.level], Prediction = trial@meta.data[,"prediction"])
    
    result.saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/SingleR/", prop,"/prediction_results_", i, "_", length(index), ".csv")
    
    write.csv(pred_results, file = result.saving.path)
    
    end_time <- Sys.time()
    running_time <- difftime(end_time, start_time, units = 'mins')
    
    running_time_df[i,"Time"] <- as.vector(running_time)
    
    rm(start_time, testing_index, trial, ref, sce, singler.annotation, pred_results,result.saving.path, end_time, running_time)
    print(paste0("done_", i))
  }
  write.csv(running_time_df, paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/SingleR/", prop, "/running_time.csv"))
  print("Done!")
  return()
}
benchmarking_pbmc_SingleR(query = pbmc.covid, reference = pbmc.68k, prop="P1")

# scClassify #######
library(Seurat)
library(scClassify)
benchmarking_pbmc_scClassify <- function(query, reference, prop){
  index <- readRDS(paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/01_datasets/index/Index_", prop,".rds"))
  running_time_df <- data.frame(Iteration = 1:length(index), Time = NA)
  for (ind in 1:length(index)) {
    gc()
    
    testing_index <- index[[ind]][[1]]
    ref_index <- index[[ind]][[2]]
    trial <- query[,testing_index]
    ref <- reference[,ref_index]
    
    trial  <- SCTransform(trial, verbose = TRUE, vst.flavor="v2", return.only.var.genes = FALSE)
    ref <- SCTransform(ref, verbose = TRUE, vst.flavor="v2", return.only.var.genes = FALSE)
    commonGenes = intersect(rownames(trial), rownames(ref))
    
    trial <- trial[commonGenes,]
    ref <- ref[commonGenes,]
    
    start_time <- Sys.time()
    
    trial <- DietSeurat(trial, assays = "SCT")
    ref <- DietSeurat(ref, assays = "SCT")
    
    super_rare = names(which(table(ref$cell_type_1)==1))
    
    if(length(super_rare) != 0){
      ref <- subset(ref, subset = cell_type_1 %in% setdiff(unique(ref$cell_type_1), super_rare))
    }
    
    scClassify_res <- scClassify(exprsMat_train = ref@assays$SCT@data,
                                 cellTypes_train = ref$cell_type_1,
                                 exprsMat_test = list(test = trial@assays$SCT@data),
                                 tree = "HOPACH",
                                 algorithm = "WKNN",
                                 selectFeatures = c("limma"),
                                 similarity = c("pearson"),
                                 returnList = FALSE,
                                 verbose = TRUE)
    
    pred_results <- data.frame(Barcodes = colnames(trial), True_label = trial$cell_type_1, Prediction = scClassify_res$testRes$test$pearson_WKNN_limma$predRes)
    result.saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/scClassify/", prop,"/prediction_results_", ind, "_", length(index), ".csv")
    write.csv(pred_results, file = result.saving.path)
    end_time <- Sys.time()
    running_time <- difftime(end_time, start_time, units = 'mins')
    running_time_df[ind,"Time"] <- as.vector(running_time)
    print(paste0("done_", ind))
    rm(testing_index, ref, trial, scClassify_res, start_time, result.saving.path, end_time, running_time)
  }
  
  write.csv(running_time_df, paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/scClassify/", prop,"/running_time.csv"))
  print("Done!")
  return()
}
benchmarking_pbmc_scClassify(query = pbmc.covid, reference = pbmc.68k, prop="P1")

# scPred ########
library(scPred)
library(Seurat)
library(magrittr)
library(harmony)
library(SingleCellExperiment)
library(tidyverse)
library(dplyr)
benchmarking_pbmc_scPred <- function(query, reference, prop){
  index <- readRDS(paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/01_datasets/index/Index_", prop,".rds"))
  running_time_df <- data.frame(Iteration = 1:length(index), Time = NA)
  for (ind in 2:5) {
    gc()
    
    testing_index <- index[[ind]][[1]]
    ref_index <- index[[ind]][[2]]
    trial <- query[,testing_index]
    ref <- reference[,ref_index]
    
    trial  <- SCTransform(trial, verbose = TRUE, vst.flavor="v2", return.only.var.genes = FALSE)
    ref <- SCTransform(ref, verbose = TRUE, vst.flavor="v2", return.only.var.genes = FALSE)
    commonGenes = intersect(rownames(trial), rownames(ref))
    
    trial <- trial[commonGenes,]
    ref <- ref[commonGenes,]
    
    start_time <- Sys.time()
    
    # trial <- DietSeurat(trial, assays = "SCT")
    # ref <- DietSeurat(ref, assays = "SCT")
    
    super_rare = names(which(table(ref$cell_type_1)==1))
    
    if(length(super_rare) != 0){
      ref <- subset(ref, subset = cell_type_1 %in% setdiff(unique(ref$cell_type_1), super_rare))
    }
    
    super_rare = names(which(table(ref$cell_type_1)==2))
    
    if(length(super_rare) != 0){
      ref <- subset(ref, subset = cell_type_1 %in% setdiff(unique(ref$cell_type_1), super_rare))
    }
    
    ref <- RunPCA(ref, features = VariableFeatures(ref), assay = "SCT")
    ref <- RunUMAP(ref, dims = 1:30, reduction = "pca")
    ref <- getFeatureSpace(ref, "cell_type_1")
    
    ref <- trainModel(ref) #SVM is the default
    
    trial <- scPredict(trial, ref)
    
    trial <- RunUMAP(trial, reduction = "scpred", dims = 1:30)
    
    pred_results <- data.frame(Barcodes = colnames(trial), True_label = trial$cell_type_1, Prediction = trial$scpred_prediction)
    result.saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/scPred/", prop,"/prediction_results_", ind, "_", length(index), ".csv")
    write.csv(pred_results, file = result.saving.path)
    end_time <- Sys.time()
    running_time <- difftime(end_time, start_time, units = 'mins')
    running_time_df[ind,"Time"] <- as.vector(running_time)
    print(paste0("done_", ind))
    rm(testing_index, ref, trial, scClassify_res, start_time, result.saving.path, end_time, running_time)
  }
  
  write.csv(running_time_df, paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/scPred/", prop,"/running_time.csv"))
  print("Done!")
  return()
}
benchmarking_pbmc_scPred(query = pbmc.covid, reference = pbmc.68k, prop="P2")
# hpc is fine with 1,2,3 and my commputer is fine with 2, 3, 4, 5
# iteration = 4, Time difference of 5.217075 mins
# iteration = 5, Time difference of 5.31078 mins
# iteration = 1, Time difference of 5.259438 mins (hpc)