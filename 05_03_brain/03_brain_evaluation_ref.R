setwd("R:/ANNCOM2023-Q6025")
library(Seurat)
library(tidyverse)
library(dplyr)
library(caret)

method.ref <- c("CHETAH", "scmapcluster", "SingleCellNet", "SingleR", "scClassify", "scPred")

evaluation_brain_function_ref_ct1 <- function(method, prop){
  index <- readRDS(file = paste0("08_StartOver_NewStory/05_interdatasets/03_brain/01_datasets/index/Index_",prop,".rds"))

  q_cell_types <- c("Pax6", "L5/6_NP", "IT", "L6_CT", "L4_IT", "Astrocyte",      
                    "Vip", "Sst_Chodl", "L6_IT_Car3", "Sncg", "Pvalb", "Oligodendrocyte",
                    "Sst", "VLMC", "Lamp5", "Microglia", "OPC", "L6b",            
                    "Endothelial", "L5_ET", "Chandelier")
  r_cell_types <- c("Vip", "Lamp5", "IT", "Pax6", "Oligodendrocyte", "Astrocyte",      
                    "L5/6_IT_Car3", "L5/6_NP", "Sst", "L6_CT", "OPC", "Pvalb",          
                    "L6b", "Microglia", "L5_ET", "Pericyte", "Endothelial", "L4_IT", "VLMC")
  
  shared_cell_types <- intersect(q_cell_types, r_cell_types)
  q_only_cell_types <- setdiff(q_cell_types, shared_cell_types)
  r_only_cell_types <- setdiff(r_cell_types, shared_cell_types)
  union_cell_types <- unique(c(q_cell_types, r_cell_types))
  
  for (i in 1:length(method)){
    k_index = length(index)
    eval.tbl <- data.frame(Method=rep(method[i], k_index), Iteration=1:k_index, KFold = k_index, 
                           assignedCounts=NA, assignedRate=NA,unassignedCounts=NA, unassignedRate=NA, 
                           Accuracy=NA, WeightedAveF1=NA, ARI=NA)
    
    eval.tbl$CellCounts <- sapply(1:5, function(ind){length(index[[ind]][[1]])})
    
    running.time <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/03_brain/02_labelling/", method[i],"/", prop,"/running_time.csv"), row.names = 1)
    
    res.ls <- vector("list", k_index)
    
    for (j in 1:k_index) {
      res <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/03_brain/02_labelling/", method[i],"/", prop,"/prediction_results_",j, "_5.csv"), row.names=1)
      names(res.ls)[j] <- paste0("res",j)
      res.ls[[j]] <- res
      rm(res)
    }
    
    res.ls <- lapply(res.ls, function(res){
      rownames(res) <- res[,1] 
      res <- res[,-1]
      res <- res[,c("True_label", "Prediction")]
    })
    
    res.common.ls <- lapply(res.ls, function(res){
      res_notna <- res %>% drop_na(Prediction)
      res_common <- res_notna[res_notna$True_label %in% shared_cell_types,]
      res_common <- res_common[res_common$Prediction %in% r_cell_types,]
      res_common$True_label <- factor(res_common$True_label, levels = shared_cell_types)
      res_common$Prediction <- factor(res_common$Prediction, levels = c(shared_cell_types, r_only_cell_types))
      return(res_common)
    })
    
    res.qonly.ls <- lapply(res.ls, function(res){
      res_notna <- res %>% replace_na(list(Prediction="missing"))
      res_qonly <- res_notna[res_notna$True_label %in% q_only_cell_types,]
      return(res_qonly)
    })
    
    res.id.ls <- lapply(res.ls, function(res){
      res_notna <- res %>% drop_na(Prediction)
      res_identified <- res_notna[res_notna$Prediction %in% r_cell_types,]
      res_identified$True_label <- factor(res_identified$True_label, levels = union_cell_types)
      res_identified$Prediction <- factor(res_identified$Prediction, levels = union_cell_types)
      return(res_identified)
    })
    
    eval.tbl$assignedCounts <- sapply(res.id.ls, nrow)
    eval.tbl$unassignedCounts <- eval.tbl$CellCounts - eval.tbl$assignedCounts
    eval.tbl$assignedRate <- eval.tbl$assignedCounts/eval.tbl$CellCounts *100
    eval.tbl$unassignedRate <- eval.tbl$unassignedCounts/eval.tbl$CellCounts *100
    
    cm.ls <- lapply(res.common.ls, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
    })
    
    cm.ls <- cm.ls[!sapply(cm.ls,is.null)]
    
    maxlength <- length(index)
    
    acc <- sapply(cm.ls, function(cm){
      n = sum(cm)
      diag = diag(cm)
      accuracy = sum(diag) / n 
    })
    
    acc = c(acc, rep(NA, maxlength - length(acc)))  
    eval.tbl$Accuracy <- acc
    
    wa.f1 <- sapply(cm.ls, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      
      precision = diag / colsums[rownames(cm)]
      recall = diag / rowsums 
      f1 = 2 * precision * recall / (precision + recall) 
      
      f1[is.na(f1)] <- 0
      WA.f1 <- weighted.mean(f1,p)
    })
    wa.f1 = c(wa.f1, rep(NA, maxlength - length(wa.f1)))
    eval.tbl$WeightedAveF1 <- wa.f1
    
    precision <- sapply(cm.ls, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      precision = diag / colsums[rownames(cm)]
      precision[is.na(precision)] <- 0
      precision = weighted.mean(precision, p)
    })
    eval.tbl$Precision <- precision
    
    recall <- sapply(cm.ls, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      recall = diag / rowsums
      recall[is.na(recall)] <- 0
      recall = weighted.mean(recall,p)
    })
    eval.tbl$Recall <- recall
    
    
    rare_cell_type = c("Astrocyte", "Chandelier", "Endothelial", "L5/6_NP", "L5_ET", "L6_CT", "L6_IT_Car3",     
                         "L6b", "Lamp5", "Microglia", "Oligodendrocyte", "OPC", "Pax6", "Sncg",           
                         "Sst_Chodl", "VLMC")
    abundant_cell_type = c("IT", "L4_IT", "Pvalb", "Sst", "Vip")
    
    res.rare <- lapply(res.common.ls, function(df){
      df.temp <- df[df$True_label %in% intersect(shared_cell_types, rare_cell_type),]
      return(df.temp)
    })
    
    res.rich <- lapply(res.common.ls, function(df){
      df.temp <- df[df$True_label %in% intersect(shared_cell_types, abundant_cell_type),]
      return(df.temp)
    })
    
    cm.ls.rare <- lapply(res.rare, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
    })
    cm.ls.rich <- lapply(res.rich, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
    })
    
    
    acc.rare <- sapply(cm.ls.rare, function(cm){
      n = sum(cm)
      diag = diag(cm)
      accuracy = sum(diag) / n 
    })
    eval.tbl$Accuracy.rare <- acc.rare
    
    acc.rich <- sapply(cm.ls.rich, function(cm){
      n = sum(cm)
      diag = diag(cm)
      accuracy = sum(diag) / n 
    })
    eval.tbl$Accuracy.rich <- acc.rich
    
    wa.f1.rare <- sapply(cm.ls.rare, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      
      precision = diag / colsums[rownames(cm)] 
      recall = diag / rowsums 
      f1 = 2 * precision * recall / (precision + recall) 
      
      f1[is.na(f1)] <- 0
      WA.f1 <- weighted.mean(f1,p)
    })
    eval.tbl$WeightedAveF1.rare <- wa.f1.rare
    
    wa.f1.rich <- sapply(cm.ls.rich, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      
      precision = diag / colsums[rownames(cm)] 
      recall = diag / rowsums 
      f1 = 2 * precision * recall / (precision + recall) 
      
      f1[is.na(f1)] <- 0
      WA.f1 <- weighted.mean(f1,p)
    })
    eval.tbl$WeightedAveF1.rich <- wa.f1.rich
    
    ari <- sapply(res.common.ls, function(res_identified){
      ARI <- pdfCluster::adj.rand.index(res_identified$True_label, res_identified$Prediction)
      return(ARI)
    })
    eval.tbl$ARI <- ari
    
    vscore <- lapply(res.common.ls, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
      
      cm.types <- rownames(cm)
      
      res.num <- data.frame(True_label=rep(NA, nrow(res_identified)), Prediction=rep(NA, nrow(res_identified)))
      res.num$True_label <- c(1:length(cm.types), res_identified$True_label)[match(res_identified$True_label, c(cm.types, res_identified$True_label))]
      res.num$Prediction <- c(1:length(cm.types), res_identified$Prediction)[match(res_identified$Prediction, c(cm.types, res_identified$Prediction))]
      v.score <- sabre::vmeasure(res.num$True_label, res.num$Prediction)
      return(unlist(v.score))
    })
    eval.tbl <- cbind(eval.tbl, data.frame(t(sapply(vscore,c))))
    
    honesty <- sapply(res.qonly.ls, function(res){
      honesty <- nrow(res[res$Prediction%in%c("unassigned", "unassinged", "missing"),])/nrow(res) * 100
    })
    
    eval.tbl$honesty <- honesty
    
    eval.tbl <- cbind(eval.tbl, running.time$Time)
    names(eval.tbl)[ncol(eval.tbl)] <- "Time"
    
    saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/03_brain/03_evaluation/", method[i], "/", prop,"/evaluating_result_K5.csv")
    write.csv(eval.tbl, file = saving.path)
  }
}
evaluation_brain_function_ref_ct2 <- function(method, prop){
  
  index <- readRDS(file = paste0("08_StartOver_NewStory/05_interdatasets/03_brain/01_datasets/index/Index_",prop,".rds"))
  
  q_cell_types <- c("Pax6", "L5/6_NP", "IT", "L6_CT", "Astrocyte", "Vip",            
                    "Sst", "Sncg", "Pvalb", "Oligodendrocyte", "VLMC", "Lamp5",          
                    "Microglia", "OPC", "L6b", "Endothelial", "L5_ET", "Chandelier")
  r_cell_types <- c("Vip", "Lamp5", "IT", "Pax6", "Oligodendrocyte", "Astrocyte",      
                    "L5/6_NP", "Sst", "L6_CT", "OPC", "Pvalb", "L6b",            
                    "Microglia", "L5_ET", "Pericyte", "Endothelial", "VLMC")
  
  shared_cell_types <- intersect(q_cell_types, r_cell_types)
  q_only_cell_types <- setdiff(q_cell_types, shared_cell_types)
  r_only_cell_types <- setdiff(r_cell_types, shared_cell_types)
  union_cell_types <- unique(c(q_cell_types, r_cell_types))
  
  for (i in 1:length(method)) {
    k_index = length(index)
    eval.tbl <- data.frame(Method=rep(method[i], k_index), Iteration=1:k_index, KFold = k_index, 
                           assignedCounts=NA, assignedRate=NA,unassignedCounts=NA, unassignedRate=NA, 
                           Accuracy=NA, WeightedAveF1=NA, ARI=NA)
    
    eval.tbl$CellCounts <- sapply(1:5, function(ind){length(index[[ind]][[1]])})
    
    running.time <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/03_brain/02_labelling/", method[i],"/", prop,"/running_time.csv"), row.names = 1)
    
    res.ls <- vector("list", length(index))
    
    for (j in 1:length(index)) {
      res <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/03_brain/02_labelling/", method[i],"/", prop,"/prediction_results_",j, "_", length(index),".csv"), row.names=1)
      names(res.ls)[j] <- paste0("res",j)
      res.ls[[j]] <- res
      rm(res)
    }
    
    # res.ls <- res.ls[!sapply(res.ls,is.null)]
    
    res.ls <- lapply(res.ls, function(res){
      rownames(res) <- res[,1] 
      res <- res[,-1]
      res <- res[,c("True_label", "Prediction")]
    })
    
    res.ls <- lapply(res.ls, function(res){
      res <- res %>% mutate(True_label = case_when(True_label %in% c("L4_IT", "L6_IT_Car3", "IT") ~ "IT",
                                                   True_label %in% c("Sst", "Sst_Chodl") ~ "Sst",
                                                   .default = True_label))
      
      res <- res %>% mutate(Prediction = case_when(Prediction %in% c("L4_IT", "IT", "L5/6_IT_Car3") ~ "IT",
                                                   .default = Prediction))
    })
    
    res.common.ls <- lapply(res.ls, function(res){
      res_notna <- res %>% drop_na(Prediction)
      res_common <- res_notna[res_notna$True_label %in% shared_cell_types,]
      res_common <- res_common[res_common$Prediction %in% r_cell_types,]
      res_common$True_label <- factor(res_common$True_label, levels = shared_cell_types)
      res_common$Prediction <- factor(res_common$Prediction, levels = c(shared_cell_types, r_only_cell_types))
      return(res_common)
    })
    
    res.qonly.ls <- lapply(res.ls, function(res){
      res_notna <- res %>% replace_na(list(Prediction="missing"))
      res_qonly <- res_notna[res_notna$True_label %in% q_only_cell_types,]
      return(res_qonly)
    })
    
    res.id.ls <- lapply(res.ls, function(res){
      res_notna <- res %>% drop_na(Prediction)
      res_identified <- res_notna[res_notna$Prediction %in% r_cell_types,]
      res_identified$True_label <- factor(res_identified$True_label, levels = union_cell_types)
      res_identified$Prediction <- factor(res_identified$Prediction, levels = union_cell_types)
      return(res_identified)
    })
    
    eval.tbl$assignedCounts <- sapply(res.id.ls, nrow)
    eval.tbl$unassignedCounts <- eval.tbl$CellCounts - eval.tbl$assignedCounts
    eval.tbl$assignedRate <- eval.tbl$assignedCounts/eval.tbl$CellCounts *100
    eval.tbl$unassignedRate <- eval.tbl$unassignedCounts/eval.tbl$CellCounts *100
    
    cm.ls <- lapply(res.common.ls, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
    })
    maxlength = 5
    cm.ls <- cm.ls[!sapply(cm.ls,is.null)]
    acc <- sapply(cm.ls, function(cm){
      n = sum(cm)
      diag = diag(cm)
      accuracy = sum(diag) / n 
    })
    acc = c(acc, rep(NA, maxlength - length(acc)))
    eval.tbl$Accuracy <- acc
    
    wa.f1 <- sapply(cm.ls, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      
      precision = diag / colsums[rownames(cm)]
      recall = diag / rowsums 
      f1 = 2 * precision * recall / (precision + recall) 
      
      f1[is.na(f1)] <- 0
      WA.f1 <- weighted.mean(f1,p)
    })
    
    wa.f1 = c(wa.f1, rep(NA, maxlength - length(wa.f1)))
    eval.tbl$WeightedAveF1 <- wa.f1
    
    precision <- sapply(cm.ls, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      precision = diag / colsums[rownames(cm)]
      precision[is.na(precision)] <- 0
      precision = weighted.mean(precision, p)
    })
    eval.tbl$Precision <- precision
    
    recall <- sapply(cm.ls, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      recall = diag / rowsums
      recall[is.na(recall)] <- 0
      recall = weighted.mean(recall,p)
    })
    eval.tbl$Recall <- recall
  
    rare_cell_type = c("Astrocyte", "Chandelier", "Endothelial", "L5/6_NP", "L5_ET", "L6_CT", "L6b",            
                         "Lamp5", "Microglia", "Oligodendrocyte", "OPC", "Pax6", "Sncg", "VLMC")
    abundant_cell_type = c("IT", "Pvalb", "Sst", "Vip")
    
    res.rare <- lapply(res.common.ls, function(df){
      df.temp <- df[df$True_label %in% intersect(shared_cell_types, rare_cell_type),]
      return(df.temp)
    })
    res.rich <- lapply(res.common.ls, function(df){
      df.temp <- df[df$True_label %in% intersect(shared_cell_types, abundant_cell_type),]
      return(df.temp)
    })
    
    cm.ls.rare <- lapply(res.rare, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
    })
    
    
    cm.ls.rich <- lapply(res.rich, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
    })
    
    acc.rare <- sapply(cm.ls.rare, function(cm){
      n = sum(cm)
      diag = diag(cm)
      accuracy = sum(diag) / n 
    })
    eval.tbl$Accuracy.rare <- acc.rare
    
    acc.rich <- sapply(cm.ls.rich, function(cm){
      n = sum(cm)
      diag = diag(cm)
      accuracy = sum(diag) / n 
    })
    eval.tbl$Accuracy.rich <- acc.rich
    
    wa.f1.rare <- sapply(cm.ls.rare, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      
      precision = diag / colsums[rownames(cm)] 
      recall = diag / rowsums 
      f1 = 2 * precision * recall / (precision + recall) 
      
      f1[is.na(f1)] <- 0
      WA.f1 <- weighted.mean(f1,p)
    })
    eval.tbl$WeightedAveF1.rare <- wa.f1.rare
    
    wa.f1.rich <- sapply(cm.ls.rich, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      
      precision = diag / colsums[rownames(cm)] 
      recall = diag / rowsums 
      f1 = 2 * precision * recall / (precision + recall) 
      
      f1[is.na(f1)] <- 0
      WA.f1 <- weighted.mean(f1,p)
    })
    eval.tbl$WeightedAveF1.rich <- wa.f1.rich
    
    ari <- sapply(res.common.ls, function(res_identified){
      ARI <- pdfCluster::adj.rand.index(res_identified$True_label, res_identified$Prediction)
      return(ARI)
    })
    eval.tbl$ARI <- ari
    
    vscore <- lapply(res.common.ls, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
      
      cm.types <- rownames(cm)
      
      res.num <- data.frame(True_label=rep(NA, nrow(res_identified)), Prediction=rep(NA, nrow(res_identified)))
      res.num$True_label <- c(1:length(cm.types), res_identified$True_label)[match(res_identified$True_label, c(cm.types, res_identified$True_label))]
      res.num$Prediction <- c(1:length(cm.types), res_identified$Prediction)[match(res_identified$Prediction, c(cm.types, res_identified$Prediction))]
      v.score <- sabre::vmeasure(res.num$True_label, res.num$Prediction)
      return(unlist(v.score))
    })
    eval.tbl <- cbind(eval.tbl, data.frame(t(sapply(vscore,c))))
    
    honesty <- sapply(res.qonly.ls, function(res){
      honesty <- nrow(res[res$Prediction%in%c("unassigned", "unassinged", "missing"),])/nrow(res) * 100
    })
    
    eval.tbl$honesty <- honesty
    
    eval.tbl <- cbind(eval.tbl, running.time$Time)
    names(eval.tbl)[ncol(eval.tbl)] <- "Time"
    
    saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/03_brain/03_evaluation/", method[i], "/", prop,"/evaluating_result_K5_CT2.csv")
    write.csv(eval.tbl, file = saving.path)
  }
}
evaluation_brain_function_ref_ct3 <- function(method, prop){
  
  index <- readRDS(file = paste0("08_StartOver_NewStory/05_interdatasets/03_brain/01_datasets/index/Index_",prop,".rds"))
  
  q_cell_types <- c("GABAergic", "Glutamatergic", "Astrocyte", "Oligodendrocyte", "VLMC", "Microglia",      
                    "OPC", "Endothelial")
  r_cell_types <- c("GABAergic", "Glutamatergic", "Oligodendrocyte", "Astrocyte", "OPC", "Microglia",      
                    "Pericyte", "Endothelial", "VLMC")
  
  shared_cell_types <- intersect(q_cell_types, r_cell_types)
  q_only_cell_types <- setdiff(q_cell_types, shared_cell_types)
  r_only_cell_types <- setdiff(r_cell_types, shared_cell_types)
  union_cell_types <- unique(c(q_cell_types, r_cell_types))
  
  for (i in 1:length(method)) {
    k_index = length(index)
    eval.tbl <- data.frame(Method=rep(method[i], k_index), Iteration=1:k_index, KFold = length(index), 
                           assignedCounts=NA, assignedRate=NA,unassignedCounts=NA, unassignedRate=NA, 
                           Accuracy=NA, WeightedAveF1=NA, ARI=NA)
    
    eval.tbl$CellCounts <- sapply(1:5, function(ind){length(index[[ind]][[1]])})

    running.time <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/03_brain/02_labelling/", method[i],"/", prop,"/running_time.csv"), row.names = 1)
    
    res.ls <- vector("list", length(index))
    
    for (j in 1:length(index)) {
      res <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/03_brain/02_labelling/", method[i],"/", prop,"/prediction_results_",j, "_", length(index),".csv"), row.names=1)
      names(res.ls)[j] <- paste0("res",j)
      res.ls[[j]] <- res
      rm(res)
    }
    
    res.ls <- lapply(res.ls, function(res){
      rownames(res) <- res[,1] 
      res <- res[,-1]
      res <- res[,c("True_label", "Prediction")]
    })
    
    res.ls <- lapply(res.ls, function(res){
      res <- res %>% mutate(True_label = case_when(True_label %in% c("Chandelier", "Lamp5", "Pax6", "Pvalb", "Sncg", "Sst", "Sst_Chodl", "Vip") ~ "GABAergic",
                                                   True_label %in% c("IT", "L4_IT", "L5/6_NP", "L5_ET", "L6_CT", "L6_IT_Car3", "L6b") ~ "Glutamatergic",
                                                   True_label == "Oligodendrocyte" ~ "Oligodendrocyte",
                                                   True_label == "Astrocyte" ~ "Astrocyte",
                                                   True_label == "OPC" ~ "OPC",
                                                   True_label == "Microglia" ~ "Microglia",
                                                   True_label == "Endothelial" ~ "Endothelial",
                                                   True_label == "VLMC" ~ "VLMC",
                                                   .default = True_label))
      
      res <- res %>% mutate(Prediction = case_when(Prediction %in% c("Lamp5", "Pax6", "Pvalb", "Sst", "Vip") ~ "GABAergic",
                                                   Prediction %in% c("IT", "L4_IT", "L5/6_IT_Car3", "L5/6_NP", "L5_ET", "L6_CT", "L6b") ~ "Glutamatergic", 
                                                   Prediction == "Oligodendrocyte" ~ "Oligodendrocyte", 
                                                   Prediction == "Astrocyte" ~ "Astrocyte", 
                                                   Prediction == "OPC" ~ "OPC", 
                                                   Prediction == "Microglia" ~ "Microglia", 
                                                   Prediction == "Pericyte" ~ "Pericyte", 
                                                   Prediction == "Endothelial" ~ "Endothelial",
                                                   Prediction == "VLMC" ~ "VLMC", 
                                                   .default = Prediction))
    })
    
    res.common.ls <- lapply(res.ls, function(res){
      res_notna <- res %>% drop_na(Prediction)
      res_common <- res_notna[res_notna$True_label %in% shared_cell_types,]
      res_common <- res_common[res_common$Prediction %in% r_cell_types,]
      res_common$True_label <- factor(res_common$True_label, levels = shared_cell_types)
      res_common$Prediction <- factor(res_common$Prediction, levels = c(shared_cell_types, r_only_cell_types))
      return(res_common)
    })
    
    res.qonly.ls <- lapply(res.ls, function(res){
      res_notna <- res %>% replace_na(list(Prediction="missing"))
      res_qonly <- res_notna[res_notna$True_label %in% q_only_cell_types,]
      return(res_qonly)
    })
    
    res.id.ls <- lapply(res.ls, function(res){
      res_notna <- res %>% drop_na(Prediction)
      res_identified <- res_notna[res_notna$Prediction %in% r_cell_types,]
      res_identified$True_label <- factor(res_identified$True_label, levels = union_cell_types)
      res_identified$Prediction <- factor(res_identified$Prediction, levels = union_cell_types)
      return(res_identified)
    })
    
    eval.tbl$assignedCounts <- sapply(res.id.ls, nrow)
    eval.tbl$unassignedCounts <- eval.tbl$CellCounts - eval.tbl$assignedCounts
    eval.tbl$assignedRate <- eval.tbl$assignedCounts/eval.tbl$CellCounts *100
    eval.tbl$unassignedRate <- eval.tbl$unassignedCounts/eval.tbl$CellCounts *100
    
    cm.ls <- lapply(res.common.ls, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
    })
    maxlength = 5
    cm.ls <- cm.ls[!sapply(cm.ls,is.null)]
    acc <- sapply(cm.ls, function(cm){
      n = sum(cm)
      diag = diag(cm)
      accuracy = sum(diag) / n 
    })
    acc = c(acc, rep(NA, maxlength - length(acc)))
    eval.tbl$Accuracy <- acc
    
    wa.f1 <- sapply(cm.ls, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      
      precision = diag / colsums[rownames(cm)]
      recall = diag / rowsums 
      f1 = 2 * precision * recall / (precision + recall) 
      
      f1[is.na(f1)] <- 0
      WA.f1 <- weighted.mean(f1,p)
    })
    
    wa.f1 = c(wa.f1, rep(NA, maxlength - length(wa.f1)))
    eval.tbl$WeightedAveF1 <- wa.f1
    
    precision <- sapply(cm.ls, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      precision = diag / colsums[rownames(cm)]
      precision[is.na(precision)] <- 0
      precision = weighted.mean(precision, p)
    })
    eval.tbl$Precision <- precision
    
    recall <- sapply(cm.ls, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      recall = diag / rowsums
      recall[is.na(recall)] <- 0
      recall = weighted.mean(recall,p)
    })
    eval.tbl$Recall <- recall
    
    rare_cell_type = c("Astrocyte", "Endothelial", "Microglia", "Oligodendrocyte", "OPC", "VLMC")
    abundant_cell_type = c("GABAergic", "Glutamatergic")
    
    res.rare <- lapply(res.common.ls, function(df){
      df.temp <- df[df$True_label %in% intersect(shared_cell_types, rare_cell_type),]
      return(df.temp)
    })
    res.rich <- lapply(res.common.ls, function(df){
      df.temp <- df[df$True_label %in% intersect(shared_cell_types, abundant_cell_type),]
      return(df.temp)
    })
    
    cm.ls.rare <- lapply(res.rare, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
    })
    cm.ls.rich <- lapply(res.rich, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
    })
    
    acc.rare <- sapply(cm.ls.rare, function(cm){
      n = sum(cm)
      diag = diag(cm)
      accuracy = sum(diag) / n 
    })
    eval.tbl$Accuracy.rare <- acc.rare
    
    acc.rich <- sapply(cm.ls.rich, function(cm){
      n = sum(cm)
      diag = diag(cm)
      accuracy = sum(diag) / n 
    })
    eval.tbl$Accuracy.rich <- acc.rich
    
    wa.f1.rare <- sapply(cm.ls.rare, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      
      precision = diag / colsums[rownames(cm)]
      recall = diag / rowsums 
      f1 = 2 * precision * recall / (precision + recall) 
      
      f1[is.na(f1)] <- 0
      WA.f1 <- weighted.mean(f1,p)
    })
    eval.tbl$WeightedAveF1.rare <- wa.f1.rare
    
    wa.f1.rich <- sapply(cm.ls.rich, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      
      precision = diag / colsums[rownames(cm)] 
      recall = diag / rowsums 
      f1 = 2 * precision * recall / (precision + recall) 
      
      f1[is.na(f1)] <- 0
      WA.f1 <- weighted.mean(f1,p)
    })
    eval.tbl$WeightedAveF1.rich <- wa.f1.rich
    
    ari <- sapply(res.common.ls, function(res_identified){
      ARI <- pdfCluster::adj.rand.index(res_identified$True_label, res_identified$Prediction)
      return(ARI)
    })
    eval.tbl$ARI <- ari
    
    vscore <- lapply(res.common.ls, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
      cm.types <- rownames(cm)
      
      res.num <- data.frame(True_label=rep(NA, nrow(res_identified)), Prediction=rep(NA, nrow(res_identified)))
      res.num$True_label <- c(1:length(cm.types), res_identified$True_label)[match(res_identified$True_label, c(cm.types, res_identified$True_label))]
      res.num$Prediction <- c(1:length(cm.types), res_identified$Prediction)[match(res_identified$Prediction, c(cm.types, res_identified$Prediction))]
      v.score <- sabre::vmeasure(res.num$True_label, res.num$Prediction)
      return(unlist(v.score))
    })
    eval.tbl <- cbind(eval.tbl, data.frame(t(sapply(vscore,c))))
    
    honesty <- sapply(res.qonly.ls, function(res){
      honesty <- nrow(res[res$Prediction%in%c("unassigned", "unassinged", "missing"),])/nrow(res) * 100
    })
    
    eval.tbl$honesty <- honesty
    
    eval.tbl <- cbind(eval.tbl, running.time$Time)
    names(eval.tbl)[ncol(eval.tbl)] <- "Time"
    
    saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/03_brain/03_evaluation/", method[i], "/", prop,"/evaluating_result_K5_CT3.csv")
    write.csv(eval.tbl, file = saving.path)
  }
}
evaluation_brain_function_ref_ct4 <- function(method, prop){
  
  index <- readRDS(file = paste0("08_StartOver_NewStory/05_interdatasets/03_brain/01_datasets/index/Index_",prop,".rds"))
  
  q_cell_types <- c("Neuronal", "NonNeuronal")
  r_cell_types <- c("Neuronal", "NonNeuronal")
  shared_cell_types <- intersect(q_cell_types, r_cell_types)
  q_only_cell_types <- setdiff(q_cell_types, shared_cell_types)
  r_only_cell_types <- setdiff(r_cell_types, shared_cell_types)
  union_cell_types <- unique(c(q_cell_types, r_cell_types))
  
  for (i in 1:length(method)) {
    k_index = length(index)
    eval.tbl <- data.frame(Method=rep(method[i], k_index), Iteration=1:k_index, KFold = length(index), 
                           assignedCounts=NA, assignedRate=NA,unassignedCounts=NA, unassignedRate=NA, 
                           Accuracy=NA, WeightedAveF1=NA, ARI=NA)
    
    eval.tbl$CellCounts <- sapply(1:5, function(ind){length(index[[ind]][[1]])})
    
    running.time <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/03_brain/02_labelling/", method[i],"/", prop,"/running_time.csv"), row.names = 1)
    
    res.ls <- vector("list", length(index))
    
    for (j in 1:length(index)) {
      res <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/03_brain/02_labelling/", method[i],"/", prop,"/prediction_results_",j, "_5.csv"), row.names=1)
      names(res.ls)[j] <- paste0("res",j)
      res.ls[[j]] <- res
      rm(res)
    }
    
    res.ls <- lapply(res.ls, function(res){
      rownames(res) <- res[,1] 
      res <- res[,-1]
      res <- res[,c("True_label", "Prediction")]
    })
    
    res.ls <- lapply(res.ls, function(res){
      res <- res %>% mutate(True_label = case_when(True_label %in% c("Chandelier", "Lamp5", "Pax6", "Pvalb", "Sncg", "Sst", "Sst_Chodl", "Vip",
                                                                     "IT", "L4_IT", "L5/6_NP", "L5_ET", "L6_CT", "L6_IT_Car3", "L6b") ~ "Neuronal",
                                                   True_label %in% c("Astrocyte", "Endothelial", "Microglia", "Oligodendrocyte", "OPC", "VLMC") ~ "NonNeuronal",
                                                   .default = True_label))
      
      res <- res %>% mutate(Prediction = case_when(Prediction %in% c("Lamp5", "Pax6", "Pvalb", "Sst", "Vip", 
                                                                     "IT", "L4_IT", "L5/6_IT_Car3", "L5/6_NP", "L5_ET", "L6_CT", "L6b") ~ "Neuronal",
                                                   Prediction %in% c("Astrocyte", "Endothelial", "Microglia", "Oligodendrocyte", "OPC", "Pericyte", "VLMC") ~ "NonNeuronal",
                                                   .default = Prediction))
    })
    
    
    res.id.ls <- lapply(res.ls, function(res){
      res_notna <- res %>% drop_na(Prediction)
      res_identified <- res_notna[res_notna$Prediction %in% r_cell_types,]
      res_identified$True_label <- factor(res_identified$True_label, levels = union_cell_types)
      res_identified$Prediction <- factor(res_identified$Prediction, levels = union_cell_types)
      return(res_identified)
    })
    
    eval.tbl$assignedCounts <- sapply(res.id.ls, nrow)
    eval.tbl$unassignedCounts <- eval.tbl$CellCounts - eval.tbl$assignedCounts
    eval.tbl$assignedRate <- eval.tbl$assignedCounts/eval.tbl$CellCounts *100
    eval.tbl$unassignedRate <- eval.tbl$unassignedCounts/eval.tbl$CellCounts *100
    
    cm.ls <- lapply(res.id.ls, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
    })
    
    acc <- sapply(cm.ls, function(cm){
      n = sum(cm)
      diag = diag(cm)
      accuracy = sum(diag) / n 
    })
    eval.tbl$Accuracy <- acc
    
    wa.f1 <- sapply(cm.ls, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      
      precision = diag / colsums[rownames(cm)] 
      recall = diag / rowsums 
      f1 = 2 * precision * recall / (precision + recall) 
      
      f1[is.na(f1)] <- 0
      WA.f1 <- weighted.mean(f1,p)
    })
    eval.tbl$WeightedAveF1 <- wa.f1
    
    precision <- sapply(cm.ls, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      precision = diag / colsums[rownames(cm)]
      precision[is.na(precision)] <- 0
      precision = weighted.mean(precision, p)
    })
    eval.tbl$Precision <- precision
    
    recall <- sapply(cm.ls, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      recall = diag / rowsums
      recall[is.na(recall)] <- 0
      recall = weighted.mean(recall,p)
    })
    eval.tbl$Recall <- recall
    
    rare_cell_type = c("NonNeuronal")
    abundant_cell_type = c("Neuronal")
    
    res.rare <- lapply(res.id.ls, function(df){
      df.temp <- df[df$True_label %in% intersect(shared_cell_types, rare_cell_type),]
      return(df.temp)
    })
    res.rich <- lapply(res.id.ls, function(df){
      df.temp <- df[df$True_label %in% intersect(shared_cell_types, abundant_cell_type),]
      return(df.temp)
    })
    
    cm.ls.rare <- lapply(res.rare, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
    })
    cm.ls.rich <- lapply(res.rich, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
    })
    
    acc.rare <- sapply(cm.ls.rare, function(cm){
      n = sum(cm)
      diag = diag(cm)
      accuracy = sum(diag) / n 
    })
    eval.tbl$Accuracy.rare <- acc.rare
    
    acc.rich <- sapply(cm.ls.rich, function(cm){
      n = sum(cm)
      diag = diag(cm)
      accuracy = sum(diag) / n 
    })
    eval.tbl$Accuracy.rich <- acc.rich
    
    wa.f1.rare <- sapply(cm.ls.rare, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      
      precision = diag / colsums[rownames(cm)] 
      recall = diag / rowsums 
      f1 = 2 * precision * recall / (precision + recall) 
      
      f1[is.na(f1)] <- 0
      WA.f1 <- weighted.mean(f1,p)
    })
    eval.tbl$WeightedAveF1.rare <- wa.f1.rare
    
    wa.f1.rich <- sapply(cm.ls.rich, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      
      precision = diag / colsums[rownames(cm)] 
      recall = diag / rowsums 
      f1 = 2 * precision * recall / (precision + recall) 
      
      f1[is.na(f1)] <- 0
      WA.f1 <- weighted.mean(f1,p)
    })
    eval.tbl$WeightedAveF1.rich <- wa.f1.rich
    
    ari <- sapply(res.id.ls, function(res_identified){
      ARI <- pdfCluster::adj.rand.index(res_identified$True_label, res_identified$Prediction)
      return(ARI)
    })
    eval.tbl$ARI <- ari
    
    vscore <- lapply(res.id.ls, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
      cm.types <- rownames(cm)
      
      res.num <- data.frame(True_label=rep(NA, nrow(res_identified)), Prediction=rep(NA, nrow(res_identified)))
      res.num$True_label <- c(1:length(cm.types), res_identified$True_label)[match(res_identified$True_label, c(cm.types, res_identified$True_label))]
      res.num$Prediction <- c(1:length(cm.types), res_identified$Prediction)[match(res_identified$Prediction, c(cm.types, res_identified$Prediction))]
      v.score <- sabre::vmeasure(res.num$True_label, res.num$Prediction)
      return(unlist(v.score))
    })
    eval.tbl <- cbind(eval.tbl, data.frame(t(sapply(vscore,c))))
    
    eval.tbl <- cbind(eval.tbl, running.time$Time)
    names(eval.tbl)[ncol(eval.tbl)] <- "Time"
    
    saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/03_brain/03_evaluation/", method[i], "/", prop,"/evaluating_result_K5_CT4.csv")
    write.csv(eval.tbl, file = saving.path)
  }
}