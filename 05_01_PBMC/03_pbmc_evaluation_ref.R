setwd("R:/ANNCOM2023-Q6025")
library(Seurat)
library(dplyr)
library(tidyverse)

method <- c("CHETAH", "scmapcluster", "SingleCellNet", "SingleR", "scClassify", "scPred")

evaluation_pbmc_function_ref_ct1 <- function(method, prop){
  index <- readRDS(file = paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/01_datasets/index/Index_",prop,".rds"))

  q_cell_types <- c("RBC","B","CD14_monocyte","CD8_mem_T","CD4_mem_T","CD4_naive_T",    
                    "Platelet","CD4_T","NK","Neutrophil","CD16_monocyte","CD8_cytotoxic_T",
                    "gd_T","DC","Granulocyte","Eosinophil")
  r_cell_types <- c("CD8_cytotoxic_T","CD8_naive_T","CD4_mem_T","B","CD4_T","NK","CD4_naive_T","CD34","DC","CD14_monocyte")
  
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
    
    running.time <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/", method[i],"/", prop,"/running_time.csv"), row.names = 1)
    
    res.ls <- vector("list", length(index))
    
    for (j in 1:length(index)) {
      res <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/", method[i],"/", prop,"/prediction_results_",j, "_5.csv"), row.names=1)
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
    
    
    rare_cell_type = c("CD16_monocyte", "CD4_T", "CD8_cytotoxic_T", "DC", "Eosinophil", 
                       "gd_T", "Granulocyte", "Neutrophil", "Platelet", "RBC")
    abundant_cell_type = c("B", "CD14_monocyte", "CD8_mem_T", "CD4_mem_T", "CD4_naive_T", "NK")
    
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
    
    saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/03_evaluation/", method[i],"/", prop,"/evaluating_result_K5.csv")
    write.csv(eval.tbl, file = saving.path)
  }
}
evaluation_pbmc_function_ref_ct2 <- function(method, prop){
  index <- readRDS(file = paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/01_datasets/index/Index_",prop,".rds"))
  
  q_cell_types <- c("RBC","B","Monocyte","CD8_T","CD4_T","Platelet","NK","Granulocyte","gd_T","DC")
  r_cell_types <- c("CD8_T","CD4_T","B","NK","CD34","DC","Monocyte")
  
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
    
    running.time <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/", method[i],"/", prop,"/running_time.csv"), row.names = 1)
    
    res.ls <- vector("list", length(index))
    
    for (j in 1:length(index)) {
      res <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/", method[i],"/", prop,"/prediction_results_",j, "_5.csv"), row.names=1)
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
      res <- res %>% mutate(True_label = case_when(True_label == "B" ~ "B",
                                                   True_label %in% c("CD14_monocyte", "CD16_monocyte") ~ "Monocyte",
                                                   True_label %in% c("CD4_mem_T", "CD4_naive_T", "CD4_T") ~ "CD4_T",
                                                   True_label %in% c("CD8_cytotoxic_T", "CD8_mem_T") ~ "CD8_T",
                                                   True_label == "DC" ~ "DC",
                                                   True_label == "gd_T" ~ "gd_T",
                                                   True_label %in% c("Eosinophil", "Granulocyte", "Neutrophil") ~ "Granulocyte",
                                                   True_label == "NK" ~ "NK",
                                                   True_label == "Platelet" ~ "Platelet",
                                                   True_label == "RBC" ~ "RBC"))
      
      res <- res %>% mutate(Prediction = case_when(Prediction == "B" ~ "B",
                                                   Prediction == "CD14_monocyte" ~ "Monocyte",
                                                   Prediction == "CD34" ~ "CD34",
                                                   Prediction %in% c("CD4_mem_T", "CD4_naive_T", "CD4_T") ~ "CD4_T",
                                                   Prediction %in% c("CD8_cytotoxic_T", "CD8_naive_T") ~ "CD8_T",
                                                   Prediction == "DC" ~ "DC",
                                                   Prediction == "NK" ~ "NK",
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
    
    rare_cell_type = c("DC", "gd_T", "Granulocyte", "Platelet", "RBC")
    abundant_cell_type = c("B", "Monocyte", "CD8_T", "CD4_T", "NK")
   
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
    
    saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/03_evaluation/", method[i],"/", prop,"/evaluating_result_K5_CT2.csv")
    write.csv(eval.tbl, file = saving.path)
  }
}
evaluation_pbmc_function_ref_ct3 <- function(method, prop){
  index <- readRDS(file = paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/01_datasets/index/Index_",prop,".rds"))
  
  q_cell_types <- c("Lymphoid", "Myeloid")
  r_cell_types <- c("Lymphoid", "Myeloid")
  
  shared_cell_types <- intersect(q_cell_types, r_cell_types)
  union_cell_types <- unique(c(q_cell_types, r_cell_types))
  q_only_cell_types <- setdiff(q_cell_types, shared_cell_types)
  r_only_cell_types <- setdiff(r_cell_types, shared_cell_types)
  
  for (i in 1:length(method)) {
    k_index = length(index)
    eval.tbl <- data.frame(Method=rep(method[i], k_index), Iteration=1:k_index, KFold = length(index), 
                           assignedCounts=NA, assignedRate=NA,unassignedCounts=NA, unassignedRate=NA, 
                           Accuracy=NA, WeightedAveF1=NA, ARI=NA)
    
    eval.tbl$CellCounts <- sapply(1:5, function(ind){length(index[[ind]][[1]])})
    
    running.time <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/", method[i],"/", prop,"/running_time.csv"), row.names = 1)
    
    res.ls <- vector("list", length(index))
    
    for (j in 1:length(index)) {
      res <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/02_labelling/", method[i],"/", prop,"/prediction_results_",j, "_5.csv"), row.names=1)
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
      res <- res %>% mutate(True_label = case_when(True_label == "B" ~ "B",
                                                   True_label %in% c("CD14_monocyte", "CD16_monocyte") ~ "Monocyte",
                                                   True_label %in% c("CD4_mem_T", "CD4_naive_T", "CD4_T") ~ "CD4_T",
                                                   True_label %in% c("CD8_cytotoxic_T", "CD8_mem_T") ~ "CD8_T",
                                                   True_label == "DC" ~ "DC",
                                                   True_label == "gd_T" ~ "gd_T",
                                                   True_label %in% c("Eosinophil", "Granulocyte", "Neutrophil") ~ "Granulocyte",
                                                   True_label == "NK" ~ "NK",
                                                   True_label == "Platelet" ~ "Platelet",
                                                   True_label == "RBC" ~ "RBC"))
      
      res <- res %>% mutate(Prediction = case_when(Prediction == "B" ~ "B",
                                                   Prediction == "CD14_monocyte" ~ "Monocyte",
                                                   Prediction == "CD34" ~ "CD34",
                                                   Prediction %in% c("CD4_mem_T", "CD4_naive_T", "CD4_T") ~ "CD4_T",
                                                   Prediction %in% c("CD8_cytotoxic_T", "CD8_naive_T") ~ "CD8_T",
                                                   Prediction == "DC" ~ "DC",
                                                   Prediction == "NK" ~ "NK",
                                                   .default = Prediction))
      res <- res %>% mutate(True_label = case_when(True_label %in% c("CD4_T", "CD8_T", "B", "gd_T", "NK") ~ "Lymphoid",
                                                   True_label %in% c("DC", "Granulocyte", "Monocyte", "Platelet", "RBC") ~ "Myeloid"))
      
      res <- res %>% mutate(Prediction = case_when(Prediction %in% c("B", "CD4_T", "CD8_T", "NK") ~ "Lymphoid",
                                                   Prediction %in% c("Monocyte", "CD34", "DC") ~ "Myeloid",
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
      
      precision = diag / colsums 
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
    
    
    rare_cell_type = c("Myeloid")
    abundant_cell_type = c("Lymphoid")
    
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
    
    saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/03_evaluation/", method[i],"/", prop,"/evaluating_result_K", length(index),"_CT3.csv")
    write.csv(eval.tbl, file = saving.path)
  }
}

for(i in 1:14){
  prp = paste0("P", i)
  evaluation_pbmc_function_ref_ct1(method, prp)
  evaluation_pbmc_function_ref_ct2(method, prp)
  evaluation_pbmc_function_ref_ct3(method, prp)
}

