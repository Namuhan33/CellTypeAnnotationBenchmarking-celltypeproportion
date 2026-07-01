setwd("R:/ANNCOM2023-Q6025")
library(Seurat)
library(dplyr)
library(tidyverse)

method <- c("CHETAH", "scmapcluster", "SingleCellNet", "SingleR", "scClassify", "scPred")


evaluation_pancreas_function_ref_ct1 <- function(method, prop){
  index <- readRDS(file = paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/01_datasets/index/Index_",prop,".rds"))
  
  q_cell_types <- c("Delta","Alpha","Gamma","Ductal","Acinar","Beta","Unclassified_endocrine","Endocrine_coexpression",
                    "MHC_class_II","PSC","Endothelial","Epsilon","Mast","Unclassified_exocrine")
  r_cell_types <- c("Acinar","Beta","Delta","PSC","Ductal","Alpha","Epsilon","Gamma",      
                    "Endothelial","Macrophage","Schwann","Mast","T_cell")
  
  shared_cell_types <- intersect(q_cell_types, r_cell_types)
  q_only_cell_types <- setdiff(q_cell_types, shared_cell_types)
  r_only_cell_types <- setdiff(r_cell_types, shared_cell_types)
  union_cell_types <- unique(c(q_cell_types, r_cell_types))
  
  for (i in 1:length(method)) {
    k_index = length(index)
    eval.tbl <- data.frame(Method=rep(method[i], k_index), Iteration=1:k_index, KFold = k_index, 
                           assignedCounts=NA, assignedRate=NA,unassignedCounts=NA, unassignedRate=NA, 
                           Accuracy=NA, WeightedAveF1=NA, ARI=NA)
    
    running.time <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/02_labelling/", method[i],"/", prop,"/running_time.csv"), row.names = 1)
    
    eval.tbl$CellCounts <- sapply(1:5, function(ind){length(index[[ind]][[1]])})
    
    res.ls <- vector("list", length(index))
    
    for (j in 1:length(index)) {
      res <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/02_labelling/", method[i],"/", prop,"/prediction_results_",j, "_5.csv"), row.names=1)
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
    
    rare_cell_type = c("Endocrine_coexpression", "Endothelial", "Epsilon", "Mast", "MHC_class_II", "PSC", "Unclassified_endocrine", "Unclassified_exocrine")
    abundant_cell_type = c("Acinar", "Alpha", "Beta", "Delta", "Ductal", "Gamma")
    
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
    
    saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/03_evaluation/", method[i], "/", prop,"/evaluating_result_K5.csv")
    write.csv(eval.tbl, file = saving.path)
  }
}
evaluation_pancreas_function_ref_ct2 <- function(method, prop){
  index <- readRDS(file = paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/01_datasets/index/Index_",prop,".rds"))

  q_cell_types <- c("Endocrine","Ductal","Exocrine","Myeloid","PSC","Endothelial")
  r_cell_types <- c("Exocrine","Endocrine","PSC","Ductal","Endothelial","Myeloid","Schwann","Lymphoid")
  
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
    
    running.time <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/02_labelling/", method[i],"/", prop,"/running_time.csv"), row.names = 1)
    
    res.ls <- vector("list", length(index))
    
    for (j in 1:length(index)) {
      res <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/02_labelling/", method[i],"/", prop,"/prediction_results_",j, "_5.csv"), row.names=1)
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
      res <- res %>% mutate(True_label = case_when(True_label %in% c("Alpha", "Beta" ,"Gamma", "Delta", "Epsilon", "Unclassified_endocrine", "Endocrine_coexpression") ~ "Endocrine",
                                                   True_label %in% c("Unclassified_exocrine", "Acinar") ~ "Exocrine",
                                                   True_label %in% c("MHC_class_II" ,"Mast") ~ "Myeloid",
                                                   True_label == "Endothelial" ~ "Endothelial",
                                                   True_label == "PSC" ~ "PSC",
                                                   True_label == "Ductal" ~ "Ductal"))
      
      res <- res %>% mutate(Prediction = case_when(Prediction %in% c("Alpha", "Beta" ,"Gamma", "Delta", "Epsilon") ~ "Endocrine",
                                                   Prediction %in% c("Macrophage", "Mast") ~ "Myeloid",
                                                   Prediction == "Acinar" ~ "Exocrine",
                                                   Prediction == "Endothelial" ~ "Endothelial",
                                                   Prediction == "PSC" ~ "PSC",
                                                   Prediction == "Ductal" ~ "Ductal",
                                                   Prediction == "T_cell" ~ "Lymphoid",
                                                   Prediction == "Schwann" ~ "Schwann",
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
    
    rare_cell_type = c("Endothelial", "Myeloid", "PSC")
    abundant_cell_type = c("Ductal", "Endocrine", "Exocrine")
    
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
    
    saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/03_evaluation/", method[i], "/", prop,"/evaluating_result_K5_CT2.csv")
    write.csv(eval.tbl, file = saving.path)
  }
}
evaluation_pancreas_function_ref_ct3 <- function(method, prop){
  index <- readRDS(file = paste0("08_StartOver_NewStory/05_interdatasets/01_pbmc/01_datasets/index/Index_",prop,".rds"))
  
  q_cell_types <- c("Endocrine", "Non-Endocrine")
  r_cell_types <- c("Endocrine", "Non-Endocrine")
  
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
    
    running.time <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/02_labelling/", method[i],"/", prop,"/running_time.csv"), row.names = 1)
    
    res.ls <- vector("list", length(index))
    
    for (j in 1:length(index)) {
      res <- read.csv(file=paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/02_labelling/", method[i],"/", prop,"/prediction_results_",j, "_5.csv"), row.names=1)
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
      res <- res %>% mutate(True_label = case_when(True_label %in% c("Alpha", "Beta" ,"Gamma", "Delta", "Epsilon", "Unclassified_endocrine", "Endocrine_coexpression") ~ "Endocrine",
                                                   True_label %in% c("Acinar", "Ductal", "Endothelial", "Mast", "MHC_class_II", "PSC", "Unclassified_exocrine") ~ "Non-Endocrine"))
      
      res <- res %>% mutate(Prediction = case_when(Prediction %in% c("Alpha", "Beta", "Delta","Gamma", "Epsilon") ~ "Endocrine",
                                                   Prediction %in% c("Acinar", "Ductal", "Macrophage", "Mast", "Endothelial", "PSC", "Schwann", "T_cell") ~ "Non-Endocrine",
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
    
    rare_cell_type = c("Non-Endocrine")
    abundant_cell_type = c("Endocrine")
    
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
    
    saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/03_evaluation/", method[i], "/", prop,"/evaluating_result_K5_CT3.csv")
    write.csv(eval.tbl, file = saving.path)
  }
}