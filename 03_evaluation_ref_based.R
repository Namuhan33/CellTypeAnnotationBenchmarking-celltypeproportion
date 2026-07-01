setwd("R:/ANNCOM2023-Q6025")
library(Seurat)
library(tidyverse)
library(dplyr)

method <- c("CHETAH", "scmapcluster", "SingleCellNet", "SingleR", "scClassify", "scPred")

evaluation_function_ref <- function(method, prop){
  index <- readRDS(file = paste0("08_StartOver_NewStory/01_datasets/index/Index_", prop, ".rds"))
  
  for (i in 1:length(method)) {
    k_index = length(index)
    
    eval.tbl <- data.frame(Method=rep(method[i], k_index), Iteration=1:k_index, KFold = k_index, 
                           assignedCounts=NA, assignedRate=NA,unassignedCounts=NA, unassignedRate=NA, 
                           Accuracy=NA, WeightedAveF1=NA, ARI=NA)
    
    eval.tbl$CellCounts <- sapply(index, lengths)[1,]
    
    running.time <- read.csv(file=paste0("08_StartOver_NewStory/02_labelling/", method[i],"/", prop,"/running_time.csv"), row.names = 1)
    
    res.ls <- vector("list", length(index))
    
    for (j in 1:length(index)) {
      res <- read.csv(file=paste0("08_StartOver_NewStory/02_labelling/", method[i],"/", prop,"/prediction_results_",j, "_", length(index),".csv"), row.names=1)
      names(res.ls)[j] <- paste0("res",j)
      res.ls[[j]] <- res
      rm(res)
    }
    
    res.ls <- res.ls[!sapply(res.ls,is.null)]
    
    res.ls <- lapply(res.ls, function(res){
      rownames(res) <- res[,1] 
      res <- res[,-1]
      res <- res[,c("True_label", "Prediction")]
    })
    
    res.id.ls <- lapply(res.ls, function(res){
      
      res_notna <- res %>% drop_na(Prediction)
      res_identified <- res_notna[res_notna$Prediction %in% c("B_cell","Basal","DC","G2MS","LC1","LC2","Macrophage","Neutrophil","NK","T_cell"),]

      res_identified$True_label <- factor(res_identified$True_label, levels = c("B_cell","Basal","DC","G2MS","LC1","LC2","Macrophage","Neutrophil","NK","T_cell"))
      res_identified$Prediction <- factor(res_identified$Prediction, levels = c("B_cell","Basal","DC","G2MS","LC1","LC2","Macrophage","Neutrophil","NK","T_cell"))
      return(res_identified)
    })
    
    eval.tbl$assignedCounts <- sapply(res.id.ls, nrow)
    eval.tbl$unassignedCounts <- eval.tbl$CellCounts - eval.tbl$assignedCounts
    eval.tbl$assignedRate <- eval.tbl$assignedCounts/eval.tbl$CellCounts *100
    eval.tbl$unassignedRate <- eval.tbl$unassignedCounts/eval.tbl$CellCounts *100
    
    cm.ls <- lapply(res.id.ls, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
      rm_empty_class <- function(df, classes = c("B_cell","Basal","DC","G2MS","LC1","LC2","Macrophage","Neutrophil","NK","T_cell")){
        df <- as.data.frame.matrix(df)
        removing <- c()
        for(k in 1:length(classes)){
          if(sum(df[k,])==sum(df[,k]) & sum(df[k,])==0){
            removing <- append(removing, k)
          }
        }
        if(!is.null(removing)){
          df <- df[ , -removing]
          df <- df[ -removing, ]
        }
        
        df <- data.matrix(df)
        return(df)
      }
      cm <- rm_empty_class(df=cm)
      return(cm)
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
      precision = diag / colsums
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
    
    if(prop %in% paste0("P", c(1, 2))){
      rare_cell_type = c("B_cell", "Basal", "DC", "G2MS", "Neutrophil", "NK", "T_cell")
      abundant_cell_type = c("LC1", "LC2", "Macrophage")
    }
    if(prop %in% paste0("P", c(3,4, 6:12, 13:26, 27:40))){
      rare_cell_type = c("B_cell", "Basal", "DC", "G2MS", "Neutrophil", "NK")
      abundant_cell_type = c("LC1", "LC2", "Macrophage", "T_cell")
    }
    if(prop %in% paste0("P", c(5))){
      rare_cell_type = c("B_cell", "Basal", "DC", "G2MS", "NK")
      abundant_cell_type = c("LC1", "LC2", "Macrophage", "Neutrophil", "T_cell")
    }

    
    res.rare <- lapply(res.id.ls, function(df){
      df.temp <- df[df$True_label %in% rare_cell_type,]
      return(df.temp)
    })
    res.rich <- lapply(res.id.ls, function(df){
      df.temp <- df[df$True_label %in% abundant_cell_type,]
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
    
    precision.rare <- sapply(cm.ls.rare, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      precision = diag / colsums
      precision = weighted.mean(precision, p)
    })
    eval.tbl$Precision.rare <- precision.rare
    
    precision.rich <- sapply(cm.ls.rich, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      precision = diag / colsums
      precision = weighted.mean(precision, p)
    })
    eval.tbl$Precision.rich <- precision.rich
    
    recall.rare <- sapply(cm.ls.rare, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      recall = diag / rowsums 
      recall = weighted.mean(recall,p)
    })
    eval.tbl$Recall.rare <- recall.rare
    
    recall.rich <- sapply(cm.ls.rich, function(cm){
      n = sum(cm)
      diag = diag(cm)  
      rowsums = apply(cm, 1, sum)
      colsums = apply(cm, 2, sum)
      p = rowsums / n 
      recall = diag / rowsums 
      recall = weighted.mean(recall,p)
    })
    eval.tbl$Recall.rich <- recall.rich
    
    ari <- sapply(res.id.ls, function(res_identified){
      ARI <- pdfCluster::adj.rand.index(res_identified$True_label, res_identified$Prediction)
      return(ARI)
    })
    eval.tbl$ARI <- ari
    
    vscore <- lapply(res.id.ls, function(res_identified){
      cm <- as.matrix(table(Actual = res_identified$True_label, Predicted = res_identified$Prediction)) # create the confusion matrix
      rm_empty_class <- function(df, classes = c("B_cell","Basal","DC","G2MS","LC1","LC2","Macrophage","Neutrophil","NK","T_cell")){
        df <- as.data.frame.matrix(df)
        removing <- c()
        for(k in 1:length(classes)){
          if(sum(df[k,])==sum(df[,k]) & sum(df[k,])==0){
            removing <- append(removing, k)
          }
        }
        if(!is.null(removing)){
          df <- df[ , -removing]
          df <- df[ -removing, ]
        }
        
        df <- data.matrix(df)
        return(df)
      }
      cm <- rm_empty_class(df=cm)
      
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
    
    saving.path <- paste0("08_StartOver_NewStory/03_evaluation/", method[i], "/", prop,"/evaluating_result_K", length(index),".csv")
    write.csv(eval.tbl, file = saving.path)
  }
}
evaluation_function_ref.G <- function(method, prop){
  index <- readRDS(file = paste0("08_StartOver_NewStory/01_datasets/index/Index_", prop, ".rds"))
  
  for (i in 1:length(method)) {
    k_index = length(index)
  
    eval.tbl <- data.frame(Method=rep(method[i], k_index), Iteration=1:k_index, KFold = length(index), 
                           assignedCounts=NA, assignedRate=NA,unassignedCounts=NA, unassignedRate=NA, 
                           Accuracy=NA, WeightedAveF1=NA, ARI=NA)
    
    eval.tbl$CellCounts <- sapply(index, lengths)[1,]
    
    running.time <- read.csv(file=paste0("08_StartOver_NewStory/02_labelling/", method[i],"/", prop,"/running_time.csv"), row.names = 1)
    
    res.ls <- vector("list", length(index))
    
    for (j in 1:length(index)) {
      res <- read.csv(file=paste0("08_StartOver_NewStory/02_labelling/", method[i],"/", prop,"/prediction_results_",j, "_", length(index),".csv"), row.names=1)
      names(res.ls)[j] <- paste0("res",j)
      res.ls[[j]] <- res
      rm(res)
    }
    
    res.ls <- res.ls[!sapply(res.ls,is.null)]
    
    res.ls <- lapply(res.ls, function(res){
      rownames(res) <- res[,1] 
      res <- res[,-1]
      res <- res[,c("True_label", "Prediction")]
    })
    
    res.id.ls <- lapply(res.ls, function(res){
      res_notna <- res %>% drop_na(Prediction)
      res_identified <- res_notna[res_notna$Prediction %in% c("B_cell","Basal","DC","G2MS","LC1","LC2","Macrophage","Neutrophil","NK","T_cell"), ]
      res_identified$True_label <- factor(res_identified$True_label, levels = c("B_cell","Basal","DC","G2MS","LC1","LC2","Macrophage","Neutrophil","NK","T_cell"))
      res_identified$Prediction <- factor(res_identified$Prediction, levels = c("B_cell","Basal","DC","G2MS","LC1","LC2","Macrophage","Neutrophil","NK","T_cell"))
      
      res_identified <- res_identified %>% mutate(T.G = case_when(True_label %in% c("B_cell","DC","Macrophage","Neutrophil","NK","T_cell") ~ "Immune",
                                                                  True_label %in% c("Basal","G2MS","LC1","LC2") ~ "Epithelial"),
                                                  P.G = case_when(Prediction %in% c("B_cell","DC","Macrophage","Neutrophil","NK","T_cell") ~ "Immune",
                                                                  Prediction %in% c("Basal","G2MS","LC1","LC2") ~ "Epithelial"))
      res.G <- data.frame(True_label = res_identified$T.G, Prediction = res_identified$P.G, row.names = rownames(res_identified))
      
      return(res.G)
    })
    
    eval.tbl$assignedCounts <- sapply(res.id.ls, nrow)
    eval.tbl$unassignedCounts <- eval.tbl$CellCounts - eval.tbl$assignedCounts
    eval.tbl$assignedRate <- eval.tbl$assignedCounts/eval.tbl$CellCounts *100
    eval.tbl$unassignedRate <- eval.tbl$unassignedCounts/eval.tbl$CellCounts *100
    
    cm.ls <- lapply(res.id.ls, function(res_identified){
      cm <- as.matrix(table(Actual = factor(res_identified$True_label, levels = c("Epithelial", "Immune")), Predicted = factor(res_identified$Prediction, levels = c("Epithelial", "Immune")))) # create the confusion matrix
      return(cm)
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
    
    if(prop %in% paste0("P", c(1, 2, 3))){
      rare_cell_type = c("IMMUNE")
      abundant_cell_type = c("EPITHELIAL")
    }
    if(prop %in% paste0("P", c(4, 5, 6:12, 13:26, 27:40))){
      rare_cell_type = c("EPITHELIAL")
      abundant_cell_type = c("IMMUNE")
    }
    
    res.rare <- lapply(res.id.ls, function(df){
      df.temp <- df[df$True_label %in% rare_cell_type,]
      return(df.temp)
    })
    res.rich <- lapply(res.id.ls, function(df){
      df.temp <- df[df$True_label %in% abundant_cell_type,]
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
    
    saving.path <- paste0("08_StartOver_NewStory/03_evaluation/", method[i], "/", prop,"/evaluating_result_K", length(index),"_General.csv")
    write.csv(eval.tbl, file = saving.path)
  }
}
