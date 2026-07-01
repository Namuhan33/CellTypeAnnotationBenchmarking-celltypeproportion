setwd("R:/ANNCOM2023-Q6025")
library(Seurat)
library(dplyr)
library(tidyverse)

method <- c("CHETAH", "scmapcluster", "SingleCellNet", "SingleR", "scClassify", "scPred")

evaluation_pancreas_rough_accuracy <- function(method, prop){
  index <- readRDS(file = paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/01_datasets/index/Index_",prop,".rds"))
  
  for (i in 1:length(method)) {
    k_index = length(index)
    eval.tbl <- data.frame(Method=rep(method[i], k_index), Iteration=1:k_index, KFold = k_index, 
                           RoughAccuracy=NA)
    
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
    
    res.rough.acc.ls <- lapply(res.ls, function(res){
      cells_common_cell_type <- res[res$True_label==res$Prediction,]
    })
    
    eval.tbl$CorrectLabels <- sapply(res.rough.acc.ls, nrow)
    eval.tbl$RoughAccuracy <- eval.tbl$CorrectLabels/eval.tbl$CellCounts*100
    
    saving.path <- paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/03_evaluation/", method[i],"/", prop,"/evaluating_rough_accuracy.csv")
    write.csv(eval.tbl, file = saving.path)
  }
}

for(i in 1:14){
  prp = paste0("P", i)
  evaluation_pancreas_rough_accuracy(method, prp)
}

combine_rough_acc <- function(method){
  results.list <- list()
  for(method in method){
    for (prp in 1:14) {
      path <- paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/03_evaluation/", method, "/P", prp, "/")
      df.1 <- read.csv(paste0(path, "evaluating_rough_accuracy.csv"), row.names = 1)
      df.1$proportion=paste0("Prop", prp)
      results.list <- rlist::list.append(results.list, df.1)
    }
  }
  return(results.list)
}
RoughAcc <- combine_rough_acc(method)
RoughAcc <- do.call(rbind, RoughAcc)

prop.labs <- c("Endo:Non(10:1)[#Q=#R]", "Endo:Non(5:1)[#Q=#R]", "Endo:Non(2:1)[#Q=#R]", "Endo:Non(1:1)[#Q=#R]", "Endo:Non(0.5:1)[#Q=#R]", "Endo:Non(0.2:1)[#Q=#R]", "Endo:Non(0.1:1)[#Q=#R]",
               "Endo:Non(10:1)[#Q<#R]", "Endo:Non(5:1)[#Q<#R]", "Endo:Non(2:1)[#Q<#R]", "Endo:Non(1:1)[#Q<#R]", "Endo:Non(0.5:1)[#Q<#R]", "Endo:Non(0.2:1)[#Q<#R]", "Endo:Non(0.1:1)[#Q<#R]")
names(prop.labs) <- paste0("Prop", 1:14)


ggplot(RoughAcc, aes(x=proportion, y=RoughAccuracy, fill=Method)) + 
  geom_boxplot() +
  facet_wrap(~Method) +
  theme_bw() +
  scale_fill_manual(values = c("#9d92d6", "#faa9b4", "#fffd7c", "#c69189", "#ffb169", "#83d0ff")) +
  scale_x_discrete(labels = prop.labs,
                   guide = guide_axis(angle = 90),
                   limits=factor(paste0("Prop", 1:14), levels = paste0("Prop", 1:14))) +
  xlab(label = "") +
  ylab(label = "Rough Accuracy") +
  theme(strip.text.x = element_text(size = 15), 
        text = element_text(size=20), 
        axis.text=element_text(size=12),
        axis.ticks.y=element_blank(),
        panel.spacing = unit(0.1, "cm"),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank())
