setwd("R:/ANNCOM2023-Q6025")
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ggh4x)
library(tibble)
library(ggforce)
library(ggbreak)
library(ggrepel)

# method.ref <- c("CHETAH", "scmapcell", "singleCellNet", "SingleR", "scAnno", "scClassify")
method <- c("CHETAH", "scmapcluster", "SingleCellNet", "SingleR", "scClassify", "scPred")

vector.labels <- c("QR(10:1)", "QR(5:1)", "QR(2:1)", "QR(1:1)", "QR(0.5:1)",
                   "R(10:1)", "R(5:1)","R(2:1)", "R(1:1)", "R(0.5:1)", "R(0.2:1)", "R(0.1:1)",
                   "R(10:1)*1.5", "R(5:1)*1.5","R(2:1)*1.5", "R(1:1)*1.5", "R(0.5:1)*1.5", "R(0.2:1)*1.5", "R(0.1:1)*1.5",
                   "R(10:1)*2", "R(5:1)*2","R(2:1)*2", "R(1:1)*2", "R(0.5:1)*2", "R(0.2:1)*2", "R(0.1:1)*2",
                   "LC2:LC1(4:1)", "LC2:LC1(2.5:1)", "LC2:LC1(2:1)", "LC2:LC1(1:1)", "LC2:LC1(0.5:1)", "LC2:LC1(0.4:1)", "LC2:LC1(0.25:1)",
                   "Mac:T(4:1)", "Mac:T(2.5:1)", "Mac:T(2:1)", "Mac:T(1:1)", "Mac:T(0.5:1)", "Mac:T(0.4:1)", "Mac:T(0.27:1)")
names(vector.labels) <- paste0("Prop", 1:40)

combine_everything <- function(methods){
  results.list <- list()
  for(method in methods){
    path.1 <- paste0("08_StartOver_NewStory/03_evaluation/", method)
    for (prp in 1:length(list.files(path.1))) {
      path.2 <- paste0(path.1, "/P", prp, "/")

      df.1 <- read.csv(paste0(path.2, "evaluating_result_K5.csv"), row.names = 1)
      df.2 <- read.csv(paste0(path.2, "evaluating_result_K5_General.csv"), row.names = 1)

      df.1 <- df.1[,c("assignedRate", "WeightedAveF1", "CellCounts", "Precision", "Recall", "WeightedAveF1.rare", "WeightedAveF1.rich", "Time")]
      df.2 <- df.2[,c("assignedRate", "WeightedAveF1", "CellCounts", "WeightedAveF1.rare", "WeightedAveF1.rich", "Time")]
      temp.df <- data.frame(Precision = rep(0,nrow(df.2)), Recall = rep(0,nrow(df.2)))
      df.2 <- do.call("cbind", list(df.2[,1:3], temp.df, df.2[,4:ncol(df.2)]))

      vec.1 <- colMeans(df.1)
      vec.2 <- colMeans(df.2)

      vec.1 <- append(vec.1, c(method=method, proportion=paste0("Prop", prp), level="D"))
      vec.2 <- append(vec.2, c(method=method, proportion=paste0("Prop", prp), level="G"))

      results.list <- rlist::list.append(results.list, vec.1)
      results.list <- rlist::list.append(results.list, vec.2)
    }
  }
  return(results.list)
}

ref.list <- combine_everything(methods = method)

results <- data.frame(t(sapply(ref.list, c)))
results$proportion <- factor(results$proportion, levels= paste0("Prop",c(1:5, 6:26, 27:40)))

results <- results %>%
  mutate(method = replace(method, method == "scmapcluster", "scmap-cluster"))

results$level <- factor(results$level, levels= c("G", "D"))

# bar/col plot for F1 #####
metrics <- c("WeightedAveF1")
level.labs <- c("Cell-subtype", "Main cell-type")
names(level.labs) <- c("D", "G")
df <- results[,c(metrics, "method", "proportion", "level")]
df[,metrics] <- as.numeric(df[,metrics])

# session 1
df <- df[df$proportion %in% paste0("Prop", 1:5),] # 1:5, 6:26, 27:40

ggplot(data=df, mapping=aes(x=method, y=WeightedAveF1, fill=proportion)) +
  geom_bar(stat = "identity", position=position_dodge(0.8), width=0.8) +
  xlab(label = "") +
  ylab(label = "Weighted Ave F1 Score") +
  labs(fill= "Proportions") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size = 10),
        text = element_text(size=20))+
  geom_text(aes(label = round(WeightedAveF1, 3)), position = position_dodge(.8), size = 4, hjust=1) +
  facet_grid(.~level , scales = "free_x", space = "free_x", labeller = labeller(level = level.labs)) +
  # facet_grid(level~ , scales = "free_x", space = "free_x", labeller = labeller(level = level.labs)) +
  coord_flip() +
  scale_fill_manual(values=c("#9d92d6", "#66b0f0", "#fffd7c", "#c69189", "#ffb169", "#83d0ff", "#046C9A",
                             "#faa9b4", "#FAD510","#EABE94", "#9986A5", "#CCC591", "#F98400", "#c3abdb",
                             "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#E79805", "#EC7A05"),
                    guide=guide_legend(reverse = TRUE),
                    labels = vector.labels)

# session 2
df <- results[,c(metrics, "method", "proportion", "level")]
df[,metrics] <- as.numeric(df[,metrics])
df <- df[df$proportion %in% paste0("Prop", 6:26),] # 1:5, 6:12, 13:26, 27:40
df <- df %>% mutate(cellcounts = case_when(proportion %in% paste0("Prop", 6:12) ~ "Same cell counts",
                                     proportion %in% paste0("Prop", 13:19) ~ "1.5 times cell counts",
                                     proportion %in% paste0("Prop", 20:26) ~ "2 times cell counts"))
df$cellcounts <- factor(df$cellcounts, levels = c("Same cell counts", "1.5 times cell counts","2 times cell counts"))
ggplot(data=df, mapping=aes(x=method, y=WeightedAveF1, fill= proportion)) +
  geom_bar(stat = "identity", position=position_dodge(0.8), width=0.8) +
  xlab(label = "") +
  ylab(label = "Weighted Ave F1 Score") +
  labs(fill= "Proportions") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size = 10),
        text = element_text(size=20))+
  geom_text(aes(label = round(WeightedAveF1, 3)), position = position_dodge(.8), size = 4, hjust=1) +
  facet_grid(level~cellcounts, scales = "free_x", space = "free_x", labeller = labeller(level = level.labs)) + # wide version
  coord_flip() +
  scale_fill_manual(values=c("#9d92d6", "#66b0f0", "#fffd7c", "#c69189", "#ffb169", "#83d0ff", "#046C9A",
                             "#faa9b4", "#FAD510","#EABE94", "#9986A5", "#CCC591", "#F98400", "#c3abdb",
                             "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#E79805", "#EC7A05"),
                    guide=guide_legend(reverse = TRUE),
                    labels = vector.labels) +
  guides(fill = guide_legend(ncol = 1))

ggplot(data=df, mapping=aes(x=method, y=WeightedAveF1, fill= proportion)) +
  geom_bar(stat = "identity", position=position_dodge(0.8), width=0.8) +
  xlab(label = "") +
  ylab(label = "Weighted Ave F1 Score") +
  labs(fill= "Proportions") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size = 10),
        text = element_text(size=20))+
  geom_text(aes(label = round(WeightedAveF1, 3)), position = position_dodge(.8), size = 4, hjust=1) +
  facet_grid(cellcounts~level , scales = "free_x", space = "free_x", labeller = labeller(level = level.labs)) + # long version
  coord_flip() +
  scale_fill_manual(values=c("#9d92d6", "#66b0f0", "#fffd7c", "#c69189", "#ffb169", "#83d0ff", "#046C9A",
                             "#faa9b4", "#FAD510","#EABE94", "#9986A5", "#CCC591", "#F98400", "#c3abdb",
                             "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#E79805", "#EC7A05"),
                    guide=guide_legend(reverse = TRUE),
                    labels = vector.labels)

# session 3
df <- results[,c(metrics, "method", "proportion", "level")]
df[,metrics] <- as.numeric(df[,metrics])
df <- df[df$proportion %in% paste0("Prop", 27:40),] # 1:5, 6:12, 13:26, 27:40
df <- df %>% mutate(Type = case_when(proportion %in% paste0("Prop", 27:33) ~ "Epithelial-subtypes",
                                     proportion %in% paste0("Prop", 34:40) ~ "Immune-subtypes"))
df$Type <- factor(df$Type, levels = c("Epithelial-subtypes", "Immune-subtypes"))
df$proportion <- factor(df$proportion, levels = paste0("Prop", c(34:40,27:33)))
ggplot(data=df, mapping=aes(x=method, y=WeightedAveF1, fill= proportion)) +
  geom_bar(stat = "identity", position=position_dodge(0.8), width=0.8) +
  xlab(label = "") +
  ylab(label = "Weighted Ave F1 Score") +
  labs(fill= "Proportions") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size = 10),
        text = element_text(size=20))+
  geom_text(aes(label = round(WeightedAveF1, 3)), position = position_dodge(.8), size = 4, hjust=1) +
  facet_grid(Type~level , scales = "free_x", space = "free_x", labeller = labeller(level = level.labs)) +
  coord_flip() +
  scale_fill_manual(values=c("#9d92d6", "#66b0f0", "#fffd7c", "#c69189", "#ffb169", "#83d0ff", "#046C9A",
                             "#faa9b4", "#FAD510","#EABE94", "#9986A5", "#CCC591", "#F98400", "#c3abdb",
                             "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#E79805", "#EC7A05"),
                    guide=guide_legend(reverse = TRUE),
                    labels = vector.labels)

# heatmap for F1-score
library(hrbrthemes)
library(ggh4x)

metrics <- c("WeightedAveF1")
level.labs <- c("Cell-subtype", "Main cell-type")
names(level.labs) <- c("D", "G")
df <- results[,c(metrics, "method", "proportion", "level")]
df[,metrics] <- as.numeric(df[,metrics])

# session 1
df1 <- df[df$proportion %in% paste0("Prop", 1:5),]
ggplot(df1, aes(proportion, method, fill= WeightedAveF1)) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  facet_grid(level~., labeller = labeller(level = level.labs)) +
  theme_bw() +
  xlab(label = "") +
  ylab(label = "Weighted Ave F1 Score") +
  labs(fill= "Proportions") +
  scale_x_discrete(labels=vector.labels)

# session 2
df2 <- df[df$proportion %in% paste0("Prop", 6:26),]
df2 <- df2 %>% mutate(cellcounts = case_when(proportion %in% paste0("Prop", 6:12) ~ "Same cell counts",
                                             proportion %in% paste0("Prop", 13:19) ~ "1.5 times cell counts",
                                             proportion %in% paste0("Prop", 20:26) ~ "2 times cell counts"))
df2$cellcounts <- factor(df2$cellcounts, levels = c("Same cell counts", "1.5 times cell counts","2 times cell counts"))

ggplot(df2, aes(proportion, method, fill= WeightedAveF1)) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  facet_grid(level~cellcounts, scales = "free", space = "free", labeller = labeller(level = level.labs)) +
  theme_bw() +
  xlab(label = "") +
  ylab(label = "Weighted Ave F1 Score") +
  labs(fill= "Proportions") +
  scale_x_discrete(labels=vector.labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))

# session 3
df3 <- df[df$proportion %in% paste0("Prop", 27:40),] # 1:5, 6:12, 13:26, 27:40
df3 <- df3 %>% mutate(Type = case_when(proportion %in% paste0("Prop", 27:33) ~ "Epithelial-subtypes",
                                     proportion %in% paste0("Prop", 34:40) ~ "Immune-subtypes"))
df3$Type <- factor(df3$Type, levels = c("Epithelial-subtypes", "Immune-subtypes"))
df3$proportion <- factor(df3$proportion, levels = paste0("Prop", c(34:40,27:33)))

ggplot(df3, aes(proportion, method, fill= WeightedAveF1)) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  facet_grid(level~Type, scales = "free", space = "free", labeller = labeller(level = level.labs)) +
  theme_bw() +
  xlab(label = "") +
  ylab(label = "Weighted Ave F1 Score") +
  labs(fill= "Proportions") +
  scale_x_discrete(labels=vector.labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.9))

# accuracy for rare cell type and abundant cell types ######
metrics <- c("WeightedAveF1.rare", "WeightedAveF1.rich")
df <- results[results$level=="D",]
df <- df[,c(metrics, "method", "proportion")]
df$WeightedAveF1.rare <- as.numeric(df$WeightedAveF1.rare)
df$WeightedAveF1.rich <- as.numeric(df$WeightedAveF1.rich)

# session 1
df <- df[df$proportion %in% paste0("Prop", 1:5),] # 1:5, 6:26, 27:40

df <- gather(df, CellType, WeightedAveF1, c(WeightedAveF1.rare,WeightedAveF1.rich), factor_key=TRUE)
df$CellType <- factor(df$CellType, levels = c("WeightedAveF1.rich","WeightedAveF1.rare"))

level.labs <- c("Abundance < 5%", "Abundance > 5%")
names(level.labs) <- c("WeightedAveF1.rare", "WeightedAveF1.rich")

ggplot(data=df, mapping=aes(x=method, y=WeightedAveF1, fill= proportion)) +
  geom_bar(stat = "identity", position=position_dodge(0.8), width=0.8) +
  xlab(label = "") +
  ylab(label = "Weighted Ave F1 Score") +
  labs(fill= "Proportions") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 18),
        text = element_text(size=18),
        axis.text=element_text(size=18),
        axis.text.x = element_text(size = 10)) +
  geom_text(aes(label = round(WeightedAveF1, 3)), position = position_dodge(.8), size = 4, hjust=1) +
  facet_grid(CellType ~. , scales = "free_x", space = "free_x", labeller = labeller(CellType  = level.labs)) +
  # facet_grid(CellType~ly_temp , scales = "free_x", space = "free_x", labeller = labeller(CellType = level.labs)) +
  coord_flip() +
  scale_fill_manual(values=c("#9d92d6", "#66b0f0", "#fffd7c", "#c69189", "#ffb169", "#83d0ff", "#046C9A",
                             "#faa9b4", "#FAD510","#EABE94", "#9986A5", "#CCC591", "#F98400", "#c3abdb",
                             "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#E79805", "#EC7A05"),
                    guide=guide_legend(reverse = TRUE),
                    labels = vector.labels)
# session 2
df <- df[df$proportion %in% paste0("Prop", 6:26),] # 1:5, 6:12 , 13:26, 27:40

df <- df %>% mutate(cellcounts = case_when(proportion %in% paste0("Prop", 6:12) ~ "Same cell counts",
                                           proportion %in% paste0("Prop", 13:19) ~ "1.5 times cell counts",
                                           proportion %in% paste0("Prop", 20:26) ~ "2 times cell counts"))
df$cellcounts <- factor(df$cellcounts, levels = c("Same cell counts", "1.5 times cell counts","2 times cell counts"))

df <- gather(df, CellType, WeightedAveF1, c(WeightedAveF1.rare,WeightedAveF1.rich), factor_key=TRUE)
df$CellType <- factor(df$CellType, levels = c("WeightedAveF1.rich","WeightedAveF1.rare"))

level.labs <- c("Abundance < 5%", "Abundance > 5%")
names(level.labs) <- c("WeightedAveF1.rare", "WeightedAveF1.rich")
ggplot(data=df, mapping=aes(x=method, y=WeightedAveF1, fill= proportion)) +
  geom_bar(stat = "identity", position=position_dodge(0.8), width=0.8) +
  xlab(label = "") +
  ylab(label = "Weighted Ave F1 Score") +
  labs(fill= "Proportions") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 18),
        text = element_text(size=18),
        axis.text=element_text(size=18),
        axis.text.x = element_text(size = 10))+
  geom_text(aes(label = round(WeightedAveF1, 3)), position = position_dodge(.8), size = 4, hjust=1) +
  # facet_grid(CellType ~. , scales = "free_x", space = "free_x", labeller = labeller(CellType  = level.labs)) +
  facet_grid(CellType~cellcounts , scales = "free_x", space = "free_x", labeller = labeller(CellType = level.labs)) +
  coord_flip() +
  scale_fill_manual(values=c("#9d92d6", "#66b0f0", "#fffd7c", "#c69189", "#ffb169", "#83d0ff", "#046C9A",
                             "#faa9b4", "#FAD510","#EABE94", "#9986A5", "#CCC591", "#F98400", "#c3abdb",
                             "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#E79805", "#EC7A05"),
                    guide=guide_legend(reverse = TRUE),
                    labels = vector.labels)

# session 3
df <- df[df$proportion %in% paste0("Prop", 27:40),] # 1:5, 6:12 , 13:26, 27:40

df <- df %>% mutate(Type = case_when(proportion %in% paste0("Prop", 27:33) ~ "Epithelial-subtype",
                                           proportion %in% paste0("Prop", 34:40) ~ "Immune-subtype"))
df$Type <- factor(df$Type, levels = c("Epithelial-subtype", "Immune-subtype"))
df$proportion <- factor(df$proportion, levels = paste0("Prop", c(34:40, 27:33)))

df <- gather(df, CellType, WeightedAveF1, c(WeightedAveF1.rare,WeightedAveF1.rich), factor_key=TRUE)
df$CellType <- factor(df$CellType, levels = c("WeightedAveF1.rich","WeightedAveF1.rare"))
# df$proportion <- factor(df$proportion, levels = paste0("Prop", c(34:40,27:33)))

level.labs <- c("Abundance < 5%", "Abundance > 5%")
names(level.labs) <- c("WeightedAveF1.rare", "WeightedAveF1.rich")
ggplot(data=df, mapping=aes(x=method, y=WeightedAveF1, fill= proportion)) +
  geom_bar(stat = "identity", position=position_dodge(0.8), width=0.8) +
  xlab(label = "") +
  ylab(label = "Weighted Ave F1 Score") +
  labs(fill= "Proportions") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 18),
        text = element_text(size=18),
        axis.text=element_text(size=18),
        axis.text.x = element_text(size = 10))+
  geom_text(aes(label = round(WeightedAveF1, 3)), position = position_dodge(.8), size = 4, hjust=1) +
  # facet_grid(CellType ~. , scales = "free_x", space = "free_x", labeller = labeller(CellType  = level.labs)) +
  facet_grid(Type~CellType, scales = "free_x", space = "free_x", labeller = labeller(CellType = level.labs)) +
  coord_flip() +
  scale_fill_manual(values=c("#9d92d6", "#66b0f0", "#fffd7c", "#c69189", "#ffb169", "#83d0ff", "#046C9A",
                             "#faa9b4", "#FAD510","#EABE94", "#9986A5", "#CCC591", "#F98400", "#c3abdb",
                             "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#E79805", "#EC7A05"),
                    guide=guide_legend(reverse = TRUE),
                    labels = vector.labels)


# Weighted average precision and recall ####
metrics <- c("Precision", "Recall")
df <- results[results$level=="D",]
df <- df[,c(metrics, "method", "proportion")]
df$Precision <- as.numeric(df$Precision)
df$Recall <- as.numeric(df$Recall)

# session 1
df <- df[df$proportion %in% paste0("Prop", 1:5),] # 1:5, 6:26, 27:40

df <- gather(df, CellType, P_R, c(Precision,Recall), factor_key=TRUE)

level.labs <- c("WeigtedAve_Precision", "WeigtedAve_Recall")
names(level.labs) <- c("Precision", "Recall")

ggplot(data=df, mapping=aes(x=method, y=P_R, fill= proportion)) +
  geom_bar(stat = "identity", position=position_dodge(0.8), width=0.8) +
  xlab(label = "") +
  ylab(label = "Precision and Recall") +
  labs(fill= "Proportions") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 18),
        text = element_text(size=18),
        axis.text=element_text(size=18))+
  geom_text(aes(label = round(P_R, 3)), position = position_dodge(.8), size = 4, hjust=1) +
  facet_grid(CellType ~. , scales = "free_x", space = "free_x", labeller = labeller(CellType  = level.labs)) +
  # facet_grid(CellType~ly_temp , scales = "free_x", space = "free_x", labeller = labeller(CellType = level.labs)) +
  coord_flip() +
  scale_fill_manual(values=c("#9d92d6", "#66b0f0", "#fffd7c", "#c69189", "#ffb169", "#83d0ff", "#046C9A",
                             "#faa9b4", "#FAD510","#EABE94", "#9986A5", "#CCC591", "#F98400", "#c3abdb",
                             "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#E79805", "#EC7A05"),
                    guide=guide_legend(reverse = TRUE),
                    labels = vector.labels)

# session 2
df <- df[df$proportion %in% paste0("Prop", 6:26),] # 1:5, 6:12, 13:26, 27:40
df <- df %>% mutate(cellcounts = case_when(proportion %in% paste0("Prop", 6:12) ~ "Same cell counts",
                                           proportion %in% paste0("Prop", 13:19) ~ "1.5 times cell counts",
                                           proportion %in% paste0("Prop", 20:26) ~ "2 times cell counts"))
df$cellcounts <- factor(df$cellcounts, levels = c("Same cell counts", "1.5 times cell counts","2 times cell counts"))
df <- gather(df, CellType, P_R, c(Precision,Recall), factor_key=TRUE)

level.labs <- c("WeigtedAve_Precision", "WeigtedAve_Recall")
names(level.labs) <- c("Precision", "Recall")

ggplot(data=df, mapping=aes(x=method, y=P_R, fill= proportion)) +
  geom_bar(stat = "identity", position=position_dodge(0.8), width=0.8) +
  xlab(label = "") +
  ylab(label = "Precision and Recall") +
  labs(fill= "Proportions") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 18),
        text = element_text(size=18),
        axis.text=element_text(size=18))+
  geom_text(aes(label = round(P_R, 3)), position = position_dodge(.8), size = 4, hjust=1) +
  facet_grid(CellType ~cellcounts, scales = "free_x", space = "free_x", labeller = labeller(CellType  = level.labs)) +
  # facet_grid(CellType~ly_temp , scales = "free_x", space = "free_x", labeller = labeller(CellType = level.labs)) +
  coord_flip() +
  scale_fill_manual(values=c("#9d92d6", "#66b0f0", "#fffd7c", "#c69189", "#ffb169", "#83d0ff", "#046C9A",
                             "#faa9b4", "#FAD510","#EABE94", "#9986A5", "#CCC591", "#F98400", "#c3abdb",
                             "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#E79805", "#EC7A05"),
                    guide=guide_legend(reverse = TRUE),
                    labels = vector.labels)

# session 3
df <- df[df$proportion %in% paste0("Prop", 27:40),] # 1:5, 6:26, 27:40
df <- df %>% mutate(Type = case_when(proportion %in% paste0("Prop", 27:33) ~ "Epithelial-subtype",
                                     proportion %in% paste0("Prop", 34:40) ~ "Immune-subtype"))
df$Type <- factor(df$Type, levels = c("Epithelial-subtype", "Immune-subtype"))
df$proportion <- factor(df$proportion, levels = paste0("Prop", c(34:40, 27:33)))
df <- gather(df, CellType, P_R, c(Precision,Recall), factor_key=TRUE)

level.labs <- c("WeigtedAve_Precision", "WeigtedAve_Recall")
names(level.labs) <- c("Precision", "Recall")

ggplot(data=df, mapping=aes(x=method, y=P_R, fill= proportion)) +
  geom_bar(stat = "identity", position=position_dodge(0.8), width=0.8) +
  xlab(label = "") +
  ylab(label = "Precision and Recall") +
  labs(fill= "Proportions") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 18),
        text = element_text(size=18),
        axis.text=element_text(size=18),
        axis.text.x = element_text(size=15))+
  geom_text(aes(label = round(P_R, 3)), position = position_dodge(.8), size = 4, hjust=1) +
  facet_grid(Type~CellType, scales = "free_x", space = "free_x", labeller = labeller(CellType  = level.labs)) +
  # facet_grid(CellType~ly_temp , scales = "free_x", space = "free_x", labeller = labeller(CellType = level.labs)) +
  coord_flip() +
  scale_fill_manual(values=c("#9d92d6", "#66b0f0", "#fffd7c", "#c69189", "#ffb169", "#83d0ff", "#046C9A",
                             "#faa9b4", "#FAD510","#EABE94", "#9986A5", "#CCC591", "#F98400", "#c3abdb",
                             "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#E79805", "#EC7A05"),
                    guide=guide_legend(reverse = TRUE),
                    labels = vector.labels)

# heatmap for rare & rich cell precision and recall: session 2 #####
library(ggplot2)
library(viridis)
library(grid)
library(gtable)

method <- c("CHETAH", "scmapcluster", "SingleCellNet", "SingleR", "scClassify", "scPred")

vector.labels <- c("QR(10:1)", "QR(5:1)", "QR(2:1)", "QR(1:1)", "QR(0.5:1)",
                   "R(10:1)", "R(5:1)","R(2:1)", "R(1:1)", "R(0.5:1)", "R(0.2:1)", "R(0.1:1)",
                   "R(10:1)*1.5", "R(5:1)*1.5","R(2:1)*1.5", "R(1:1)*1.5", "R(0.5:1)*1.5", "R(0.2:1)*1.5", "R(0.1:1)*1.5",
                   "R(10:1)*2", "R(5:1)*2","R(2:1)*2", "R(1:1)*2", "R(0.5:1)*2", "R(0.2:1)*2", "R(0.1:1)*2",
                   "LC2:LC1(4:1)", "LC2:LC1(2.5:1)", "LC2:LC1(2:1)", "LC2:LC1(1:1)", "LC2:LC1(0.5:1)", "LC2:LC1(0.4:1)", "LC2:LC1(0.25:1)",
                   "Mac:T(4:1)", "Mac:T(2.5:1)", "Mac:T(2:1)", "Mac:T(1:1)", "Mac:T(0.5:1)", "Mac:T(0.4:1)", "Mac:T(0.25:1)")
names(vector.labels) <- paste0("Prop", 1:40)
combine_everything_p_r <- function(methods){
  results.list <- list()
  for(method in methods){
    path.1 <- paste0("08_StartOver_NewStory/03_evaluation/", method)
    for (prp in 1:length(list.files(path.1))) {
      path.2 <- paste0(path.1, "/P", prp, "/")

      df.1 <- read.csv(paste0(path.2, "evaluating_result_K5.csv"), row.names = 1)

      df.1 <- df.1[,c("Precision.rare", "Precision.rich", "Recall.rare", "Recall.rich")]

      vec.1 <- colMeans(df.1)

      vec.1 <- append(vec.1, c(method=method, proportion=paste0("Prop", prp), level="D"))

      results.list <- rlist::list.append(results.list, vec.1)
    }
  }
  return(results.list)
}
ref.list <- combine_everything_p_r(methods = method)
results <- data.frame(t(sapply(ref.list, c)))
results <- results[results$proportion %in% paste0("Prop", 6:26),] # 1:5, 6:12, 13:26, 27:40
results$proportion <- factor(results$proportion, levels= paste0("Prop",c(6:12, 13:26)))
results <- results %>%
  mutate(method = replace(method, method == "scmapcluster", "scmap-cluster"))
df <- results
df <- df %>% mutate(cellcounts = case_when(proportion %in% paste0("Prop", 6:12) ~ "Same cell counts",
                                           proportion %in% paste0("Prop", 13:19) ~ "1.5 times cell counts",
                                           proportion %in% paste0("Prop", 20:26) ~ "2 times cell counts"))
df$cellcounts <- factor(df$cellcounts, levels = c("Same cell counts", "1.5 times cell counts","2 times cell counts"))
df <- gather(df, Metrics, Score, c(Precision.rare, Precision.rich, Recall.rare, Recall.rich), factor_key=TRUE)
df$Score <- as.numeric(df$Score)

vector.labels.hm <- paste0("P",c(1:7, 1:7, 1:7))
names(vector.labels.hm) <- paste0("Prop", c(6:26))

ggplot(df, aes(x=proportion, y=Metrics, fill= Score)) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  facet_grid2(cellcounts ~ method, scales = "free", independent = "x") +
  # facet_wrap(cellcounts ~ method, drop=TRUE, ncol = 6, scales = "free_x") +
  theme_bw() +
  theme(text=element_text(size=15)) +
  scale_x_discrete(labels=vector.labels.hm) +
  scale_y_discrete(labels=c("Precision (minority)", "Precision (majority)", "Recall (minority)", "Recall (majority)"))


# robustness (scatter) ####
metrics <- c("WeightedAveF1")

level.labs <- c("Cell-subtype", "Main cell-type")
names(level.labs) <- c("D", "G")

df <- results[,c(metrics, "method", "proportion", "level")]
df[,metrics] <- as.numeric(df[,metrics])

# session 1
df <- df[df$proportion %in% paste0("Prop", 1:5),] # 1:5
dist.ref <- df %>%
  dplyr::group_by(method, level) %>%
  dplyr::summarise_at(vars(-c(proportion)), list(mean = mean, sd = sd))
# session 2
df <- df[df$proportion %in% paste0("Prop", 6:26),] # 1:5, 6:26, 27:40
df <- df %>% mutate(cellcounts = case_when(proportion %in% paste0("Prop", 6:12) ~ "Same cell counts",
                                           proportion %in% paste0("Prop", 13:19) ~ "1.5 times cell counts",
                                           proportion %in% paste0("Prop", 20:26) ~ "2 times cell counts"))
df$cellcounts <- factor(df$cellcounts, levels = c("Same cell counts", "1.5 times cell counts","2 times cell counts"))
dist.ref <- df %>%
  dplyr::group_by(method, level, cellcounts) %>%
  dplyr::summarise_at(vars(-c(proportion)), list(mean = mean, sd = sd))

# session 3
df <- df[df$proportion %in% paste0("Prop", 27:40),] # 1:5, 6:26, 27:40
df <- df %>% mutate(Type = case_when(proportion %in% paste0("Prop", 27:33) ~ "Epithelial-subtypes",
                                     proportion %in% paste0("Prop", 34:40) ~ "Immune-subtypes"))
df$Type <- factor(df$Type, levels = c("Epithelial-subtypes", "Immune-subtypes"))
dist.ref <- df %>%
  dplyr::group_by(method, level, Type) %>%
  dplyr::summarise_at(vars(-c(proportion)), list(mean = mean, sd = sd))

# scatter plot
ggplot(dist.ref, aes(x=mean, y=sd)) +
  labs(x= "Mean of F1-score", y="Standard Deviation of F1-score")+
  geom_point(aes(colour=method, shape = Type), size=3.5) + # Scenario 3
  # geom_point(aes(colour=method, shape = cellcounts), size=3.5) + # Scenario 2
  # geom_point(aes(colour=method), size=3.5) + # Scenario 1
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        # legend.position = "none",
        strip.text.x = element_text(size = 20)) +
  facet_wrap(level~., labeller = labeller(level = level.labs)) +
  geom_text_repel(aes(label = method)) +
  theme(legend.text = element_text(size = 12))





# for assigned and unassigned rate #####
assignment <- results[,c("assignedRate", "method", "proportion", "level")]
assignment$assignment <- "Assigned"
colnames(assignment)[1] <- "Rate"
assignment$Rate <- as.numeric(assignment$Rate)
assignment <- assignment[assignment$level=="D",]
assignment <- subset(assignment, select = -c(level, assignment))

assignment$ly[assignment$proportion %in% paste0("Prop", 1:5)] <- "Scenario 1"
assignment$ly[assignment$proportion %in% paste0("Prop", 6:26)] <- "Scenario 2"
assignment$ly[assignment$proportion %in% paste0("Prop", 27:40)] <- "Scenario 3"

ggplot(assignment, aes(y=Rate, x=proportion)) +
  geom_bar(stat="identity", fill = "#EAD3BF") +
  facet_grid2(vars(ly), vars(method), scales = "free_x", independent = "x") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 20),
        text = element_text(size=20),
        legend.position = "none",
        axis.text=element_text(size=12),
        axis.ticks.y=element_blank(),
        panel.spacing = unit(0.01, "cm")) +
  ylab("Assigned rate (%)") +
  xlab("") +
  scale_x_discrete(labels = vector.labels,
                   guide = guide_axis(angle = 90))



# for runtime #####
Time <- results[,c("Time", "method", "proportion", "level")]
Time <- Time[Time$level=="D",]
Time$Time <- as.numeric(Time$Time)

Time$ly[Time$proportion %in% paste0("Prop", 1:5)] <- "Scenario 1"
Time$ly[Time$proportion %in% paste0("Prop", 6:26)] <- "Scenario 2"
Time$ly[Time$proportion %in% paste0("Prop", 27:40)] <- "Scenario 3"


ggplot(data=Time, aes(x=proportion, y=Time, group=method, color=method)) +
  geom_line(size = 0.8) +
  geom_point(size = 3) +
  facet_grid(~ly, scales = "free_x", space="free_x") +
  labs(y="Runtime (min)", x="", color="Methods") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 15),
        text = element_text(size=20),
        axis.text=element_text(size=12),
        axis.ticks.y=element_blank(),
        panel.spacing = unit(0.1, "cm"),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank()) +
  scale_color_manual(values = c("#9d92d6", "#faa9b4", "#fffd7c", "#c69189", "#ffb169", "#83d0ff")) +
  scale_x_discrete(labels = vector.labels,
                   guide = guide_axis(angle = 45))



# cell counts and proportions (Session 1, 2, 3, 4) #####

ind.list <- lapply(c(1:40), function(x){
  readRDS(paste0("08_StartOver_NewStory/01_datasets/index/Index_P", x, ".rds"))
})

cR <- readRDS("08_StartOver_NewStory/01_datasets/seuobj/combined.p1.rds")
cR <- cR@meta.data$bm_general
c.list <- lapply(1:40, function(x){cR})
rm(cR)

c.test <- t(sapply(lapply(seq_along(c.list), function(i){
  sapply(seq_along(ind.list[[i]]), function(j){
    x <- c.list[[i]][ind.list[[i]][[j]][[1]]]
    x <- t(table(x))
  })
}), rowMeans))

c.ref <- t(sapply(lapply(seq_along(c.list), function(i){
  sapply(seq_along(ind.list[[i]]), function(j){
    x <- c.list[[i]][ind.list[[i]][[j]][[2]]]
    x <- t(table(x))
  })
}), rowMeans))

c.final <- list(c.test, c.ref)
c.final <- lapply(seq_along(c.final), function(i){
  c.final[[i]] <- data.frame(c.final[[i]])
  colnames(c.final[[i]]) <- c("Epithelial", "Immune")
  c.final[[i]]$Proportion <- factor(paste0("Prop", c(1:40)), levels = paste0("Prop", c(1:40)))
  c.final[[i]] <- gather(c.final[[i]], celltype, cellcounts, c(Epithelial, Immune), factor_key = TRUE)
})

c.final[[1]]$group <- "Query"
c.final[[2]]$group <- "Reference"

c.final <- Reduce(rbind, c.final)
rm(c.ref, c.test)
c.final$group <- factor(c.final$group, levels = c("Query","Reference"))

c.final.S1 <- c.final[c.final$Proportion %in% paste0("Prop", 1:5),]

c.final.S1 <- c.final.S1 %>%
  dplyr::group_by(Proportion, group) %>%
  dplyr::mutate(Value=cellcounts/sum(cellcounts))



c.final.S2 <- c.final[c.final$Proportion %in% paste0("Prop", 6:26),]
c.final.S2 <- c.final.S2 %>%
  dplyr::group_by(Proportion, group) %>%
  dplyr::mutate(Value=cellcounts/sum(cellcounts))
c.final.S2 <- c.final.S2 %>% mutate(grouping = case_when(Proportion %in% paste0("Prop", 6:12) ~ "Same cell counts",
                                                             Proportion %in% paste0("Prop", 13:19) ~ "1.5 times cell counts",
                                                             Proportion %in% paste0("Prop", 20:26) ~ "2 times cell counts"))
c.final.S2$grouping <- factor(c.final.S2$grouping, levels = c("Same cell counts", "1.5 times cell counts","2 times cell counts"))




c.final.S3 <- c.final[c.final$Proportion %in% paste0("Prop", 27:40),]
c.final.S3 <- c.final.S3 %>%
  dplyr::group_by(Proportion, group) %>%
  dplyr::mutate(Value=cellcounts/sum(cellcounts))




library(wesanderson)
library(gridExtra)
library(grid)


p1 <- ggplot(data=c.final.S1, aes(x=Proportion, y=Value, fill=celltype)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(group~., ncol = 1, nrow = 2) +
  geom_text(aes(label=round(cellcounts, digits = 0)), position = position_stack(vjust = .5), size=3.5, angle=0) +
  scale_fill_manual(values=c("#F2AD00", "#5BBCD6")) +
  xlab("") +
  ylab("Cell proportions") +
  labs(fill = "Cell Type") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 18),
        text = element_text(size=18),
        # legend.position = "none",
        axis.text=element_text(size=13),
        axis.text.x=element_text(size=12, angle=0, hjust=1),
        axis.ticks.y=element_blank()) +
  scale_x_discrete(labels=vector.labels)

library(ggh4x)


p2 <- ggplot(data=c.final.S2, aes(x=Proportion, y=Value, fill=celltype)) +
  geom_bar(position="stack", stat="identity") +
  facet_nested(~ group + grouping, scales = "free_x") +
  # facet_wrap(group~grouping, ncol = 3, nrow = 2, scales = "free_x") +
  geom_text(aes(label=round(cellcounts, digits = 0)), position = position_stack(vjust = .5), size=3.5, angle=0) +
  scale_fill_manual(values=c("#F2AD00", "#5BBCD6")) +
  xlab("") +
  ylab("Cell proportions") +
  labs(fill = "Cell Type") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 20),
        text = element_text(size=20),
        # legend.position = "none",
        axis.text=element_text(size=12),
        axis.text.x=element_text(size=12, angle=45, hjust=1),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank()) +
  scale_x_discrete(labels=vector.labels)


p3 <- ggplot(data=c.final.S3, aes(x=Proportion, y=Value, fill=celltype)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(group~., ncol = 1, nrow = 2) +
  geom_text(aes(label=round(cellcounts, digits = 0)), position = position_stack(vjust = .5), size=3.5, angle=0) +
  scale_fill_manual(values=c("#F2AD00", "#5BBCD6")) +
  xlab("") +
  ylab("Cell proportions") +
  labs(fill = "Cell Type") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 20),
        text = element_text(size=20),
        # legend.position = "none",
        axis.text=element_text(size=12),
        axis.text.x=element_text(size=15, angle=90, hjust=1),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank()) +
  scale_x_discrete(labels=vector.labels)



# "#78B7C5", "#3B9AB2", "#F3DF6C", "#CEAB07"
# "#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00"

grid.arrange(p1, p2, p3, nrow=1, widths=c(1.2, 1.5, 3), top=textGrob(
  "The Counts and Proportions of Cells Used in Four Sessions of Comparison",
  gp = gpar(fontface = 1, fontsize = 25),
  hjust = -0.5,
  vjust = 0.5,
  x = 0
))

# session 2


# sub-types

ind.list.R3 <- lapply(27:40, function(x){
  readRDS(paste0("08_StartOver_NewStory/01_datasets/index/Index_P",x,".rds"))
})

cR2 <- readRDS("08_StartOver_NewStory/01_datasets/seuobj/combined.p1.rds")
cR2 <- cR2@meta.data$bm_cell_type
c.list.R3 <- lapply(1:14, function(x){cR2})

c.test.R3 <- t(sapply(lapply(seq_along(c.list.R3), function(i){
  sapply(seq_along(ind.list.R3[[i]]), function(j){
    x <- c.list.R3[[i]][ind.list.R3[[i]][[j]][[2]]]
    x <- t(table(factor(x, levels = c("B_cell", "Basal", "DC", "G2MS", "LC1", "LC2", "Macrophage", "Neutrophil", "NK", "T_cell"))))
  })
}), rowMeans))

colnames(c.test.R3) <- c("B_cell", "Basal", "DC", "Cycling", "LC1", "LC2", "Macrophage", "Neutrophil", "NK", "T_cell")
rownames(c.test.R3) <- paste0("Prop", c(27:40))
rm(c.list.R3, ind.list.R3)

c.test.R3 <- c.test.R3[,c("LC1", "LC2", "Macrophage", "T_cell")]

c.test.R3 <- data.frame(c.test.R3)


c.test.R3$Proportion <- rownames(c.test.R3)
c.test.R3 <- gather(c.test.R3, celltype, cellcounts, c(LC1, LC2, Macrophage, T_cell), factor_key = TRUE)

c.test.R3$Proportion <- factor(c.test.R3$Proportion, levels=paste0("Prop", c(27:40)))

c.test.R3$celltype <- factor(c.test.R3$celltype, levels = c("LC2", "LC1", "Macrophage", "T_cell"))
c.test.R3 <- c.test.R3 %>% mutate(group=case_when(
  celltype %in% c("LC2", "LC1") ~ "Epithelial",
  celltype %in% c("Macrophage", "T_cell") ~ "Immune"
))

c.test.R3 <- c.test.R3 %>%
  dplyr::group_by(Proportion, group) %>%
  dplyr::mutate(Value=cellcounts/sum(cellcounts))
ggplot(data=c.test.R3, aes(x=Proportion, y=Value, fill=celltype)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(group~.,nrow=2, strip.position="right") +
  geom_text(aes(label=round(cellcounts, digits = 0)), position = position_stack(vjust = .5), size=3.5, angle=0) +
  scale_fill_manual(values=c("#E1AF00", "#EBCC2A", "#3B9AB2", "#78B7C5")) +
  xlab("") +
  ylab("Cell proportions") +
  labs(fill = "Cell Sub-type") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 20),
        text = element_text(size=20),
        axis.text=element_text(size=12),
        axis.ticks.y=element_blank(),
        panel.spacing = unit(0.01, "cm"),
        axis.text.x = element_text(size=15, angle = 90, hjust=1)) +
  scale_x_discrete(labels=vector.labels)



# cell counts and proportions (S4, specifically showing the subtype proportions) #####

ind.list.S4 <- lapply(27:40, function(x){
  readRDS(paste0("08_StartOver_NewStory/01_datasets/index/Index_P",x,".rds"))
})

cR <- readRDS("01_dataset/combined.p1.rds")
cR <- cR@meta.data$bm_cell_type
c.list.S4 <- lapply(1:14, function(x){cR})

c.ref.S4 <- t(sapply(lapply(seq_along(c.list.S4), function(i){
  sapply(seq_along(ind.list.S4[[i]]), function(j){
    x <- c.list.S4[[i]][ind.list.S4[[i]][[j]][[2]]]
    x <- t(table(factor(x, levels = c("B_cell", "Basal", "DC", "G2MS", "LC1", "LC2", "Macrophage", "Neutrophil", "NK", "T_cell"))))
  })
}), rowMeans))

colnames(c.ref.S4) <- c("B_cell", "Basal", "DC", "Cycling", "LC1", "LC2", "Macrophage", "Neutrophil", "NK", "T_cell")
rownames(c.ref.S4) <- paste0("Prop", c(27:40))
rm(c.list.S4, ind.list.S4)

c.ref.S4 <- c.ref.S4[,c("LC1", "LC2", "Macrophage", "T_cell")]

c.ref.S4 <- data.frame(c.ref.S4)


c.ref.S4$Proportion <- rownames(c.ref.S4)
c.ref.S4 <- gather(c.ref.S4, celltype, cellcounts, c(LC1, LC2, Macrophage, T_cell), factor_key = TRUE)

c.ref.S4$Proportion <- factor(c.ref.S4$Proportion, levels=paste0("Prop", c(27:40)))

c.ref.S4$celltype <- factor(c.ref.S4$celltype, levels = c("LC2", "LC1", "Macrophage", "T_cell"))
c.ref.S4 <- c.ref.S4 %>% mutate(group=case_when(
  celltype %in% c("LC2", "LC1") ~ "Epithelial",
  celltype %in% c("Macrophage", "T_cell") ~ "Immune"
))

c.ref.S4 <- c.ref.S4 %>%
  dplyr::group_by(Proportion, group) %>%
  dplyr::mutate(Value=cellcounts/sum(cellcounts))


ggplot(data=c.ref.S4, aes(x=Proportion, y=Value, fill=celltype)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(group~.,nrow=2, strip.position="right") +
  geom_text(aes(label=round(cellcounts, digits = 0)), position = position_stack(vjust = .5), size=3.5, angle=0) +
  scale_fill_manual(values=c("#E1AF00", "#EBCC2A", "#3B9AB2", "#78B7C5")) +
  xlab("") +
  ylab("Cell proportions") +
  labs(fill = "Cell Sub-type") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 20),
        text = element_text(size=20),
        axis.text.x=element_text(size=12, angle = 45, hjust = 1),
        axis.text=element_text(size=12),
        axis.ticks.y=element_blank(),
        panel.spacing = unit(0.01, "cm")) +
  scale_x_discrete(labels=vector.labels)


# complexity visualization #####
complexity <- read.csv(file = "03_evaluating/complexity_abdelaal.csv", row.names = 1)
complexity <- gather(complexity, Type, complexity_score, c(Query, Ref), factor_key=TRUE)

# strip plot
complexity$Proportion <- factor(complexity$Proportion, levels = paste0("Prop", c(1:4, 5:14, 25:34)))

levels(complexity$Type) <- c(levels(complexity$Type), "Reference")
complexity$Type[complexity$Type == 'Ref'] <- 'Reference'
levels(complexity$Type)[match("Ref",levels(complexity$Type))] <- "Reference"


complexity <- complexity %>% mutate(ly = case_when(
  Proportion %in% paste0("Prop", 1:4) ~ "Scenario 1",
  Proportion %in% paste0("Prop", 5:14) ~ "Scenario 2",
  Proportion %in% paste0("Prop", 25:34) ~ "Scenario 3"
))

complexity.1 <- complexity[complexity$ly == "Scenario 1", ]
complexity.2 <- complexity[complexity$ly == "Scenario 2", ]
complexity.3 <- complexity[complexity$ly == "Scenario 3", ]


ggplot(complexity.2, aes(x=Proportion, y=complexity_score, color=Type)) + # complexity.3, complexity.2, complexity.1
  geom_jitter(position=position_dodge(0.01)) +
  facet_grid(.~ly, scales = "free_x", space="free_x") +
  theme(strip.text.x = element_text(size = 20),
        text = element_text(size=20),
        axis.text=element_text(size=12),
        axis.ticks.y=element_blank()) +
  xlab("") +
  ylab("Complexity score") +
  labs(color=NULL) +
  theme_bw() +
  scale_color_manual(values=c(Query="#FD6467", Reference="#5B1A18")) +
  scale_x_discrete(labels=c("Prop1" = "QR(10:1)", "Prop2" = "QR(5:1)", "Prop3" = "QR(1:1)", "Prop4" = "QR(0.5:1)",
                            "Prop5" = "R(10:1)", "Prop6" = "R(5:1)", "Prop7" = "R(2:1)", "Prop8" = "R(1:1)", "Prop9" = "R(0.5:1)",
                            "Prop10"= "Q(10:1)", "Prop11" = "Q(5:1)", "Prop12" = "Q(2:1)", "Prop13" = "Q(1:1)", "Prop14" = "Q(0.5:1)",
                            "Prop25" = "LC2(3:1)", "Prop26" = "LC2(2:1)", "Prop27" = "LC2(1:1)", "Prop28" = "LC2(0.5:1)", "Prop29" = "LC2(0.1:1)",
                            "Prop30"= "M(3:1)", "Prop31" = "M(2:1)", "Prop32" = "M(1:1)", "Prop33" = "M(0.5:1)", "Prop34" = "M(0.1:1)"))


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(median = median(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      max = max(x[[col]], na.rm=TRUE),
      min = min(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("median" = varname))
  return(data_sum)
}

df1 <- data_summary(complexity, varname="complexity_score",
                    groupnames=c("Proportion", "Type"))
df1$ly[df1$Proportion %in% paste0("Prop", 1:4)] <- "Scenario 1"
df1$ly[df1$Proportion %in% paste0("Prop", 5:14)] <- "Scenario 2"
df1$ly[df1$Proportion %in% paste0("Prop", 25:34)] <- "Scenario 3"

df1.1 <- df1[df1$ly=="Scenario 1", ]
df1.2 <- df1[df1$ly=="Scenario 2", ]
df1.3 <- df1[df1$ly=="Scenario 3", ]

ggplot(df1.3, aes(x=Proportion, y=complexity_score,  color=Type, group = ly))+
  geom_line(aes(group = Type), linewidth=1) +
  geom_point()+
  geom_errorbar(aes(ymin = min, ymax = max), width = 0.2) +
  # facet_grid(.~ly, scales = "free_x", space="free_x") +
  labs(x="", y = "Complexity score", color=NULL)+
  theme_bw() +
  theme(strip.text.x = element_text(size = 20),
        text = element_text(size=20),
        axis.text=element_text(size=10),
        axis.ticks.y=element_blank()) +
  scale_color_manual(values = c("#FD6467", "#5B1A18"))+
  scale_x_discrete(labels=c("Prop1" = "QR(10:1)", "Prop2" = "QR(5:1)", "Prop3" = "QR(1:1)", "Prop4" = "QR(0.5:1)",
                            "Prop5" = "R(10:1)", "Prop6" = "R(5:1)", "Prop7" = "R(2:1)", "Prop8" = "R(1:1)", "Prop9" = "R(0.5:1)",
                            "Prop10"= "Q(10:1)", "Prop11" = "Q(5:1)", "Prop12" = "Q(2:1)", "Prop13" = "Q(1:1)", "Prop14" = "Q(0.5:1)",
                            "Prop25" = "LC2(3:1)", "Prop26" = "LC2(2:1)", "Prop27" = "LC2(1:1)", "Prop28" = "LC2(0.5:1)", "Prop29" = "LC2(0.1:1)",
                            "Prop30"= "M(3:1)", "Prop31" = "M(2:1)", "Prop32" = "M(1:1)", "Prop33" = "M(0.5:1)", "Prop34" = "M(0.1:1)"))

# sparsity #######
sparsity <- read.csv(file = "03_evaluating/sparsity.csv", row.names = 1)
sparsity <- gather(sparsity, Type, sparsity_score, c(Query, Ref), factor_key=TRUE)

# strip plot
sparsity$Proportion <- factor(sparsity$Proportion, levels = paste0("Prop", c(1:4, 5:14, 25:34)))

levels(sparsity$Type) <- c(levels(sparsity$Type), "Reference")
sparsity$Type[sparsity$Type == 'Ref'] <- 'Reference'
levels(sparsity$Type)[match("Ref",levels(sparsity$Type))] <- "Reference"

sparsity <- sparsity %>% mutate(ly = case_when(
  Proportion %in% paste0("Prop", 1:4) ~ "Scenario 1",
  Proportion %in% paste0("Prop", 5:14) ~ "Scenario 2",
  Proportion %in% paste0("Prop", 25:34) ~ "Scenario 3"
))

ggplot(sparsity, aes(x=Proportion, y=sparsity_score, color=Type)) +
  geom_jitter(position=position_dodge(0.01)) +
  facet_grid(.~ly, scales = "free_x", space="free_x") +
  theme(strip.text.x = element_text(size = 20),
        text = element_text(size=20),
        axis.text=element_text(size=12),
        axis.ticks.y=element_blank()) +
  xlab("") +
  ylab("Sparsity") +
  labs(color=NULL) +
  theme_bw() +
  scale_color_manual(values=c("#FD6467", "#5B1A18")) +
  scale_x_discrete(labels=c("Prop1" = "QR(10:1)", "Prop2" = "QR(5:1)", "Prop3" = "QR(1:1)", "Prop4" = "QR(0.5:1)",
                            "Prop5" = "R(10:1)", "Prop6" = "R(5:1)", "Prop7" = "R(2:1)", "Prop8" = "R(1:1)", "Prop9" = "R(0.5:1)",
                            "Prop10"= "Q(10:1)", "Prop11" = "Q(5:1)", "Prop12" = "Q(2:1)", "Prop13" = "Q(1:1)", "Prop14" = "Q(0.5:1)",
                            "Prop25" = "LC2(3:1)", "Prop26" = "LC2(2:1)", "Prop27" = "LC2(1:1)", "Prop28" = "LC2(0.5:1)", "Prop29" = "LC2(0.1:1)",
                            "Prop30"= "M(3:1)", "Prop31" = "M(2:1)", "Prop32" = "M(1:1)", "Prop33" = "M(0.5:1)", "Prop34" = "M(0.1:1)"))

# robustness #####
combine_everything <- function(methods){
  results.list <- list()
  for(method in methods){
    path.1 <- paste0("03_evaluating/", method)
    for (prp in 1:length(list.files(path.1))) {
      path.2 <- paste0(path.1, "/prop", prp, "/")
      if(all(methods==method.marker)){
        for (k in 1:length(list.files(path.2))) {
          tops <- c("top10", "top25", "top50")
          if(length(list.files(path.2))==5){
            tops <- c("top5", "top10", "top25", "top50", "top75")
          }
          path.3 <- paste0(path.2, tops[k], "/")
          Markers <- tops[k]

          df.1 <- read.csv(paste0(path.3, "evaluating_result_K5.csv"), row.names = 1)
          df.2 <- read.csv(paste0(path.3, "evaluating_result_K5_General.csv"), row.names = 1)

          vec.1 <- colMeans(df.1[,c(5,7:10,12:15)])
          vec.2 <- colMeans(df.2[,c(5,7:10,12:15)])

          vec.1 <- append(vec.1, c(method=method, proportion=paste0("Prop", prp), level="D", markers=Markers))
          vec.2 <- append(vec.2, c(method=method, proportion=paste0("Prop", prp), level="G", markers=Markers))

          results.list <- rlist::list.append(results.list, vec.1)
          results.list <- rlist::list.append(results.list, vec.2)
        }
      }
      else{
        path.3 <- path.2
        Markers <- "Ref"

        df.1 <- read.csv(paste0(path.3, "evaluating_result_K5.csv"), row.names = 1)
        df.2 <- read.csv(paste0(path.3, "evaluating_result_K5_General.csv"), row.names = 1)

        vec.1 <- colMeans(df.1[,c(5,7:10,12:15)])
        vec.2 <- colMeans(df.2[,c(5,7:10,12:15)])

        vec.1 <- append(vec.1, c(method=method, proportion=paste0("Prop", prp), level="D", markers=Markers))
        vec.2 <- append(vec.2, c(method=method, proportion=paste0("Prop", prp), level="G", markers=Markers))

        results.list <- rlist::list.append(results.list, vec.1)
        results.list <- rlist::list.append(results.list, vec.2)
      }
    }
  }
  return(results.list)
}

ref.list <- combine_everything(methods = method.ref)
marker.list <- combine_everything(methods = method.marker)

results <- rbind(data.frame(t(sapply(ref.list, c))), data.frame(t(sapply(marker.list, c))))
results <- results[!results$proportion %in% paste0("Prop",15:24),]

results$markers <- factor(results$markers, levels= c("Ref", "top5", "top10", "top25", "top50", "top75"))
results$proportion <- factor(results$proportion, levels= c(paste0("Prop", c(1:14,25:34))))

results <- results[, c("WeightedAveF1", "method", "proportion", "level", "markers")]

results <- results %>%
  mutate(method = replace(method, method == "singleCellNet", "SingleCellNet")) %>%
  mutate(method = replace(method, method == "scmapcell", "scmap-cell"))

results$WeightedAveF1 <- as.numeric(results[,"WeightedAveF1"])
results <- results %>% mutate(ly = case_when(
  proportion %in% paste0("Prop",1:4) ~ "Scenario 1",
  proportion %in% paste0("Prop",10:14) ~ "Scenario 2.1",
  proportion %in% paste0("Prop",5:9) ~ "Scenario 2.2",
  proportion %in% paste0("Prop",25:29) ~ "Scenario 3.1",
  proportion %in% paste0("Prop",30:34) ~ "Scenario 3.2"))

results$distance <- 1-results$WeightedAveF1
results <- results[,-1]

results.ref <- results[results$markers=="Ref",]
results.marker <- results[results$markers!="Ref",]

dist.ref <- results.ref %>%
  dplyr::group_by(method, ly, level) %>%
  dplyr::summarise_at(vars(-c(proportion, markers)), list(mean = mean, SD = sd))

dist.marker <- results.marker %>%
  dplyr::group_by(method, ly, level, markers) %>%
  dplyr::summarise_at(vars(-c(proportion)), list(mean = mean, SD = sd))

# scatter plot (not in use)
ggplot(dist.ref, aes(x=mean, y=SD)) +
  geom_point(aes(colour=method)) +
  facet_wrap(level~ly)

# barplot
dist.ref <- dist.ref[dist.ref$level=="D",]
ggplot(dist.ref[dist.ref$ly=="Scenario 3.2",], aes(x=reorder(method, -mean), y=mean, fill=method)) +
  geom_bar(stat="identity") +
  # facet_wrap(~ly) +
  coord_flip() +
  theme_bw() +
  theme(strip.text.x = element_text(size = 18),
        text = element_text(size=18),
        axis.text=element_text(size=18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")+
  ylab("Rank of Accuracy (1 - F1)") +
  xlab(NULL)

ggplot(dist.ref[dist.ref$ly=="Scenario 3.2",], aes(x=reorder(method, -SD), y=SD, fill=method)) +
  geom_bar(stat="identity") +
  # facet_wrap(~ly) +
  coord_flip() +
  theme_bw() +
  theme(strip.text.x = element_text(size = 18),
        text = element_text(size=18),
        axis.text=element_text(size=18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")+
  ylab("Rank of Robustness (SD)") +
  xlab(NULL)


library(tidyverse)
library(tidytext)

dist.marker <- dist.marker[dist.marker$level=="D",]
ggplot(dist.marker[dist.marker$ly=="Scenario 3.2",], aes(x=reorder_within(method, -mean, markers), y=mean, fill=method)) +
  geom_bar(stat="identity") +
  # facet_wrap(~markers, scales="free_y", ncol = 5) +
  facet_wrap(~markers, scales="free_y", nrow = 5) +
  scale_x_reordered() +
  coord_flip() +
  theme_bw() +
  theme(strip.text.x = element_text(size = 18),
        text = element_text(size=18),
        axis.text=element_text(size=18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")+
  ylab("Rank of Accuracy (1 - F1)") +
  xlab(NULL)

ggplot(dist.marker[dist.marker$ly=="Scenario 3.2",], aes(x=reorder_within(method, -SD, markers), y=SD, fill=method)) +
  geom_bar(stat="identity") +
  # facet_wrap(~markers,  scales="free_y", ncol = 5) +
  facet_wrap(~markers,  scales="free_y", nrow = 5) +
  scale_x_reordered() +
  coord_flip() +
  theme_bw() +
  theme(strip.text.x = element_text(size = 18),
        text = element_text(size=18),
        axis.text=element_text(size=18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")+
  ylab("Rank of Robustness (SD)") +
  xlab(NULL)







