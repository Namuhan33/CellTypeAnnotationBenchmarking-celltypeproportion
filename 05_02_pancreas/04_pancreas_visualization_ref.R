setwd("R:/ANNCOM2023-Q6025")
library(Seurat)
library(tidyverse)
library(dplyr)
library(caret)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ggh4x)
library(tibble)
library(ggforce)
library(ggbreak)


method.ref <- c("CHETAH", "scmapcluster", "SingleCellNet", "SingleR", "scClassify", "scPred")

combine_everything_ct1 <- function(methods){
  results.list <- list()
  for(method in methods){
    for (prp in 1:14) {
      path <- paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/03_evaluation/", method, "/P", prp, "/")
      df.1 <- read.csv(paste0(path, "evaluating_result_K5.csv"), row.names = 1)
      df.1 <- df.1[,c(5,9:ncol(df.1))]
      vec.1 <- colMeans(df.1, na.rm = TRUE)
      vec.1 <- append(vec.1, c(method=method, proportion=paste0("Prop", prp), celltype="Cell-type-1"))
      results.list <- rlist::list.append(results.list, vec.1)
    }
  }
  return(results.list)
}
combine_everything_ct2 <- function(methods){
  results.list <- list()
  for(method in methods){
    for (prp in 1:14) {
      path <- paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/03_evaluation/", method, "/P", prp, "/")
      df.1 <- read.csv(paste0(path, "evaluating_result_K5_CT2.csv"), row.names = 1)
      df.1 <- df.1[,c(5,9:ncol(df.1))]
      vec.1 <- colMeans(df.1, na.rm = TRUE)
      vec.1 <- append(vec.1, c(method=method, proportion=paste0("Prop", prp), celltype="Cell-type-2"))
      results.list <- rlist::list.append(results.list, vec.1)
    }
  }
  return(results.list)
}
combine_everything_ct3 <- function(methods){
  results.list <- list()
  for(method in methods){
    for (prp in 1:14) {
      path <- paste0("08_StartOver_NewStory/05_interdatasets/02_pancreas/03_evaluation/", method, "/P", prp, "/")
      df.1 <- read.csv(paste0(path, "evaluating_result_K5_CT3.csv"), row.names = 1)
      df.1 <- df.1[,c(5,9:ncol(df.1))]
      vec.1 <- colMeans(df.1, na.rm = TRUE)
      vec.1 <- append(vec.1, c(method=method, proportion=paste0("Prop", prp), celltype="Cell-type-3"))
      results.list <- rlist::list.append(results.list, vec.1)
    }
  }
  return(results.list)
}

res.1 <- combine_everything_ct1(methods = method.ref)
res.2 <- combine_everything_ct2(methods = method.ref)
res.3 <- combine_everything_ct3(methods = method.ref)

pancreas.result <- do.call("rbind", list(data.frame(t(sapply(res.1, c))), data.frame(t(sapply(res.2, c))),
                                      data.frame(t(sapply(res.3, c)))))

pancreas.result$proportion <- factor(pancreas.result$proportion, levels= paste0("Prop",c(8:14, 1:7)))
pancreas.result <- pancreas.result %>% mutate(ly = case_when(
  proportion %in% paste0("Prop", 1:7) ~ "Scenario 1",
  proportion %in% paste0("Prop", 8:14) ~ "Scenario 2"
))
pancreas.result$celltype <- factor(pancreas.result$celltype, levels = c("Cell-type-3", "Cell-type-2", "Cell-type-1"))
pancreas.result$method[pancreas.result$method=="scmapcluster"] <- "scmap-cluster"

prop.labs <- c("Endo:Non(10:1)[#Q=#R]", "Endo:Non(5:1)[#Q=#R]", "Endo:Non(2:1)[#Q=#R]", "Endo:Non(1:1)[#Q=#R]", "Endo:Non(0.5:1)[#Q=#R]", "Endo:Non(0.2:1)[#Q=#R]", "Endo:Non(0.1:1)[#Q=#R]",
                 "Endo:Non(10:1)[#Q<#R]", "Endo:Non(5:1)[#Q<#R]", "Endo:Non(2:1)[#Q<#R]", "Endo:Non(1:1)[#Q<#R]", "Endo:Non(0.5:1)[#Q<#R]", "Endo:Non(0.2:1)[#Q<#R]", "Endo:Non(0.1:1)[#Q<#R]")
names(prop.labs) <- paste0("Prop", 1:14)

ly.labs <- c("#Q = #R", "#Q < #R")
names(ly.labs) <- c("Scenario 1", "Scenario 2")

celltype.labs <- c("Cell-type Level 1", "Cell-type Level 2", "Cell-type Level 3")
names(celltype.labs) <- c("Cell-type-1", "Cell-type-2", "Cell-type-3")
# F1
metrics <- c("WeightedAveF1")
df <- pancreas.result[,c(metrics, "method", "proportion", "ly", "celltype")]
df[,metrics] <- as.numeric(df[,metrics])



ggplot(data=df, mapping=aes(x=method, y=get(metrics), fill=proportion)) +
  geom_bar(stat = "identity", position=position_dodge(0.8), width=0.8) +
  xlab(label = "Endocrine and non-endocrine cells varying in reference") +
  ylab(label = "Weighted Ave F1 Score") +
  labs(fill= "Proportions") +
  theme_bw()+
  facet_grid( ly ~ celltype , labeller = labeller(ly=ly.labs, celltype=celltype.labs)) +
  theme(strip.text.x = element_text(size = 20),
        text = element_text(size=20),
        axis.text=element_text(size=20),
        axis.text.x = element_text(size=15))+
  geom_text(aes(label = round(get(metrics), 3)), position = position_dodge(.8), size = 4, hjust=1) +
  coord_flip() +
  scale_fill_manual(values=c("#9d92d6", "#66b0f0", "#fffd7c", "#c69189", "#ffb169", "#83d0ff", "#046C9A",
                             "#faa9b4", "#FAD510","#EABE94", "#9986A5", "#CCC591", "#F98400", "#c3abdb"),
                    guide=guide_legend(reverse = TRUE),
                    labels = prop.labs)

# Honesty
metrics <- c("honesty")
df <- pancreas.result[,c(metrics, "method", "proportion", "ly", "celltype")]
df[,metrics] <- as.numeric(df[,metrics])
df <- df[df$celltype %in% c("Cell-type-1"),]

ggplot(data=df, mapping=aes(x=method, y=get(metrics), fill=proportion)) +
  geom_bar(stat = "identity", position=position_dodge(0.8), width=0.8) +
  xlab(label = "Endocrine:non-endocrine cells varying in reference") +
  ylab(label = "Honesty") +
  labs(fill= "Proportions") +
  theme_bw()+
  facet_grid(ly ~ celltype,labeller = labeller(ly=ly.labs, celltype=celltype.labs)) +
  theme(strip.text.x = element_text(size = 20),
        text = element_text(size=20),
        axis.text=element_text(size=20))+
  geom_text(aes(label = round(get(metrics), 2)), position = position_dodge(.8), size = 3, hjust=-0.1) +
  coord_flip() +
  scale_fill_manual(values=c("#9d92d6", "#66b0f0", "#fffd7c", "#c69189", "#ffb169", "#83d0ff", "#046C9A",
                             "#faa9b4", "#FAD510","#EABE94", "#9986A5", "#CCC591", "#F98400", "#c3abdb"),
                    guide=guide_legend(reverse = TRUE),
                    labels = prop.labs)


# Precision-Recall
metrics <- c("Precision", "Recall")
df <- pancreas.result[,c(metrics, "method", "proportion", "ly", "celltype")]
# df <- df[df$ly=="Scenario 1",]
df <- df %>% pivot_longer(cols=c('Precision', 'Recall'), names_to='Precision_Recall', values_to='p_r')
df$p_r <- as.numeric(df$p_r)


ggplot(data=df, mapping=aes(x=method, y=p_r, fill=proportion)) +
  geom_bar(stat = "identity", position=position_dodge(0.8), width=0.8) +
  xlab(label = "Endocrine:non-endocrine cells varying in reference") +
  ylab(label = "Weighted Ave Precision & Recall") +
  labs(fill= "Proportions") +
  theme_bw()+
  facet_grid(Precision_Recall~celltype, labeller = labeller(celltype=celltype.labs)) +
  theme(strip.text.x = element_text(size = 20),
        text = element_text(size=20),
        axis.text=element_text(size=20))+
  geom_text(aes(label = round(p_r, 3)), position = position_dodge(.8), size = 4, hjust=1) +
  coord_flip() +
  scale_fill_manual(values=c("#9d92d6", "#66b0f0", "#fffd7c", "#c69189", "#ffb169", "#83d0ff", "#046C9A",
                             "#faa9b4", "#FAD510","#EABE94", "#9986A5", "#CCC591", "#F98400", "#c3abdb"),
                    guide=guide_legend(reverse = TRUE),
                    labels =prop.labs)



# rare cell type
metrics <- c("WeightedAveF1.rare", "WeightedAveF1.rich")

df <- pancreas.result[,c(metrics, "method", "proportion", "ly", "celltype")]
df$WeightedAveF1.rare <- as.numeric(df$WeightedAveF1.rare)
df$WeightedAveF1.rich <- as.numeric(df$WeightedAveF1.rich)
df <- df[df$ly=="Scenario 1",]

df <- df %>% pivot_longer(cols=c('WeightedAveF1.rare', 'WeightedAveF1.rich'), names_to='Abundance', values_to='F1')

level.labs <- c("Abundance < 5%", "Abundance > 5%")
names(level.labs) <- c("WeightedAveF1.rare", "WeightedAveF1.rich")

ggplot(data=df, mapping=aes(x=method, y=F1, fill= proportion)) +
  geom_bar(stat = "identity", position=position_dodge(0.8), width=0.8) +
  xlab(label = "") +
  ylab(label = "Weighted Ave F1 Score") +
  labs(fill= "Proportions") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 20),
        text = element_text(size=20),
        axis.text=element_text(size=20),
        axis.text.x=element_text(size=12))+
  geom_text(aes(label = round(F1, 3), angle = 0), position = position_dodge(.8), size = 4, hjust=0.5) +
  coord_flip() +
  facet_grid(Abundance~celltype, scales = "free_x", space = "free_x", labeller = labeller(Abundance = level.labs, celltype=celltype.labs)) +
  scale_fill_manual(values=c("#faa9b4", "#FAD510","#EABE94", "#9986A5", "#CCC591", "#F98400", "#c3abdb"),
                    guide=guide_legend(reverse = TRUE),
                    labels = prop.labs)

metrics <- c("WeightedAveF1.rare", "WeightedAveF1.rich")

df <- pancreas.result[,c(metrics, "method", "proportion", "ly", "celltype")]
df$WeightedAveF1.rare <- as.numeric(df$WeightedAveF1.rare)
df$WeightedAveF1.rich <- as.numeric(df$WeightedAveF1.rich)
df <- df[df$ly=="Scenario 2",]

df <- df %>% pivot_longer(cols=c('WeightedAveF1.rare', 'WeightedAveF1.rich'), names_to='Abundance', values_to='F1')

level.labs <- c("Abundance < 5%", "Abundance > 5%")
names(level.labs) <- c("WeightedAveF1.rare", "WeightedAveF1.rich")

ggplot(data=df, mapping=aes(x=method, y=F1, fill= proportion)) +
  geom_bar(stat = "identity", position=position_dodge(0.8), width=0.8) +
  xlab(label = "") +
  ylab(label = "Weighted Ave F1 Score") +
  labs(fill= "Proportions") +
  theme_bw()+
  theme(strip.text.x = element_text(size = 20),
        text = element_text(size=20),
        axis.text=element_text(size=20),
        axis.text.x=element_text(size=12))+
  geom_text(aes(label = round(F1, 3), angle = 0), position = position_dodge(.8), size = 4, hjust=0.5) +
  coord_flip() +
  facet_grid(Abundance~celltype, scales = "free_x", space = "free_x", labeller = labeller(Abundance = level.labs, celltype=celltype.labs)) +
  scale_fill_manual(values=c("#9d92d6", "#66b0f0", "#fffd7c", "#c69189", "#ffb169", "#83d0ff", "#046C9A"),
                    guide=guide_legend(reverse = TRUE),
                    labels = prop.labs)

# assignment rate
assignment<- pancreas.result[,c("assignedRate", "method", "proportion", "ly", "celltype")]
colnames(assignment)[1] <- "Rate"
assignment$Rate <- as.numeric(assignment$Rate)
assignment <- assignment[assignment$celltype == "Cell-type-1",]

ggplot(assignment, aes(y=Rate, x=proportion)) +
  geom_bar(stat="identity", fill = "#EAD3BF") +
  facet_wrap(~ method, nrow=3) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 20),
        text = element_text(size=20),
        legend.position = "none",
        axis.text=element_text(size=12),
        axis.ticks.y=element_blank(),
        panel.spacing = unit(0.01, "cm")) +
  ylab("Assigned rate (%)") +
  xlab("") +
  scale_x_discrete(labels=prop.labs,
                   guide = guide_axis(angle = 90),
                   limits=factor(paste0("Prop", 1:14), levels = paste0("Prop", 1:14)))

# robustness: scatter plot
metrics <- c("WeightedAveF1")

df <- pancreas.result[,c(metrics, "method", "proportion", "celltype", "ly")]
df[,metrics] <- as.numeric(df[,metrics])

colnames(df)[2] <- "Method"
colnames(df)[5] <- "Scenario"

df$Scenario[df$Scenario=="Session 1"] <- "#Q = #R"
df$Scenario[df$Scenario=="Session 2"] <- "#Q > #R"

df <- df[df$proportion %in% paste0("Prop", 1:14),]
dist.ref <- df %>%
  dplyr::group_by(Method, celltype, Scenario) %>%
  dplyr::summarise_at(vars(-c(proportion)), list(mean = mean, sd = sd))

ggplot(dist.ref, aes(x=mean, y=sd)) +
  labs(x= "Mean of F1-score", y="Standard Deviation of F1-score")+
  geom_point(aes(colour=Method, shape=Scenario), size=5) +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 20)) +
  facet_wrap(celltype~., labeller = labeller(celltype = celltype.labs)) +
  geom_text_repel(aes(label = Method), size=7)


# runtime
Time <- pancreas.result[,c("Time", "method", "proportion", "celltype", "ly")]
Time <- Time[Time$celltype=="Cell-type-1",]
Time$Time <- as.numeric(Time$Time)

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
  scale_x_discrete(labels = prop.labs,
                   guide = guide_axis(angle = 45))





