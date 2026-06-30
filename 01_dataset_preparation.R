setwd("R:/HBMSC2022-A7264")
library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)
library(harmony)
library(SingleCellExperiment)
library(SingleR)
library(celldex)

load("08_PNAS_0.4R_23k/pnas_qc.rda")
load("11_NC_qc_batch/batch1_qc.rda")
load("11_NC_qc_batch/batch2_qc.rda")
load("11_NC_qc_batch/batch3_qc.rda")

pnas_sc <- merge(pnas_qc[[1]], y=c(pnas_qc[[2]], pnas_qc[[3]], pnas_qc[[4]], pnas_qc[[5]], pnas_qc[[6]], pnas_qc[[7]]))
remove(pnas_qc)

universe <- Reduce(intersect, list(rownames(batch1_qc),rownames(batch2_qc), rownames(batch3_qc), rownames(pnas_sc)))
length(universe)

batch1_qc <- batch1_qc[universe, ]
batch2_qc <- batch2_qc[universe, ]
batch3_qc <- batch3_qc[universe, ]
pnas_sc <- pnas_sc[universe, ]

batch1_qc <- subset(batch1_qc, nFeature_RNA < 5000)
batch2_qc <- subset(batch2_qc, nFeature_RNA < 5000)
batch3_qc <- subset(batch3_qc, nFeature_RNA < 5000)
pnas_sc <- subset(pnas_sc, nFeature_RNA < 5000)

pnas_sc_meta <- pnas_sc@meta.data[,c("biosample_id","donor_id","Age", "time_post_partum_weeks", "percent.mt")]
colnames(pnas_sc_meta) <- c("Source", "Individual", "Age", "Infant_age","percent.mt")
pnas_sc_meta$batch="batch0"
pnas_sc_meta <- pnas_sc_meta %>% mutate(Parity = case_when(Individual == 'BM07' | Individual == 'BM09' ~ 3,
                                                           Individual == 'BM19' ~ 2,
                                                           TRUE ~ 1))

pnas_sc_meta <- pnas_sc_meta %>% mutate(Age = case_when(Source %in% c('B2', 'BM5', 'Bfresh', 'BT_t3_lot1', 'BT_t3_old') ~ 30,
                                                        Source == 'BM05_26wk_r1' ~ 31,
                                                        Source %in% c('K1', 'K2', 'Kfresh') ~ 25,
                                                        TRUE ~ as.numeric(Age)))

pnas_sc_meta <- pnas_sc_meta %>% mutate(Infant_age = case_when(Source %in% c('B2', 'BM5') ~ 27.29,
                                                               Source %in% c('BM05_26wk_r1', 'Bfresh')  ~ 33,
                                                               Source %in% c('BT_t3_lot1', 'BT_t3_old')  ~ 46.29,
                                                               Source %in% c('K1', 'K2')  ~ 84.57,
                                                               Source == 'Kfresh' ~ 90.29,
                                                               TRUE ~ as.numeric(Infant_age)))

batch0 <- CreateSeuratObject(counts=pnas_sc[['RNA']], project = "batch0", meta.data = pnas_sc_meta)

remove(pnas_sc_meta, pnas_sc)

batch1_qc@meta.data <- batch1_qc@meta.data %>% mutate(Infant_age = Infant_age*4)
batch2_qc@meta.data <- batch2_qc@meta.data %>% mutate(Infant_age = as.numeric(Infant_age)*4)
batch3_qc@meta.data <- batch3_qc@meta.data %>% mutate(Infant_age = Infant_age*4)

batch0@meta.data$datasets <- "PNAS"
batch1_qc@meta.data$datasets <- "NC"
batch2_qc@meta.data$datasets <- "NC"
batch3_qc@meta.data$datasets<- "NC"
gc()

batch0 <- SCTransform(batch0, vst.flavor = "v2", vars.to.regress = "percent.mt")
batch0 <- RunPCA(batch0)
batch1_qc <- SCTransform(batch1_qc, vst.flavor = "v2", vars.to.regress = "percent.mt")
batch1_qc <- RunPCA(batch1_qc)
batch2_qc <- SCTransform(batch2_qc, vst.flavor = "v2", vars.to.regress = "percent.mt")
batch2_qc <- RunPCA(batch2_qc)
batch3_qc <- SCTransform(batch3_qc, vst.flavor = "v2", vars.to.regress = "percent.mt")
batch3_qc <- RunPCA(batch3_qc)
gc()

combined <- merge(batch0, y=c(batch1_qc, batch2_qc, batch3_qc))
rm(batch0, batch1_qc, batch2_qc, batch3_qc, universe)

VariableFeatures(combined[["SCT"]]) <- rownames(combined[["SCT"]]@scale.data)
combined <- RunPCA(combined, features = VariableFeatures(combined), assay = "SCT")
gc()

combined <- RunHarmony(combined, c("batch", "datasets"), plot_convergence = TRUE, assay.use = "SCT")
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:35, assay = "SCT")
# Warning: The following arguments are not used: k.param
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:35, k.param = 315, assay = "SCT", prune.SNN = 1/15)
combined <- FindClusters(combined, resolution = 0.8)

DefaultAssay(combined) <- "RNA"
sce <- as.SingleCellExperiment(DietSeurat(combined))
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
combined@meta.data$hpca.main <- hpca.main$pruned.labels
DefaultAssay(combined) <- "SCT"
rm(hpca.main, hpca.ref, sce)

combined <- CellCycleScoring(combined, 
                             s.features = cc.genes.updated.2019$s.genes, 
                             g2m.features = cc.genes.updated.2019$g2m.genes, 
                             set.ident = FALSE)

DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.001, raster = FALSE)
DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.001, raster = FALSE, group.by = "Phase")
DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.001, raster = FALSE, group.by = "hpca.main")


# cluster and label in epithelial cells by re-clustering
epi <- subset(combined, idents = c(0,1,2,3,4,6,7,8,9,11,12,13))
epi <- RunPCA(epi, features = VariableFeatures(epi), assay = "SCT")
epi <- RunUMAP(epi, reduction = "harmony", dims = 1:30, assay = "SCT")
epi <- FindNeighbors(epi, reduction = "harmony", dims = 1:30, k.param = 311, assay = "SCT")
epi <- FindClusters(epi, resolution = 0.5)

DimPlot(epi, label = TRUE, group.by = "hpca.main")
DimPlot(epi, label = TRUE, split.by = "hpca.main")

DefaultAssay(epi) <- "RNA"
sce <- as.SingleCellExperiment(DietSeurat(epi))
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
hpca.main <- SingleR(test = sce,assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.main)
epi@meta.data$hpca.main.epi <- hpca.main$pruned.labels
DefaultAssay(epi) <- "SCT"
rm(hpca.main, hpca.ref, sce)
as.data.frame.matrix(table(epi@meta.data$hpca.main.epi, epi@meta.data$hpca.main)) # labeling for combined and epi cells are same

epi@meta.data$Phase_combined <- epi@meta.data$Phase
epi <- CellCycleScoring(epi, 
                        s.features = cc.genes.updated.2019$s.genes, 
                        g2m.features = cc.genes.updated.2019$g2m.genes, 
                        set.ident = FALSE)

epi <- AddModuleScore(epi, features = list("KRT5"), name = "KRT5_level")

DimPlot(epi, label = T)
DimPlot(epi, group.by = "Phase")
DimPlot(epi, group.by = "hpca.main.epi")
DimPlot(epi, group.by = "hpca.main")
DimPlot(epi, group.by = "Phase_combined")
DimPlot(subset(epi, subset=seurat_clusters %in% c(3,6,8)),label = T, group.by = "Phase")

epi.LC1.barcodes <- rownames(epi@meta.data)[epi@meta.data$seurat_clusters %in% c(3,6,8)]
epi.basal.barcodes <- rownames(epi@meta.data)[epi@meta.data$KRT5_level1 >= 0 & epi@meta.data$seurat_clusters == 9]
epi.cellcycle.barcodes <- rownames(epi@meta.data)[epi@meta.data$seurat_clusters==10]
epi.LC2.barcodes <- setdiff(rownames(epi@meta.data), c(epi.LC1.barcodes, epi.basal.barcodes, epi.cellcycle.barcodes))

epi.cycling.barcodes <- rownames(epi@meta.data)[epi@meta.data$seurat_clusters == 10]
int.cycling.barcodes <- rownames(combined@meta.data)[combined@meta.data$seurat_clusters == 13]


ggvenn::ggvenn(list(Integrated=int.cycling.barcodes, Epithelial_subset=epi.cycling.barcodes), 
               fill_color = c("#0073C2FF", "#EFC000FF"),
               stroke_size = 0.5, set_name_size = 4)


epi@meta.data <- epi@meta.data %>%
  mutate(bm_cell_type = case_when(rownames(epi@meta.data) %in% epi.LC1.barcodes ~ "LC1",
                                  rownames(epi@meta.data) %in% epi.LC2.barcodes ~ "LC2",
                                  rownames(epi@meta.data) %in% epi.basal.barcodes ~ "Basal",
                                  rownames(epi@meta.data) %in% epi.cellcycle.barcodes ~ "Cycling",
                                  .default = "Unknown"))

DimPlot(epi, label = TRUE, group.by = "bm_cell_type")

# cluster and label immune cells
imm <- subset(x = combined, idents = c(5, 10, 14))
imm <- RunPCA(imm, features = VariableFeatures(imm), assay = "SCT")
DefaultAssay(imm) <- "SCT"
imm <- RunUMAP(imm, reduction = "harmony", dims = 1:30, assay = "SCT")
imm <- FindNeighbors(imm, reduction = "harmony", dims = 1:30, k.param = 99, assay = "SCT")
imm <- FindClusters(imm, resolution = 0.5)

DefaultAssay(imm) <- "RNA"
sce <- as.SingleCellExperiment(DietSeurat(imm))
monaco.ref <- celldex::MonacoImmuneData()
monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
imm@meta.data$monaco.main <- monaco.main$pruned.labels
imm@meta.data$monaco.fine <- monaco.fine$pruned.labels
DefaultAssay(imm) <- "SCT"
rm(monaco.fine, monaco.main, monaco.ref, sce)

DimPlot(imm, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.001, group.by = "seurat_clusters")
DimPlot(imm, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.001, raster = FALSE, group.by = "hpca.main")
DimPlot(imm, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.001, raster = FALSE, group.by = "monaco.main")

table(imm@meta.data$hpca.main, imm@meta.data$seurat_clusters)
table(imm@meta.data$hpca.main, imm@meta.data$monaco.main) # this combination is better

imm.neutrophil.barcodes <- rownames(imm@meta.data)[imm@meta.data$monaco.main == "Neutrophils" & imm@meta.data$hpca.main == "Neutrophils"]
imm.B.barcodes <- rownames(imm@meta.data)[imm@meta.data$monaco.main == "B cells" & imm@meta.data$hpca.main == "B_cell"]

T.NK.sub <- subset(imm, subset = seurat_clusters == 3)
T.NK.sub <- subset(T.NK.sub, cells = setdiff(colnames(T.NK.sub),c(imm.neutrophil.barcodes, imm.B.barcodes)))
T.NK.sub <- RunPCA(T.NK.sub, features = VariableFeatures(T.NK.sub), assay = "SCT")
DefaultAssay(T.NK.sub) <- "SCT"
ElbowPlot(object = imm, ndims = 50, reduction = "pca")
T.NK.sub <- RunUMAP(T.NK.sub, reduction = "harmony", dims = 1:20, assay = "SCT")
T.NK.sub <- FindNeighbors(T.NK.sub, reduction = "harmony", dims = 1:20, k.param = 40, assay = "SCT")
T.NK.sub <- FindClusters(T.NK.sub, resolution = 0.8)
DimPlot(T.NK.sub, label = T, group.by = "monaco.main")
DimPlot(T.NK.sub, label = T)
table(T.NK.sub@meta.data$hpca.main, T.NK.sub@meta.data$monaco.main)

T.NK.sub@meta.data <- T.NK.sub@meta.data %>%
  mutate(temp_type = case_when(rownames(T.NK.sub@meta.data) %in% imm.nk.barcodes ~ "NK",
                               rownames(T.NK.sub@meta.data) %in% imm.T.barcodes ~ "T_cell",
                               .default = "Unknown"))

T.NK.sub.clean <- subset(T.NK.sub, subset=temp_type != "Unknown")
Idents(T.NK.sub.clean) <- "temp_type"
T.NK.markers <- FindAllMarkers(T.NK.sub.clean, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1, assay = "SCT", recorrect_umi = FALSE)
T.NK.markers <- T.NK.markers %>% filter(p_val_adj < 0.05) %>% arrange(cluster, desc(pct.1))

T_marker <- T.NK.markers[T.NK.markers$cluster == "T_cell",]$gene
NK_marker <- T.NK.markers[T.NK.markers$cluster == "NK",]$gene

intersect(T_marker, NK_marker)

write.csv(T.NK.markers, file = "")


# library(tidymodels)
# library(tidyverse)
# library(workflows)
# library(tune)

imm.nk.barcodes <- rownames(T.NK.sub@meta.data)[T.NK.sub@meta.data$monaco.main == "NK cells" & T.NK.sub@meta.data$hpca.main == "NK_cell"]
imm.T.barcodes <- rownames(T.NK.sub@meta.data)[T.NK.sub@meta.data$monaco.main %in% c("CD4+ T cells", "CD8+ T cells", "T cells") & T.NK.sub@meta.data$hpca.main == "T_cells"]
imm.TBD.barcodes <- setdiff(colnames(T.NK.sub), c(imm.nk.barcodes, imm.T.barcodes))

write.csv(imm.T.barcodes, file = "33_benchmark/01_dataset/ML_classification/T.barcodes.csv")
saveRDS(T.NK.sub, file = "33_benchmark/01_dataset/ML_classification/T_NK_dataset.rds")


M.M.DC.sub <- subset(imm, subset = seurat_clusters %in% c(0,1,2,5))
M.M.DC.sub <- subset(M.M.DC.sub, cells = setdiff(colnames(M.M.DC.sub), c(imm.nk.barcodes,imm.T.barcodes,imm.TBD.barcodes,imm.neutrophil.barcodes, imm.B.barcodes )))
M.M.DC.sub <- RunPCA(M.M.DC.sub, features = VariableFeatures(M.M.DC.sub), assay = "SCT")
DefaultAssay(M.M.DC.sub) <- "SCT"
ElbowPlot(object = M.M.DC.sub, ndims = 50, reduction = "pca")
M.M.DC.sub <- RunUMAP(M.M.DC.sub, reduction = "harmony", dims = 1:25, assay = "SCT")
M.M.DC.sub <- FindNeighbors(M.M.DC.sub, reduction = "harmony", dims = 1:25, k.param = 15, assay = "SCT")
M.M.DC.sub <- FindClusters(M.M.DC.sub, resolution = 0.8)
DimPlot(M.M.DC.sub, label = TRUE)
DimPlot(M.M.DC.sub, label = TRUE, group.by = "hpca.main")
DimPlot(M.M.DC.sub, label = TRUE, group.by = "monaco.main")

table(M.M.DC.sub@meta.data$monaco.main)
table(M.M.DC.sub@meta.data$hpca.main, M.M.DC.sub@meta.data$monaco.main)
table(M.M.DC.sub@meta.data$hpca.main, M.M.DC.sub@meta.data$seurat_clusters)
table(M.M.DC.sub@meta.data$monaco.main, M.M.DC.sub@meta.data$seurat_clusters)




dc.markers <- c("CCR7", "CD1C", "FLT3", "CLEC4C", "LY75", "TCF4")
M.M.DC.sub <- AddModuleScore(M.M.DC.sub, features = list(dc.markers), name = "dc.markers")
FeaturePlot(M.M.DC.sub, features = "dc.markers1")
M.M.DC.sub@meta.data <- M.M.DC.sub@meta.data %>% 
  mutate(DC_refs.1 = case_when(hpca.main == "DC" & monaco.main == "Dendritic cells" ~ "DC",
                              hpca.main != "DC" | monaco.main != "Dendritic cells" ~ "Not_DC"))
DimPlot(M.M.DC.sub, label = TRUE, group.by = "DC_refs.1")
summary(M.M.DC.sub@meta.data$dc.markers1)
plot(density(M.M.DC.sub@meta.data$dc.markers1))
table(M.M.DC.sub@meta.data$dc.markers1 > 0, M.M.DC.sub@meta.data$DC_refs.1)
M.M.DC.sub@meta.data <- M.M.DC.sub@meta.data %>% 
  mutate(DC_refs.2 = case_when(seurat_clusters == 4 & monaco.main == "Dendritic cells" ~ "DC",
                               monaco.main != "Dendritic cells" & monaco.main != "Monocytes" ~ "unknown",
                               .default = "Macrophage"))
DimPlot(M.M.DC.sub, label = TRUE, group.by = "DC_refs.2")

imm.macrophage.barcodes <- colnames(M.M.DC.sub)[M.M.DC.sub@meta.data$DC_refs.2 == "Macrophage"]
imm.dc.barcodes <- colnames(M.M.DC.sub)[M.M.DC.sub@meta.data$DC_refs.2 == "DC"]


write.csv(epi.cellcycle.barcodes, file = "33_benchmark/01_dataset/cell_type_barcodes/epi_cellcycle_barcodes.csv")
write.csv(imm.T.barcodes, file = "33_benchmark/01_dataset/cell_type_barcodes/imm_T_barcodes.csv")

imm@meta.data <- imm@meta.data %>%
  mutate(bm_cell_type = case_when(rownames(imm@meta.data) %in% imm.B.barcodes ~ "B_cell",
                                  rownames(imm@meta.data) %in% imm.dc.barcodes ~ "DC",
                                  rownames(imm@meta.data) %in% imm.macrophage.barcodes ~ "Macrophage/Monocyte",
                                  rownames(imm@meta.data) %in% imm.neutrophil.barcodes ~ "Neutrophil",
                                  rownames(imm@meta.data) %in% imm.nk.barcodes ~ "NK",
                                  rownames(imm@meta.data) %in% imm.T.barcodes ~ "T_cell",
                                  .default = "Unknown"))

DimPlot(imm, label = TRUE, group.by = "bm_cell_type")


combined@meta.data <- combined@meta.data %>%
  mutate(bm_cell_type = case_when(rownames(combined@meta.data) %in% imm.B.barcodes ~ "B_cell",
                                  rownames(combined@meta.data) %in% imm.dc.barcodes ~ "DC",
                                  rownames(combined@meta.data) %in% imm.macrophage.barcodes ~ "Macrophage/Monocyte",
                                  rownames(combined@meta.data) %in% imm.neutrophil.barcodes ~ "Neutrophil",
                                  rownames(combined@meta.data) %in% imm.nk.barcodes ~ "NK",
                                  rownames(combined@meta.data) %in% imm.T.barcodes ~ "T_cell",
                                  rownames(combined@meta.data) %in% epi.basal.barcodes ~ "Basal",
                                  rownames(combined@meta.data) %in% epi.cellcycle.barcodes ~ "Cycling",
                                  rownames(combined@meta.data) %in% epi.LC1.barcodes ~ "LC1",
                                  rownames(combined@meta.data) %in% epi.LC2.barcodes ~ "LC2",
                                  .default = "Unknown"))

combined@meta.data <- combined@meta.data %>%
  mutate(bm_general = case_when(bm_cell_type %in% c("Macrophage/Monocyte", "Neutrophil", "T_cell", "NK", "B_cell", "DC") ~ "IMMUNE",
                                bm_cell_type %in% c("LC2", "LC1", "Cycling", "Basal") ~ "EPITHELIAL",
                                .default = "Unknown"))

combined@meta.data <- combined@meta.data %>%
  mutate(Stage.1 = case_when(Infant_age < 4 ~ "< 1_mon",
                             Infant_age >= 4 & Infant_age < 12 ~ "1~3_mon",
                             Infant_age >= 12 & Infant_age < 24 ~ "3~6_mon",
                             Infant_age >= 24 ~ "> 6_mon"))

combined@meta.data$Stage.1 <- factor(combined@meta.data$Stage.1, levels = c("< 1_mon", "1~3_mon", "3~6_mon","> 6_mon"))
df.temp <- data.frame(table(combined@meta.data$Stage.1, combined@meta.data$bm_general))
df.temp <- df.temp[df.temp$Var2!="Unknown",]
df.temp <- df.temp %>% group_by(Var1) %>% mutate(percent = Freq/sum(Freq))

ggplot(data=df.temp, aes(x=Var1, y=percent, fill=Var2)) +
  geom_bar(stat="identity", position=position_dodge())




combined <- combined[,!is.na(combined@meta.data[,"bm_cell_type"])]
DimPlot(combined, label = TRUE, group.by = "bm_cell_type", raster=FALSE, repel = TRUE)
DimPlot(combined, label = TRUE, group.by = "bm_general", raster=FALSE)

imm.clean <- subset(imm, subset= bm_cell_type != "Unknown")
DimPlot(imm.clean, label = TRUE, group.by = "bm_cell_type", raster=FALSE, repel=TRUE)
DimPlot(epi, label = TRUE, group.by = "bm_cell_type", raster=FALSE, repel=TRUE)
Nebulosa::plot_density(combined, "KRT5")
Nebulosa::plot_density(epi, "KRT5")

saveRDS(combined, file = "33_benchmark/01_dataset/combined.rds")
saveRDS(epi, file = "33_benchmark/01_dataset/epi.rds")
saveRDS(imm, file = "33_benchmark/01_dataset/imm.rds")

M.M.DC.sub@meta.data <- M.M.DC.sub@meta.data %>% 
  mutate(DC_refs.3 = case_when(seurat_clusters == 7 & DC_refs.2 == "Macrophage" ~ "Macrophage.small",
                               seurat_clusters != 7 & DC_refs.2 == "Macrophage" ~ "Macrophage.big",
                               .default = "others"))
DimPlot(M.M.DC.sub, label = TRUE, group.by = "DC_refs.3")
DimPlot(M.M.DC.sub, label = TRUE, group.by = "hpca.main")

Idents(M.M.DC.sub) <- "DC_refs.3"
M.M.DC.sub <- PrepSCTFindMarkers(M.M.DC.sub)
xx <- FindMarkers(M.M.DC.sub, test.use = "wilcox", slot = "data",
                  ident.1 = "Macrophage.small", ident.2 = "Macrophage.big",
                  logfc.threshold = 0.1, min.pct = 0.1, recorrect_umi = TRUE, only.pos = FALSE)
Idents(M.M.DC.sub) <- "seurat_clusters"
xx <-  xx %>% 
  filter(p_val_adj < 0.05) %>%
  filter(pct.1 > 0.5 & pct.2 > 0.5) %>%
  arrange(desc(avg_log2FC))

big.gene <- rownames(xx)[xx$avg_log2FC < 0]
small.gene <- rownames(xx)[xx$avg_log2FC > 0]

write.csv(xx, file="33_benchmark/01_dataset/macrophages/deg_macrophage_clusters_smallVSbig.csv")

# nk.markers <- c("CD300A", "IL2RB", "KIT", "NCAM1", "CCL3", "KLRD1", "GZMB", "XCL1", "XCL2")
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
geneList <- xx$avg_log2FC
names(geneList) <- rownames(xx)

ego <- gseGO(geneList     = geneList,
             OrgDb        = org.Hs.eg.db,
             keyType      = "SYMBOL",
             ont          = "BP",
             pvalueCutoff = 0.05,
             minGSSize    = 10,
             maxGSSize    = 500,
             verbose      = TRUE, 
             nPerm        = 10000,
             by           = "DOSE")

kk <- gseKEGG(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              organism     = 'hsa',
              keyType      = "SYMBOL",
              minGSSize    = 120,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

goplot(ego)
dotplot(ego)
cnetplot(ego)

m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, gene_symbol)
macrophage <- GSEA(geneList = geneList, TERM2GENE = m_t2g)
barplot(macrophage, showCategory = 25)
cnetplot(macrophage)
dotplot(macrophage)

# ITGAX in CD206 neg/low cell is in DEG in big (fetal derived macophage paper)

big.CD_surface_marker <- str_subset(big.gene, "^CD[0-9]", negate = FALSE)
small.CD_surface_marker <- str_subset(small.gene, "^CD[0-9]", negate = FALSE)
