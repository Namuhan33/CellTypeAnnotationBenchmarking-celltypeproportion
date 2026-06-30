imm_cell = 9317
epi_cell = 97252

B_pct = 62/imm_cell
DC_pct = 626/imm_cell
T_pct = 1427/imm_cell
NK_pct = 64/imm_cell
Mac_pct = 6235/imm_cell
Neu_pct = 903/imm_cell

LC1_pct = 19604/epi_cell
LC2_pct = 76652/epi_cell
Basal_pct = 503/epi_cell
Cyc_pct = 493/epi_cell

imm_vec = c(B_pct, DC_pct, T_pct, NK_pct, Mac_pct, Neu_pct)
names(imm_vec) <- c("B_cell", "DC", "T_cell", "NK", "Macrophage", "Neutrophil")

epi_vec = c(LC1_pct, LC2_pct, Basal_pct, Cyc_pct)
names(epi_vec) <- c("LC1", "LC2", "Basal", "G2MS")

# Session 1
# P1 = 10:1 in Q and R
p=10

results = c(1/(1+p)*imm_vec, p/(1+p)*epi_vec) < 0.05

names(results)[results==TRUE] # rare
names(results)[results==FALSE] # abundant

# P2 = 5:1 in Q and R
p=5

results = c(1/(1+p)*imm_vec, p/(1+p)*epi_vec) < 0.05

names(results)[results==TRUE] # rare
names(results)[results==FALSE] # abundant
# P3 = 2:1 in Q and R
p=2

results = c(1/(1+p)*imm_vec, p/(1+p)*epi_vec) < 0.05

names(results)[results==TRUE] # rare
names(results)[results==FALSE] # abundant
# P4 = 1:1 in Q and R
p=1

results = c(1/(1+p)*imm_vec, p/(1+p)*epi_vec) < 0.05

names(results)[results==TRUE] # rare
names(results)[results==FALSE] # abundant
# P5 = 0.5:1 in Q and R
p=0.5

results = c(1/(1+p)*imm_vec, p/(1+p)*epi_vec) < 0.05

names(results)[results==TRUE] # rare
names(results)[results==FALSE] # abundant


# Session 2, 1:1 in Q
# P6 = 10:1 in R and 1:1 in Q
p=1

results = c(1/(1+p)*imm_vec, p/(1+p)*epi_vec) < 0.05

names(results)[results==TRUE] # rare
names(results)[results==FALSE] # abundant

# Inter-dataset, brain
# cell-type1
rare <- c("Astrocyte", "Chandelier", "Endothelial", "L5/6_NP", "L5_ET", "L6_CT", 
          "L6_IT_Car3", "L6b", "Lamp5", "Microglia", "Oligodendrocyte", "OPC", 
          "Pax6", "Sncg", "Sst_Chodl", "VLMC")
abundant <- c("IT", "L4_IT", "Vip", "Pvalb", "Sst")

# cell-type2
rare <- c("Astrocyte", "Chandelier", "Endothelial", "L5/6_NP", 
          "L5_ET", "L6_CT", "L6b", "Lamp5", "Microglia", 
          "Oligodendrocyte", "OPC", "Pax6", "Sncg", "VLMC")
abundant <- c("IT", "Vip", "Sst", "Pvalb")

# cell-type3
rare <- c("Astrocyte", "Endothelial", "Microglia", "Oligodendrocyte", "OPC", "VLMC")
abundant <- c("GABAergic", "Glutamatergic")

# cell-type4
rare <- c("NonNeuronal")
abundant <- c("Neuronal")



# Inter-dataset, PBMC
# cell-type-1
table(covid$cell_type_1)/sum(table(covid$cell_type_1))*100
table(covid$cell_type_2)/sum(table(covid$cell_type_2))*100
table(covid$cell_type_3)/sum(table(covid$cell_type_3))*100

rare <- c("CD16_monocyte", "CD4_T", "CD8_cytotoxic_T", "DC", "Eosinophil", 
          "gd_T", "Granulocyte", "Neutrophil", "Platelet", "RBC")
abundant <- c("B", "CD14_monocyte", "CD8_mem_T", "CD4_mem_T", "CD4_naive_T", "NK")
# cell-type-2
rare <- c("DC", "gd_T", "Granulocyte", "Platelet", "RBC")
abundant <- c("B", "Monocyte", "CD8_T", "CD4_T", "NK")
# cell-type-3
rare <- c("Myeloid")
abundant <- c("Lymphoid")

# Inter-dataset, pancreas
table(seg$cell_type_1)/sum(table(seg$cell_type_1))*100
table(seg$cell_type_2)/sum(table(seg$cell_type_2))*100
table(seg$cell_type_3)/sum(table(seg$cell_type_3))*100

# cell-type 1
rare <- c("Endocrine_coexpression", "Endothelial", "Epsilon",  "Mast", 
          "MHC_class_II", "PSC", "Unclassified_endocrine", "Unclassified_exocrine")
abundant <- c("Delta", "Alpha", "Gamma", "Ductal", "Acinar", "Beta")

# cell-type 2
rare <- c("Endothelial", "Myeloid", "PSC")
abundant <- c("Endocrine", "Ductal", "Exocrine")

# cell-type 3
rare <- c("Non-Endocrine")
abundant <- c("Endocrine")











