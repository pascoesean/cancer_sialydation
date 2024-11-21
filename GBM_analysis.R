# load packages

library(tidyverse)
library(Seurat)

# load data

int_data <- readRDS("data/GBM.RNA.integrated.24.rds")

var_features <- int_data@assays$integrated@var.features |>
  as_tibble()


DimPlot(int_data, reduction = "umap")

VlnPlot(int_data, features = c("ST3GAL4", "CD163", "CD44", "SIGLEC9"))


VlnPlot(int_data, features = c("CD44", "SIGLEC9"))

FeaturePlot(int_data, features = c("ST3GAL4", "CD163", "CD44", "SIGLEC9"))

FeaturePlot(int_data, features = c("CD44","NANOG"), blend = TRUE)
FeaturePlot(int_data, features = c("CD44","S100A4"), blend = TRUE)

FeaturePlot(int_data, features = c("CD44", "ST3GAL6"), blend = TRUE, split.by = 'treatment_1')

# maybe some meaningful coexpression for st3gal5?

# using which cells to look at those
sialcd44 <- WhichCells(int_data, expression = CD44 > 1 &  ("ST3GAL3" > 1 | "ST3GAL6" > 1 |
                                                               "ST3GAL4" > 1))

sig9tams <- WhichCells(int_data, expression = CD163 > 1 &  SIGLEC9 > 1)
sig9s <- WhichCells(int_data, expression = SIGLEC9 > 1)

# lets literally just look at numbers of these in each condition
total_cells <- FetchData(int_data, vars = c("treatment_3", "Pt_number")) |>
  group_by(treatment_3, Pt_number) |>
  summarize(n_total = n())

sialcd44_anno <- FetchData(int_data, vars = c("treatment_3", "Pt_number"), cells = sialcd44)|>
  group_by(treatment_3, Pt_number) |>
  summarize(n_sialcd44s = n())

sig9s_anno <- FetchData(int_data, vars = c("treatment_3", "Pt_number"), cells = sig9s)|>
  group_by(treatment_3, Pt_number) |>
  summarize(n_sig9s = n())

sig9tams_anno <- FetchData(int_data, vars = c("treatment_3", "Pt_number"), cells = sig9tams)|>
  group_by(treatment_3, Pt_number) |>
  summarize(n_sig9tams = n())

# need ntotals by patient number 
totals <- total_cells |>
  group_by(treatment_3) |>
  summarize(total_bycat = sum(n_total))

total_cells_joined <- total_cells |>
  left_join(y = sialcd44_anno) |>
  left_join(y = sig9s_anno) |>
  left_join(y = sig9tams_anno) |>
  mutate(total_forcat = case_when(
    treatment_3 == "Rec" ~ 13813,
    treatment_3 == "nonresponder" ~ 36816,
    treatment_3 == "responder" ~ 67875,
    treatment_3 == "Untreated" ~ 30544
  )) |>
  mutate(prop_sialcd44 = n_sialcd44s/total_forcat,
         prop_sig9s = n_sig9s/total_forcat,
         prop_sig9tams = n_sig9tams/total_forcat)

write_rds(total_cells_joined, file = "data/gbm_total_cells_1.0threshold.rds")

# interesting....
total_cells_joined |>
  select(treatment_3, Pt_number, starts_with("prop_")) |>
  pivot_longer(cols = starts_with("prop_"), names_to = "subset", values_to = "proportion", names_prefix = "prop_") |>
  ggplot(aes(x = treatment_3, y = proportion, fill = Pt_number)) +
  geom_col() +
  facet_wrap(~subset, scales = "free_y") +
  #scale_fill_brewer(palette = "Set3") +
  ggpubr::theme_pubr() +
  labs(y = "proportion of total cells in this category") + 
  theme(axis.text.x = element_text(angle = 45, hjust=  1, vjust = 1))

total_cells_joined |>
  ggplot(aes(x = prop_sialcd44, y = prop_sig9s, color = treatment_3, shape = treatment_3)) +
  geom_point(size = 3) + 
  ggpubr::theme_pubr() 
    
DimPlot(int_data, reduction = "umap", split.by = "treatment_3", cells = sig9tams)

dpos_subset <- subset(x = int_data, subset = SIGLEC9 > 1 & CD44 > 1)

FeaturePlot(int_data, features = c("CD44", "SIGLEC9"), split.by = 'treatment_1', blend = TRUE, cells = dpos)

VlnPlot(dpos_subset, features = c("CD44", "SIGLEC9"), split.by = 'treatment_1')

sigetal_expression <- FetchData(int_data, vars = c("CD44", "ST3GAL4", "ST3GAL3", "ST3GAL6", "ST6GALNAC3", "SIGLEC9", 
                                                   "anno_ident", "treatment_1", "PD1"))

write_rds(sigetal_expression, file = "data/sigetal_expression.rds")

sexpression <- sigetal_expression |>
  mutate(cd44_sig = CD44*SIGLEC9,
         st3_sum = ST3GAL4 + ST3GAL3 + ST3GAL6,
         cd44_stgal = CD44*st3_sum,
         category = case_when(
           treatment_1 == "TMZ" ~ 'tmz_treat',
           treatment_1 == "untreated" ~ 'untreated',
           treatment_1 == "PD1" ~ PD1
         ))

sexpression |>
  #filter(anno_ident %in% c("cDCs", "Microglial", "Monocytes", "Macrophages")) |>
  ggplot(aes(x = SIGLEC9, y = st3_sum)) +
  geom_point(alpha = 0.01) +
  facet_wrap(~anno_ident) +
  labs(y = "summed expression of ST3GAL4 + ST3GAL3 + ST3GAL6") +
  ggpubr::theme_pubr()

sexpression |>
  filter(anno_ident %in% c("cDCs", "Microglial", "Monocytes", "Macrophages")) |>
  ggplot(aes(x = CD44, y = st3_sum)) +
  geom_point(alpha = 0.01) +
  #ggpmisc::stat_poly_line(se = F) +
  #ggpmisc::stat_poly_eq() +
  facet_grid(rows = vars(anno_ident), cols = vars(category)) +
  labs(y = "summed expression of ST3GAL4 + ST3GAL5 + ST3GAL6") +
  ggpubr::theme_pubr()


sexpression |>
  filter(anno_ident %in% c("cDCs", "Microglial", "Monocytes", "Macrophages")) |>
  ggplot(aes(x = log10(cd44_stgal), color = category)) +
  #facet_grid(rows = vars(anno_ident), cols = vars(treatment_1)) +
  geom_freqpoly(stat ="density") +
  facet_wrap(~anno_ident) +
  ggpubr::theme_pubr()


### JUST LOOKING AT MONOCYTES + MICROGLIA

# mono_subset <- subset(x = int_data, subset = anno_ident %in% c("cDCs", "Microglial", "Monocytes", "Macrophages"))

# saveRDS(mono_subset, file = "data/mono_subset.rds")
mono_subset <- readRDS("data/mono_subset.rds")

mono.genes <- rownames(mono_subset)
mono_subset <- ScaleData(mono_subset, features = mono.genes)

mono_subset <- FindVariableFeatures(object = mono_subset)
mono_subset <- RunPCA(mono_subset, features = VariableFeatures(object = mono_subset))
ElbowPlot(mono_subset, ndims = 50)

DimHeatmap(mono_subset, dims = 20:40, cells = 500, balanced = TRUE)

mono_subset <- FindNeighbors(mono_subset, dims = 1:40)
mono_subset <- FindClusters(mono_subset, resolution = 0.05)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
ms.markers <- FindAllMarkers(mono_subset, only.pos = TRUE)
saveRDS(ms.markers, file = "data/mono_subset_markers.rds")

ms.markers %>%
  group_by(cluster) %>%
  #dplyr::filter(avg_log2FC > 1) |>
  View()

mono_subset <- RunUMAP(mono_subset, dims = 1:40)
DimPlot(mono_subset, reduction = "umap", group.by = "Pt_number", split.by = "treatment_3")

VlnPlot(mono_subset, features = c("SIGLEC9", "CD44", "MARCO", "ST3GAL4", "SEPP1"), split.by = 'treatment_1')
FeaturePlot(mono_subset, features = c("SIGLEC9", "CD44", "ST3GAL4", "ST6GALNAC5"))

VlnPlot(mono_subset, features = c("ST6GALNAC5", "ST6GALNAC3", "ST3GAL4-AS1", "ST3GAL5-AS1"))

cluster0.markers <- FindMarkers(mono_subset, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

write.csv(cluster0.markers, file = "data/clust0_markers.csv")
cluster1.markers <- FindMarkers(mono_subset, ident.1 = 1, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
write.csv(cluster1.markers, file = "data/clust1_markers.csv")

mono_subset <- readRDS(file = "data/mono_subset_scaled.rds")

top10 <- ms.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()

DoHeatmap(mono_subset, features = top10$gene) + NoLegend()


#saveRDS(mono_subset, file = "data/mono_subset_scaled.rds")

